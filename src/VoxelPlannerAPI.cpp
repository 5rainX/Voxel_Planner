#include "VoxelPlannerAPI.h"

#include "Module2_Morphology/PoseMaskPopulation.h"
#include "Module2_Morphology/VoxelMorphologyEngine.h"
#include "Module2_Morphology/VoxelIO.h"
#include "Module3_AStar/CoarseAStar.h"
#include "Module3_AStar/JumpAStar.h"
#include "Module2_Morphology/PoseGenerator.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <queue>
#include <new>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>

namespace {

using ProfileClock = std::chrono::high_resolution_clock;

double profileMilliseconds(ProfileClock::duration duration) {
    return std::chrono::duration<double, std::milli>(duration).count();
}

long long profileNanoseconds(ProfileClock::duration duration) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
        duration).count();
}

void logProfileMessage(const std::string& message) {
    std::cout << "[VoxelPlannerAPI::findPaths] " << message
              << '\n' << std::flush;
}

void logProfileStageBegin(const std::string& stage) {
    logProfileMessage("BEGIN " + stage);
}

void logProfileStageEnd(
    const std::string& stage,
    ProfileClock::duration duration) {
    std::cout << "[VoxelPlannerAPI::findPaths] END " << stage
              << ": " << profileMilliseconds(duration) << " ms, "
              << profileNanoseconds(duration) << " ns\n"
              << std::flush;
}

bool hasPositiveConfiguration(
    const voxel_planner::PlannerConfig& config) {
    return std::isfinite(config.busbar_width) &&
           config.busbar_width > 0.0F &&
           std::isfinite(config.busbar_thickness) &&
           config.busbar_thickness > 0.0F &&
           std::isfinite(config.flat_bend_factor) &&
           config.flat_bend_factor > 0.0F &&
           std::isfinite(config.vertical_bend_factor) &&
           config.vertical_bend_factor > 0.0F &&
           config.angle_step_deg > 0 &&
           360 % config.angle_step_deg == 0 &&
           std::isfinite(config.path_blocking_radius) &&
           config.path_blocking_radius > 0.0F &&
           std::isfinite(config.max_overlap_ratio) &&
           config.max_overlap_ratio >= 0.0F &&
           config.max_overlap_ratio <= 1.0F;
}

bool isUnitVector(const voxel_planner::Vector3D& vector) {
    if (!std::isfinite(vector.x) ||
        !std::isfinite(vector.y) ||
        !std::isfinite(vector.z)) {
        return false;
    }
    const double length = std::sqrt(
        vector.x * vector.x +
        vector.y * vector.y +
        vector.z * vector.z);
    return std::abs(length - 1.0) <= 1e-6;
}

void validateEndpointConstraint(
    const voxel_planner::EndpointConstraint& constraint) {
    if (!isUnitVector(constraint.start_tangent) ||
        !isUnitVector(constraint.end_tangent) ||
        !std::isfinite(constraint.min_begin_length) ||
        !std::isfinite(constraint.min_end_length) ||
        constraint.min_begin_length < 0.0F ||
        constraint.min_end_length < 0.0F) {
        throw std::invalid_argument(
            "Endpoint tangents must be unit vectors and terminal lengths "
            "must be finite and non-negative.");
    }
}

int shiftedCoordinate(
    int coordinate,
    double tangent,
    float length,
    int direction) {
    const double roundedOffset = std::round(tangent * length);
    const double shifted = coordinate + direction * roundedOffset;
    if (shifted < std::numeric_limits<int>::min() ||
        shifted > std::numeric_limits<int>::max()) {
        throw std::invalid_argument("Shifted endpoint coordinate overflowed.");
    }
    return static_cast<int>(shifted);
}

voxel_planner::Point3D shiftedEndpoint(
    const voxel_planner::Point3D& point,
    const voxel_planner::Vector3D& tangent,
    float length,
    int direction) {
    return {
        shiftedCoordinate(point.x, tangent.x, length, direction),
        shiftedCoordinate(point.y, tangent.y, length, direction),
        shiftedCoordinate(point.z, tangent.z, length, direction)};
}

double pointDistance(
    const voxel_planner::Point3D& lhs,
    const voxel_planner::Point3D& rhs) {
    const double dx = static_cast<double>(rhs.x) - lhs.x;
    const double dy = static_cast<double>(rhs.y) - lhs.y;
    const double dz = static_cast<double>(rhs.z) - lhs.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

int stepToward(int current, int target) noexcept {
    return current < target ? 1 : (current > target ? -1 : 0);
}

voxel_planner::Point3D pointFromVoxelIndex(
    const VoxelGrid& grid,
    std::size_t voxelIndex) {
    const std::size_t plane =
        static_cast<std::size_t>(grid.width()) * grid.height();
    const int z = static_cast<int>(voxelIndex / plane);
    const std::size_t remainder = voxelIndex % plane;
    const int y = static_cast<int>(remainder / grid.width());
    const int x = static_cast<int>(remainder % grid.width());
    return {x, y, z};
}

void appendContinuous26(
    std::vector<voxel_planner::Point3D>& path,
    const voxel_planner::Point3D& target) {
    if (path.empty()) {
        path.push_back(target);
        return;
    }
    while (!(path.back() == target)) {
        voxel_planner::Point3D next = path.back();
        next.x += stepToward(next.x, target.x);
        next.y += stepToward(next.y, target.y);
        next.z += stepToward(next.z, target.z);
        path.push_back(next);
    }
}

voxel_planner::Point3D addOffset(
    const voxel_planner::Point3D& point,
    const voxel_planner::Point3D& offset) {
    return {point.x + offset.x, point.y + offset.y, point.z + offset.z};
}

constexpr int kHintProjectionMaxRadius = 12;
constexpr double kMinimumHintProjectionSuccessRatio = 0.25;

struct HintProjectionSample {
    voxel_planner::Point3D shifted{};
    voxel_planner::Point3D projected{};
    bool success = false;
};

struct HintProjectionResult {
    std::vector<voxel_planner::Point3D> points;
    std::size_t sourcePointCount = 0U;
    std::size_t projectedPointCount = 0U;
    std::size_t filledPointCount = 0U;
    double successRate = 0.0;
};

bool isGridObstacle(
    const VoxelGrid& grid,
    const voxel_planner::Point3D& point) {
    if (!grid.isValid(point.x, point.y, point.z)) {
        return true;
    }
    return grid.isObstacle(
        grid.index(
            static_cast<std::uint32_t>(point.x),
            static_cast<std::uint32_t>(point.y),
            static_cast<std::uint32_t>(point.z)));
}

bool tryProjectHintPointToFreeSpace(
    const VoxelGrid& grid,
    const voxel_planner::Point3D& point,
    voxel_planner::Point3D& projected) {
    if (!isGridObstacle(grid, point)) {
        projected = point;
        return true;
    }

    bool found = false;
    int bestDistanceSquared = std::numeric_limits<int>::max();
    voxel_planner::Point3D best{};
    for (int radius = 1; radius <= kHintProjectionMaxRadius; ++radius) {
        for (int dz = -radius; dz <= radius; ++dz) {
            for (int dy = -radius; dy <= radius; ++dy) {
                for (int dx = -radius; dx <= radius; ++dx) {
                    if (std::max({
                            std::abs(dx),
                            std::abs(dy),
                            std::abs(dz)}) != radius) {
                        continue;
                    }
                    const voxel_planner::Point3D candidate{
                        point.x + dx,
                        point.y + dy,
                        point.z + dz};
                    if (isGridObstacle(grid, candidate)) {
                        continue;
                    }
                    const int distanceSquared =
                        dx * dx + dy * dy + dz * dz;
                    if (!found ||
                        distanceSquared < bestDistanceSquared ||
                        (distanceSquared == bestDistanceSquared &&
                         std::tie(
                             candidate.z,
                             candidate.y,
                             candidate.x) <
                         std::tie(best.z, best.y, best.x))) {
                        found = true;
                        bestDistanceSquared = distanceSquared;
                        best = candidate;
                    }
                }
            }
        }
        if (found) {
            projected = best;
            return true;
        }
    }
    return false;
}

voxel_planner::Point3D interpolateHintPoint(
    const voxel_planner::Point3D& from,
    const voxel_planner::Point3D& to,
    std::size_t numerator,
    std::size_t denominator) {
    if (denominator == 0U) {
        return from;
    }
    const double ratio =
        static_cast<double>(numerator) / static_cast<double>(denominator);
    return {
        static_cast<int>(std::lround(
            static_cast<double>(from.x) +
            static_cast<double>(to.x - from.x) * ratio)),
        static_cast<int>(std::lround(
            static_cast<double>(from.y) +
            static_cast<double>(to.y - from.y) * ratio)),
        static_cast<int>(std::lround(
            static_cast<double>(from.z) +
            static_cast<double>(to.z - from.z) * ratio))};
}

void appendHintPointDeduplicated(
    std::vector<voxel_planner::Point3D>& points,
    const voxel_planner::Point3D& point) {
    if (points.empty() || !(points.back() == point)) {
        points.push_back(point);
    }
}

bool appendProjectedHintSegment(
    const VoxelGrid& grid,
    const voxel_planner::Point3D& target,
    std::vector<voxel_planner::Point3D>& points) {
    if (points.empty()) {
        appendHintPointDeduplicated(points, target);
        return true;
    }

    voxel_planner::Point3D ideal = points.back();
    while (!(ideal == target)) {
        ideal.x += stepToward(ideal.x, target.x);
        ideal.y += stepToward(ideal.y, target.y);
        ideal.z += stepToward(ideal.z, target.z);

        voxel_planner::Point3D projected{};
        if (!tryProjectHintPointToFreeSpace(grid, ideal, projected)) {
            // The target was already projected successfully. Clamp the
            // unresolved gap to that legal endpoint rather than emitting a
            // failed or obstructed intermediate voxel.
            appendHintPointDeduplicated(points, target);
            return false;
        }
        appendHintPointDeduplicated(points, projected);
    }
    return true;
}

HintProjectionResult projectTopologyHintToExpandedGrid(
    const VoxelGrid& grid,
    const std::vector<voxel_planner::Point3D>& sourcePoints,
    const voxel_planner::Point3D& morphologyOffset) {
    HintProjectionResult result;
    result.sourcePointCount = sourcePoints.size();
    if (sourcePoints.empty()) {
        return result;
    }

    std::vector<HintProjectionSample> samples(sourcePoints.size());
    std::size_t successCount = 0U;
    for (std::size_t index = 0U; index < sourcePoints.size(); ++index) {
        HintProjectionSample& sample = samples[index];
        sample.shifted = addOffset(sourcePoints[index], morphologyOffset);
        sample.success = tryProjectHintPointToFreeSpace(
            grid,
            sample.shifted,
            sample.projected);
        if (sample.success) {
            ++successCount;
            if (!(sample.projected == sample.shifted)) {
                ++result.projectedPointCount;
            }
        }
    }

    result.successRate =
        static_cast<double>(successCount) /
        static_cast<double>(sourcePoints.size());
    if (result.successRate <= kMinimumHintProjectionSuccessRatio) {
        return result;
    }

    std::vector<int> previousSuccess(samples.size(), -1);
    std::vector<int> nextSuccess(samples.size(), -1);
    int lastSuccess = -1;
    for (std::size_t index = 0U; index < samples.size(); ++index) {
        if (samples[index].success) {
            lastSuccess = static_cast<int>(index);
        }
        previousSuccess[index] = lastSuccess;
    }
    int upcomingSuccess = -1;
    for (std::size_t reverseIndex = samples.size();
         reverseIndex > 0U;
         --reverseIndex) {
        const std::size_t index = reverseIndex - 1U;
        if (samples[index].success) {
            upcomingSuccess = static_cast<int>(index);
        }
        nextSuccess[index] = upcomingSuccess;
    }

    for (std::size_t index = 0U; index < samples.size(); ++index) {
        voxel_planner::Point3D projected{};
        if (samples[index].success) {
            projected = samples[index].projected;
        } else {
            const int previous = previousSuccess[index];
            const int next = nextSuccess[index];
            if (previous >= 0 && next >= 0 && previous != next) {
                const std::size_t numerator =
                    index - static_cast<std::size_t>(previous);
                const std::size_t denominator =
                    static_cast<std::size_t>(next - previous);
                const voxel_planner::Point3D interpolated =
                    interpolateHintPoint(
                        samples[static_cast<std::size_t>(previous)].projected,
                        samples[static_cast<std::size_t>(next)].projected,
                        numerator,
                        denominator);
                if (!tryProjectHintPointToFreeSpace(
                        grid,
                        interpolated,
                        projected)) {
                    // If the rounded interpolation still lands in a solid
                    // region, clamp to the nearest already-valid endpoint.
                    projected = numerator * 2U <= denominator
                        ? samples[static_cast<std::size_t>(previous)].projected
                        : samples[static_cast<std::size_t>(next)].projected;
                }
            } else if (previous >= 0) {
                projected =
                    samples[static_cast<std::size_t>(previous)].projected;
            } else if (next >= 0) {
                projected = samples[static_cast<std::size_t>(next)].projected;
            } else {
                return result;
            }
            ++result.filledPointCount;
        }

        if (!appendProjectedHintSegment(grid, projected, result.points)) {
            appendHintPointDeduplicated(result.points, projected);
        }
    }
    return result;
}

void translatePathToPublicCoordinates(
    std::vector<voxel_planner::Point3D>& path,
    const voxel_planner::Point3D& offset) {
    for (voxel_planner::Point3D& point : path) {
        point.x -= offset.x;
        point.y -= offset.y;
        point.z -= offset.z;
    }
}

bool sameVoxelPath(
    const std::vector<Point3D>& lhs,
    const std::vector<Point3D>& rhs) {
    return lhs.size() == rhs.size() &&
        std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

bool containsVoxelPath(
    const std::vector<::PathResult>& paths,
    const std::vector<Point3D>& candidate) {
    for (const ::PathResult& path : paths) {
        if (sameVoxelPath(path.path, candidate)) {
            return true;
        }
    }
    return false;
}

void deduplicatePlanningPaths(
    ::PlanningResult& result,
    std::size_t maximumCount) {
    std::vector<::PathResult> uniquePaths;
    uniquePaths.reserve(std::min(
        result.paths.size(),
        maximumCount));
    for (::PathResult& candidate : result.paths) {
        if (candidate.path.empty() ||
            containsVoxelPath(uniquePaths, candidate.path)) {
            continue;
        }
        uniquePaths.push_back(std::move(candidate));
        if (uniquePaths.size() >= maximumCount) {
            break;
        }
    }
    result.paths = std::move(uniquePaths);
}

void mergeFallbackPlanningPaths(
    ::PlanningResult& destination,
    ::PlanningResult&& fallback,
    std::size_t maximumCount) {
    for (::PathResult& candidate : fallback.paths) {
        if (destination.paths.size() >= maximumCount) {
            break;
        }
        if (candidate.path.empty() ||
            containsVoxelPath(destination.paths, candidate.path)) {
            continue;
        }
        destination.paths.push_back(std::move(candidate));
    }
}

struct RawEndpointResolution {
    voxel_planner::Point3D point{};
    bool shifted = false;
};

RawEndpointResolution resolveRawEndpointOutsideObstacle(
    const VoxelGrid& grid,
    const voxel_planner::Point3D& rawPoint,
    const char* endpointName) {
    if (!grid.isValid(rawPoint.x, rawPoint.y, rawPoint.z)) {
        throw std::invalid_argument(
            std::string("Raw ") + endpointName +
            " point is outside the voxel map.");
    }

    const std::size_t startIndex = grid.index(
        static_cast<std::uint32_t>(rawPoint.x),
        static_cast<std::uint32_t>(rawPoint.y),
        static_cast<std::uint32_t>(rawPoint.z));
    if (!grid.isRawObstacle(startIndex)) {
        return {rawPoint, false};
    }

    std::vector<bool> visited(grid.voxelCount(), false);
    std::queue<std::size_t> frontier;
    visited[startIndex] = true;
    frontier.push(startIndex);

    while (!frontier.empty()) {
        const std::size_t levelCount = frontier.size();
        for (std::size_t i = 0U; i < levelCount; ++i) {
            const std::size_t currentIndex = frontier.front();
            frontier.pop();
            const voxel_planner::Point3D current =
                pointFromVoxelIndex(grid, currentIndex);

            for (int dz = -1; dz <= 1; ++dz) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        if (dx == 0 && dy == 0 && dz == 0) {
                            continue;
                        }
                        const std::int64_t nextX =
                            static_cast<std::int64_t>(current.x) + dx;
                        const std::int64_t nextY =
                            static_cast<std::int64_t>(current.y) + dy;
                        const std::int64_t nextZ =
                            static_cast<std::int64_t>(current.z) + dz;
                        if (nextX < 0 || nextY < 0 || nextZ < 0 ||
                            nextX >= static_cast<std::int64_t>(grid.width()) ||
                            nextY >= static_cast<std::int64_t>(grid.height()) ||
                            nextZ >= static_cast<std::int64_t>(grid.depth())) {
                            continue;
                        }
                        const std::size_t nextIndex = grid.index(
                            static_cast<std::uint32_t>(nextX),
                            static_cast<std::uint32_t>(nextY),
                            static_cast<std::uint32_t>(nextZ));
                        if (visited[nextIndex]) {
                            continue;
                        }
                        visited[nextIndex] = true;
                        if (!grid.isRawObstacle(nextIndex)) {
                            return {
                                pointFromVoxelIndex(grid, nextIndex),
                                true};
                        }
                        frontier.push(nextIndex);
                    }
                }
            }
        }
    }

    throw std::invalid_argument(
        std::string("No raw obstacle-free voxel found near ") +
        endpointName + " point.");
}

voxel_planner::Vector3D unitPoseNormal(
    const voxel_planner::BusbarPose& pose) {
    const double length = std::sqrt(
        static_cast<double>(pose.nx) * pose.nx +
        static_cast<double>(pose.ny) * pose.ny +
        static_cast<double>(pose.nz) * pose.nz);
    return {
        pose.nx / length,
        pose.ny / length,
        pose.nz / length};
}

voxel_planner::Vector3D unitPoseTangent(
    const voxel_planner::BusbarPose& pose) {
    const double length = std::sqrt(static_cast<double>(
        pose.tx * pose.tx + pose.ty * pose.ty + pose.tz * pose.tz));
    return {
        pose.tx / length,
        pose.ty / length,
        pose.tz / length};
}

void appendPublicConditionalPose(
    const VoxelGrid& grid,
    const voxel_planner::Point3D& point,
    std::size_t waypointIndex,
    const voxel_planner::Vector3D& normal,
    const voxel_planner::Vector3D& tangent,
    std::vector<bool>& described,
    std::vector<voxel_planner::PoseDescription>& descriptions) {
    if (waypointIndex >= described.size() ||
        described[waypointIndex] ||
        !grid.isValid(point.x, point.y, point.z) ||
        grid.getState(point.x, point.y, point.z) !=
            voxel_planner::VoxelClass::POSE_CONDITIONAL) {
        return;
    }
    descriptions.push_back({normal, tangent, waypointIndex});
    described[waypointIndex] = true;
}

std::uint32_t findNearestAllowedEndpointPose(
    const VoxelGrid& grid,
    const voxel_planner::Point3D& point,
    const voxel_planner::Vector3D& requestedTangent,
    const std::vector<voxel_planner::BusbarPose>& poses,
    bool requireCollisionFree = true) {
    constexpr double kAlignmentTolerance = 1e-9;
    const std::size_t voxelIndex = grid.index(
        static_cast<std::uint32_t>(point.x),
        static_cast<std::uint32_t>(point.y),
        static_cast<std::uint32_t>(point.z));

    double bestAlignment = -std::numeric_limits<double>::infinity();
    for (const voxel_planner::BusbarPose& pose : poses) {
        const double tangentLength = std::sqrt(static_cast<double>(
            pose.tx * pose.tx + pose.ty * pose.ty + pose.tz * pose.tz));
        const double alignment =
            (requestedTangent.x * pose.tx +
             requestedTangent.y * pose.ty +
             requestedTangent.z * pose.tz) /
            tangentLength;
        bestAlignment = std::max(bestAlignment, alignment);
    }

    for (const voxel_planner::BusbarPose& pose : poses) {
        const double tangentLength = std::sqrt(static_cast<double>(
            pose.tx * pose.tx + pose.ty * pose.ty + pose.tz * pose.tz));
        const double alignment =
            (requestedTangent.x * pose.tx +
             requestedTangent.y * pose.ty +
             requestedTangent.z * pose.tz) /
            tangentLength;
        if (alignment + kAlignmentTolerance >= bestAlignment &&
            (!requireCollisionFree ||
             grid.isPoseAllowedForSearch(voxelIndex, pose.poseId))) {
            return pose.poseId;
        }
    }
    return std::numeric_limits<std::uint32_t>::max();
}

struct ResolvedEndpoint {
    voxel_planner::Point3D point{};
    std::uint32_t poseId = std::numeric_limits<std::uint32_t>::max();
};

struct FallbackEndpointSnap {
    voxel_planner::Point3D point{};
    std::uint32_t poseId = std::numeric_limits<std::uint32_t>::max();
    bool shifted = false;
};

ResolvedEndpoint resolveEndpointOutsideTerminal(
    const VoxelGrid& grid,
    const voxel_planner::Point3D& rawPoint,
    const voxel_planner::Vector3D& tangent,
    const voxel_planner::Vector3D& preferredRouteTangent,
    float minimumLength,
    int direction,
    const std::vector<voxel_planner::BusbarPose>& poses,
    const char* endpointName,
    float immunityRadius) {
    const voxel_planner::Point3D nominalPoint = shiftedEndpoint(
        rawPoint,
        tangent,
        minimumLength,
        direction);
    if (!grid.isValid(nominalPoint.x, nominalPoint.y, nominalPoint.z)) {
        throw std::invalid_argument(
            std::string("Shifted ") + endpointName +
            " point is outside the voxel map.");
    }
    const std::size_t maximumAdditionalDistance =
        2U * std::max({grid.width(), grid.height(), grid.depth()}) + 2U;
    const auto findAlongTangent =
        [&](const voxel_planner::Vector3D& searchTangent) {
            ResolvedEndpoint resolved;
            voxel_planner::Point3D previousPoint{
                std::numeric_limits<int>::min(),
                std::numeric_limits<int>::min(),
                std::numeric_limits<int>::min()};
            for (std::size_t additional = 0U;
                 additional <= maximumAdditionalDistance;
                 ++additional) {
                const voxel_planner::Point3D candidate = shiftedEndpoint(
                    rawPoint,
                    searchTangent,
                    minimumLength + static_cast<float>(additional),
                    direction);
                if (candidate == previousPoint) {
                    continue;
                }
                previousPoint = candidate;
                if (!grid.isValid(candidate.x, candidate.y, candidate.z)) {
                    break;
                }
                const std::uint32_t poseId = findNearestAllowedEndpointPose(
                    grid,
                    candidate,
                    searchTangent,
                    poses,
                    true);
                if (poseId != std::numeric_limits<std::uint32_t>::max()) {
                    return ResolvedEndpoint{candidate, poseId};
                }
            }
            return resolved;
        };

    ResolvedEndpoint resolved = findAlongTangent(tangent);
    if (resolved.poseId != std::numeric_limits<std::uint32_t>::max()) {
        return resolved;
    }

    struct RankedTangent {
        voxel_planner::Vector3D tangent{};
        double alignment = -1.0;
    };
    std::vector<RankedTangent> fallbackTangents;
    fallbackTangents.reserve(26U);
    for (int z = -1; z <= 1; ++z) {
        for (int y = -1; y <= 1; ++y) {
            for (int x = -1; x <= 1; ++x) {
                if (x == 0 && y == 0 && z == 0) {
                    continue;
                }
                const double length = std::sqrt(static_cast<double>(
                    x * x + y * y + z * z));
                voxel_planner::Vector3D candidateTangent{
                    x / length,
                    y / length,
                    z / length};
                const double requestedAlignment =
                    candidateTangent.x * tangent.x +
                    candidateTangent.y * tangent.y +
                    candidateTangent.z * tangent.z;
                if (requestedAlignment > 1.0 - 1e-9) {
                    continue;
                }
                fallbackTangents.push_back({
                    candidateTangent,
                    candidateTangent.x * preferredRouteTangent.x +
                        candidateTangent.y * preferredRouteTangent.y +
                        candidateTangent.z * preferredRouteTangent.z});
            }
        }
    }
    std::sort(
        fallbackTangents.begin(),
        fallbackTangents.end(),
        [](const RankedTangent& lhs, const RankedTangent& rhs) {
            return lhs.alignment > rhs.alignment;
        });
    for (const RankedTangent& fallback : fallbackTangents) {
        resolved = findAlongTangent(fallback.tangent);
        if (resolved.poseId != std::numeric_limits<std::uint32_t>::max()) {
            return resolved;
        }
    }

    const double nominalDx = static_cast<double>(
        nominalPoint.x - rawPoint.x);
    const double nominalDy = static_cast<double>(
        nominalPoint.y - rawPoint.y);
    const double nominalDz = static_cast<double>(
        nominalPoint.z - rawPoint.z);
    if (nominalDx * nominalDx + nominalDy * nominalDy +
            nominalDz * nominalDz <
        static_cast<double>(immunityRadius) * immunityRadius) {
        const std::uint32_t poseId = findNearestAllowedEndpointPose(
            grid,
            nominalPoint,
            tangent,
            poses,
            false);
        if (poseId != std::numeric_limits<std::uint32_t>::max()) {
            return {nominalPoint, poseId};
        }
    }

    throw std::invalid_argument(
        std::string("No collision-free pose matches the requested ") +
        endpointName + " tangent before the map boundary.");
}

FallbackEndpointSnap snapEndpointToFreeSpaceForFallback(
    const VoxelGrid& grid,
    const voxel_planner::Point3D& seed,
    const voxel_planner::Vector3D& requestedTangent,
    const std::vector<voxel_planner::BusbarPose>& poses,
    const char* endpointName) {
    if (!grid.isValid(seed.x, seed.y, seed.z)) {
        throw std::invalid_argument(
            std::string("Fallback ") + endpointName +
            " point is outside the voxel map.");
    }

    constexpr int kFallbackEndpointSnapRadius = 8;
    const auto isFree = [&](const voxel_planner::Point3D& point) {
        if (!grid.isValid(point.x, point.y, point.z)) {
            return false;
        }
        const std::size_t index = grid.index(
            static_cast<std::uint32_t>(point.x),
            static_cast<std::uint32_t>(point.y),
            static_cast<std::uint32_t>(point.z));
        return !grid.isObstacle(index);
    };

    const auto considerCandidate =
        [&](const voxel_planner::Point3D& candidate,
            FallbackEndpointSnap& best,
            bool& found) {
            if (!isFree(candidate)) {
                return;
            }
            const std::uint32_t poseId = findNearestAllowedEndpointPose(
                grid,
                candidate,
                requestedTangent,
                poses,
                true);
            if (poseId == std::numeric_limits<std::uint32_t>::max()) {
                return;
            }
            const int dx = candidate.x - seed.x;
            const int dy = candidate.y - seed.y;
            const int dz = candidate.z - seed.z;
            const int distanceSquared = dx * dx + dy * dy + dz * dz;
            const int bestDx = best.point.x - seed.x;
            const int bestDy = best.point.y - seed.y;
            const int bestDz = best.point.z - seed.z;
            const int bestDistanceSquared =
                bestDx * bestDx + bestDy * bestDy + bestDz * bestDz;
            if (!found ||
                distanceSquared < bestDistanceSquared ||
                (distanceSquared == bestDistanceSquared &&
                 std::tie(candidate.z, candidate.y, candidate.x) <
                 std::tie(best.point.z, best.point.y, best.point.x))) {
                best.point = candidate;
                best.poseId = poseId;
                best.shifted = !(candidate == seed);
                found = true;
            }
        };

    FallbackEndpointSnap best{seed, std::numeric_limits<std::uint32_t>::max(), false};
    bool found = false;
    considerCandidate(seed, best, found);
    if (found) {
        return best;
    }

    for (int radius = 1; radius <= kFallbackEndpointSnapRadius; ++radius) {
        for (int dz = -radius; dz <= radius; ++dz) {
            for (int dy = -radius; dy <= radius; ++dy) {
                for (int dx = -radius; dx <= radius; ++dx) {
                    if (std::max({
                            std::abs(dx),
                            std::abs(dy),
                            std::abs(dz)}) != radius) {
                        continue;
                    }
                    const voxel_planner::Point3D candidate{
                        seed.x + dx,
                        seed.y + dy,
                        seed.z + dz};
                    considerCandidate(candidate, best, found);
                }
            }
        }
        if (found) {
            return best;
        }
    }

    throw std::invalid_argument(
        std::string("No free-space fallback root could be found for ") +
        endpointName + ".");
}

} // namespace

namespace voxel_planner {

struct ProcessedMap::Impl {
    VoxelGrid jumpGrid;
    VoxelGrid grid;
    std::vector<BusbarPose> poses;
    PlannerConfig config;
    voxel_planner::Point3D morphologyOffset{};
    std::shared_ptr<const std::vector<std::uint64_t>> prefixSum;
    std::uint32_t prefixWidth = 0U;
    std::uint32_t prefixHeight = 0U;
    std::uint32_t prefixDepth = 0U;
    bool poseDataAugmented = false;
    std::mutex poseDataMutex;

    void ensurePoseDataAugmented() {
        if (poseDataAugmented) {
            return;
        }

        std::lock_guard<std::mutex> lock(poseDataMutex);
        if (poseDataAugmented) {
            return;
        }

        poses = generateDiscretePoses(config);
        const module2_morphology::CachedPoseFootprints footprints =
            module2_morphology::precomputePoseFootprints(
                poses,
                config);
        module2_morphology::populateVoxelPoseMasks(
            grid,
            poses,
            footprints);
        poseDataAugmented = true;
    }
};

ProcessedMap::ProcessedMap(std::shared_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}

ProcessedMap::~ProcessedMap() = default;

void ProcessedMap::exportVoxelClassificationToVtk(
    const std::string& filename) const {
    if (!impl_) {
        throw std::invalid_argument(
            "ProcessedMap is empty or has been moved from.");
    }
    if (filename.empty()) {
        throw std::invalid_argument("VTK output filename must not be empty.");
    }
    impl_->ensurePoseDataAugmented();
    module2_morphology::VoxelIO::exportVoxelClassificationToVtk(
        impl_->grid,
        filename);
}

void ProcessedMap::exportVoxelClassificationAndPosesToVtk(
    const std::string& filename) const {
    if (!impl_) {
        throw std::invalid_argument(
            "ProcessedMap is empty or has been moved from.");
    }
    if (filename.empty()) {
        throw std::invalid_argument("VTK output filename must not be empty.");
    }
    impl_->ensurePoseDataAugmented();
    module2_morphology::VoxelIO::
        exportVoxelClassificationAndPosesToVtk(
            impl_->grid,
            filename);
}

ProcessedMap loadMap(
    const std::string& filepath,
    float width,
    float thickness) {
    if (filepath.empty()) {
        throw std::invalid_argument("Map filepath must not be empty.");
    }
    if (!std::isfinite(width) || width <= 0.0F ||
        !std::isfinite(thickness) || thickness <= 0.0F) {
        throw std::invalid_argument(
            "Busbar width and thickness must be finite and positive.");
    }

    PlannerConfig config;
    config.busbar_width = width;
    config.busbar_thickness = thickness;
    if (!hasPositiveConfiguration(config)) {
        throw std::invalid_argument(
            "Physical busbar dimensions, bend factors, and angle step must "
            "be valid and positive.");
    }

    auto impl = std::make_shared<ProcessedMap::Impl>();
    VoxelGrid rawGrid = module2_morphology::VoxelIO::loadVoxelMap(filepath);
    impl->jumpGrid = rawGrid;
    const int kernelWidth = static_cast<int>(std::ceil(width));
    const int kernelThickness = static_cast<int>(std::ceil(thickness));
    module2_morphology::MorphologyResult morphology =
        module2_morphology::VoxelMorphologyEngine::buildCoarseBlockedMap(
            rawGrid,
            module2_morphology::MorphologyKernel{
                kernelWidth,
                kernelThickness,
                kernelThickness});
    impl->grid = std::move(morphology.grid);
    impl->morphologyOffset = morphology.originOffset;
    impl->prefixSum =
        std::make_shared<const std::vector<std::uint64_t>>(
            std::move(morphology.prefixSum));
    impl->prefixWidth = morphology.prefixWidth;
    impl->prefixHeight = morphology.prefixHeight;
    impl->prefixDepth = morphology.prefixDepth;
    impl->grid.attachMorphologyPrefix(
        impl->prefixSum,
        impl->prefixWidth,
        impl->prefixHeight,
        impl->prefixDepth);
    impl->config = config;
    return ProcessedMap(std::move(impl));
}

std::pair<PlanStatus, std::vector<PathResult>> findPaths(
    const ProcessedMap& map,
    const Point3D& start,
    const Point3D& goal,
    int max_paths) {
    if (!map.impl_) {
        throw std::invalid_argument(
            "ProcessedMap is empty or has been moved from.");
    }
    if (max_paths <= 0) {
        throw std::invalid_argument("max_paths must be positive.");
    }
    const EndpointConstraint constraint;
    validateEndpointConstraint(constraint);
    constexpr float kEndpointImmunityRadius = 15.0F;
    const Point3D internalStart =
        addOffset(start, map.impl_->morphologyOffset);
    const Point3D internalGoal =
        addOffset(goal, map.impl_->morphologyOffset);

    if (!map.impl_->jumpGrid.isValid(start.x, start.y, start.z) ||
        !map.impl_->jumpGrid.isValid(goal.x, goal.y, goal.z)) {
        throw std::invalid_argument(
            "Raw start or goal point is outside the voxel map.");
    }
    if (!map.impl_->grid.isValid(
            internalStart.x,
            internalStart.y,
            internalStart.z) ||
        !map.impl_->grid.isValid(
            internalGoal.x,
            internalGoal.y,
            internalGoal.z)) {
        throw std::invalid_argument(
            "Raw start or goal point is outside the voxel map.");
    }

    logProfileMessage(
        "request start=(" + std::to_string(start.x) + ", " +
        std::to_string(start.y) + ", " + std::to_string(start.z) +
        "), goal=(" + std::to_string(goal.x) + ", " +
        std::to_string(goal.y) + ", " + std::to_string(goal.z) +
        "), max_paths=" + std::to_string(max_paths));
    logProfileStageBegin("Raw endpoint resolution");
    const auto rawEndpointResolutionStart = ProfileClock::now();
    const RawEndpointResolution resolvedStart =
        resolveRawEndpointOutsideObstacle(
            map.impl_->jumpGrid,
            start,
            "start");
    const RawEndpointResolution resolvedGoal =
        resolveRawEndpointOutsideObstacle(
            map.impl_->jumpGrid,
            goal,
            "goal");
    const Point3D& jumpSafeStart = resolvedStart.point;
    const Point3D& jumpSafeGoal = resolvedGoal.point;
    const auto rawEndpointResolutionEnd = ProfileClock::now();
    logProfileStageEnd(
        "Raw endpoint resolution",
        rawEndpointResolutionEnd - rawEndpointResolutionStart);
    logProfileMessage(
        "Raw safe_start=(" + std::to_string(jumpSafeStart.x) + ", " +
        std::to_string(jumpSafeStart.y) + ", " +
        std::to_string(jumpSafeStart.z) +
        "), shifted=" + (resolvedStart.shifted ? "true" : "false") +
        ", safe_goal=(" + std::to_string(jumpSafeGoal.x) + ", " +
        std::to_string(jumpSafeGoal.y) + ", " +
        std::to_string(jumpSafeGoal.z) +
        "), shifted=" + (resolvedGoal.shifted ? "true" : "false"));
    logProfileStageBegin("JumpAStar::findPaths");
    const auto jumpStart = ProfileClock::now();
    const std::vector<std::vector<Point3D>> jumpPaths =
            module3_astar::JumpAStar::findPaths(
            map.impl_->jumpGrid,
            jumpSafeStart,
            jumpSafeGoal,
            max_paths);
    const auto jumpEnd = ProfileClock::now();
    logProfileStageEnd("JumpAStar::findPaths", jumpEnd - jumpStart);
    logProfileMessage(
        "JumpAStar paths=" + std::to_string(jumpPaths.size()));
    for (std::size_t index = 0U; index < jumpPaths.size(); ++index) {
        logProfileMessage(
            "JumpAStar path[" + std::to_string(index) +
            "] points=" + std::to_string(jumpPaths[index].size()));
    }

    std::vector<std::vector<Point3D>> topologyHintsInternal;
    topologyHintsInternal.reserve(jumpPaths.size());
    for (std::size_t pathIndex = 0U;
         pathIndex < jumpPaths.size();
         ++pathIndex) {
        const std::vector<Point3D>& centerline = jumpPaths[pathIndex];
        if (!centerline.empty() &&
            centerline.front() == jumpSafeStart &&
            centerline.back() == jumpSafeGoal) {
            std::vector<Point3D> sourceHintPoints;
            sourceHintPoints.reserve(centerline.size() + 2U);
            sourceHintPoints.push_back(start);
            for (const Point3D& point : centerline) {
                if (sourceHintPoints.empty() ||
                    !(sourceHintPoints.back() == point)) {
                    sourceHintPoints.push_back(point);
                }
            }
            if (sourceHintPoints.empty() ||
                !(sourceHintPoints.back() == goal)) {
                sourceHintPoints.push_back(goal);
            }

            HintProjectionResult projection =
                projectTopologyHintToExpandedGrid(
                    map.impl_->grid,
                    sourceHintPoints,
                    map.impl_->morphologyOffset);
            if (projection.points.size() < 2U ||
                projection.successRate <=
                    kMinimumHintProjectionSuccessRatio) {
                logProfileMessage(
                    "Discarded JumpAStar hint[" +
                    std::to_string(pathIndex) +
                    "] after morphology projection; source_points=" +
                    std::to_string(projection.sourcePointCount) +
                    ", success_rate=" +
                    std::to_string(projection.successRate * 100.0) +
                    "%");
                continue;
            }
            logProfileMessage(
                "Projected JumpAStar hint[" +
                std::to_string(pathIndex) +
                "] onto expanded grid; source_points=" +
                std::to_string(projection.sourcePointCount) +
                ", output_points=" +
                std::to_string(projection.points.size()) +
                ", success_rate=" +
                std::to_string(projection.successRate * 100.0) +
                "%, snapped_points=" +
                std::to_string(projection.projectedPointCount) +
                ", filled_points=" +
                std::to_string(projection.filledPointCount) +
                ", corrected_points=" +
                std::to_string(
                    projection.projectedPointCount +
                    projection.filledPointCount));
            topologyHintsInternal.push_back(std::move(projection.points));
        }
    }
    logProfileMessage(
        "SE3 fallback triggered=true; topology_hints=" +
        std::to_string(topologyHintsInternal.size()));
    for (std::size_t index = 0U;
         index < topologyHintsInternal.size();
         ++index) {
        logProfileMessage(
            "Topology hint[" + std::to_string(index) +
            "] points=" +
            std::to_string(topologyHintsInternal[index].size()));
    }

    map.impl_->ensurePoseDataAugmented();

    logProfileStageBegin("SE3 endpoint resolution");
    const auto endpointResolutionStart = ProfileClock::now();
    const double routeDx =
        static_cast<double>(internalGoal.x - internalStart.x);
    const double routeDy =
        static_cast<double>(internalGoal.y - internalStart.y);
    const double routeDz =
        static_cast<double>(internalGoal.z - internalStart.z);
    const double routeLength = std::sqrt(
        routeDx * routeDx + routeDy * routeDy + routeDz * routeDz);
    const Vector3D preferredRouteTangent = routeLength > 1e-9
        ? Vector3D{
            routeDx / routeLength,
            routeDy / routeLength,
            routeDz / routeLength}
        : constraint.start_tangent;
    const ResolvedEndpoint resolvedPoseStart = resolveEndpointOutsideTerminal(
        map.impl_->grid,
        internalStart,
        constraint.start_tangent,
        preferredRouteTangent,
        constraint.min_begin_length,
        1,
        map.impl_->poses,
        "start",
        kEndpointImmunityRadius);
    const ResolvedEndpoint resolvedPoseGoal = resolveEndpointOutsideTerminal(
        map.impl_->grid,
        internalGoal,
        constraint.end_tangent,
        preferredRouteTangent,
        constraint.min_end_length,
        -1,
        map.impl_->poses,
        "end",
        kEndpointImmunityRadius);
    ResolvedEndpoint selectedPoseStart = resolvedPoseStart;
    ResolvedEndpoint selectedPoseGoal = resolvedPoseGoal;
    const Point3D& safeStart = selectedPoseStart.point;
    const Point3D& safeGoal = selectedPoseGoal.point;
    const auto endpointResolutionEnd = ProfileClock::now();
    logProfileStageEnd(
        "SE3 endpoint resolution",
        endpointResolutionEnd - endpointResolutionStart);
    logProfileMessage(
        "SE3 safe_start=(" + std::to_string(safeStart.x) + ", " +
        std::to_string(safeStart.y) + ", " +
        std::to_string(safeStart.z) + "), start_pose=" +
        std::to_string(resolvedPoseStart.poseId) + ", safe_goal=(" +
        std::to_string(safeGoal.x) + ", " +
        std::to_string(safeGoal.y) + ", " +
        std::to_string(safeGoal.z) + "), goal_pose=" +
        std::to_string(resolvedPoseGoal.poseId));

    logProfileMessage(
        "BEGIN CoarseAStar::findPaths; internal work includes action catalog, "
        "backward voxel BFS, swept-volume validation, and diversity retries");
    if (topologyHintsInternal.empty()) {
        logProfileMessage("Topology hints unavailable; coarse fallback uses none");
    }
    const auto runCoarseSearch =
        [&](const ResolvedEndpoint& searchStart,
            const ResolvedEndpoint& searchGoal,
            const std::vector<Point3D>* topologyHint,
            bool allowFullSearch,
            std::size_t maxExpansionsOverride) {
            return module3_astar::CoarseAStar::findPaths(
                map.impl_->grid,
                ::PoseState{searchStart.point, searchStart.poseId},
                ::PoseState{searchGoal.point, searchGoal.poseId},
                max_paths,
                map.impl_->poses,
                map.impl_->config,
                internalStart,
                internalGoal,
                kEndpointImmunityRadius,
                topologyHint,
                allowFullSearch,
                maxExpansionsOverride);
        };
    const auto coarseStart = ProfileClock::now();
    ::PlanningResult result = runCoarseSearch(
        selectedPoseStart,
        selectedPoseGoal,
            topologyHintsInternal.empty()
                ? nullptr
                : &topologyHintsInternal.front(),
        false,
        0U);
    const auto coarseEnd = ProfileClock::now();
    logProfileStageEnd("CoarseAStar::findPaths total", coarseEnd - coarseStart);
    logProfileMessage(
        "CoarseAStar status=" +
        std::to_string(static_cast<int>(result.status)) +
        ", error_code=" +
        std::to_string(static_cast<int>(result.error_code)) +
        ", internal_paths=" + std::to_string(result.paths.size()) +
        ", message=" + result.message);

    const std::size_t requestedPathCount =
        static_cast<std::size_t>(max_paths);
    deduplicatePlanningPaths(result, requestedPathCount);
    const bool retryableCoarseResult =
        result.error_code == ::ErrorCode::NONE ||
        result.error_code == ::ErrorCode::PATH_NOT_FOUND ||
        result.error_code == ::ErrorCode::COMPUTATION_LIMIT_EXCEEDED ||
        result.error_code == ::ErrorCode::START_POINT_BLOCKED ||
        result.error_code == ::ErrorCode::END_POINT_BLOCKED;
    if ((!topologyHintsInternal.empty() || result.paths.empty()) &&
        result.paths.size() < requestedPathCount &&
        retryableCoarseResult) {
        logProfileMessage(
            "CoarseAStar produced " +
            std::to_string(result.paths.size()) +
            " unique paths; starting fallback full coarse search for " +
            std::to_string(
                requestedPathCount - result.paths.size()) +
            " additional paths.");
        const auto fallbackStart = ProfileClock::now();
        constexpr std::size_t kFallbackMaxExpansions = 150000U;
        ::PlanningResult fallbackResult = runCoarseSearch(
            selectedPoseStart,
            selectedPoseGoal,
            nullptr,
            true,
            kFallbackMaxExpansions);
        const auto fallbackEnd = ProfileClock::now();
        logProfileStageEnd(
            "CoarseAStar fallback full search",
            fallbackEnd - fallbackStart);
        const std::size_t beforeFallbackCount = result.paths.size();
        const bool fallbackProducedNoPath = fallbackResult.paths.empty();
        const ::ErrorCode fallbackErrorCode = fallbackResult.error_code;
        const ::PlannerStatus fallbackStatus = fallbackResult.status;
        const std::string fallbackMessage = fallbackResult.message;
        mergeFallbackPlanningPaths(
            result,
            std::move(fallbackResult),
            requestedPathCount);
        const std::size_t fallbackAddedPathCount =
            result.paths.size() >= beforeFallbackCount
                ? result.paths.size() - beforeFallbackCount
                : 0U;
        logProfileMessage(
            "Fallback full search result: status=" +
            std::to_string(static_cast<int>(fallbackStatus)) +
            ", error_code=" +
            std::to_string(static_cast<int>(fallbackErrorCode)) +
            ", paths=" +
            std::to_string(fallbackAddedPathCount) +
            ", message=" + fallbackMessage);

        // Do not leave the primary Hint search's stale 50k-expansion error
        // in place when the fallback is the last completed search.
        if (fallbackAddedPathCount != 0U) {
            result.status = ::PlannerStatus::OK;
            result.error_code = ::ErrorCode::NONE;
            result.message = fallbackMessage.empty()
                ? "Fallback full coarse search produced path(s)."
                : fallbackMessage;
        } else if (result.paths.empty()) {
            result.status = fallbackStatus;
            result.error_code = fallbackErrorCode;
            result.message = fallbackMessage;
        }

        const bool fallbackNeedsEndpointResnap =
            fallbackProducedNoPath ||
            fallbackErrorCode == ::ErrorCode::START_POINT_BLOCKED ||
            fallbackErrorCode == ::ErrorCode::END_POINT_BLOCKED ||
            fallbackErrorCode == ::ErrorCode::COMPUTATION_LIMIT_EXCEEDED;
        if (fallbackNeedsEndpointResnap) {
            logProfileMessage(
                "Fallback full search hit corridor rupture or endpoint block; retrying with "
                "nearest free-space endpoint roots.");
            const FallbackEndpointSnap fallbackStartRoot =
                snapEndpointToFreeSpaceForFallback(
                    map.impl_->grid,
                    selectedPoseStart.point,
                    constraint.start_tangent,
                    map.impl_->poses,
                    "start");
            const FallbackEndpointSnap fallbackGoalRoot =
                snapEndpointToFreeSpaceForFallback(
                    map.impl_->grid,
                    selectedPoseGoal.point,
                    constraint.end_tangent,
                    map.impl_->poses,
                    "goal");
            const ResolvedEndpoint snappedStart{
                fallbackStartRoot.point,
                fallbackStartRoot.poseId};
            const ResolvedEndpoint snappedGoal{
                fallbackGoalRoot.point,
                fallbackGoalRoot.poseId};
            ::PlanningResult snappedResult = runCoarseSearch(
                snappedStart,
                snappedGoal,
                nullptr,
                true,
                kFallbackMaxExpansions);
            const std::size_t snappedResultCount = snappedResult.paths.size();
            if (snappedResultCount != 0U) {
                selectedPoseStart = snappedStart;
                selectedPoseGoal = snappedGoal;
                result = std::move(snappedResult);
                deduplicatePlanningPaths(result, requestedPathCount);
                result.status = ::PlannerStatus::OK;
                result.error_code = ::ErrorCode::NONE;
                logProfileMessage(
                    "Fallback endpoint re-snap succeeded; start=(" +
                    std::to_string(selectedPoseStart.point.x) + ", " +
                    std::to_string(selectedPoseStart.point.y) + ", " +
                    std::to_string(selectedPoseStart.point.z) +
                    "), goal=(" +
                    std::to_string(selectedPoseGoal.point.x) + ", " +
                    std::to_string(selectedPoseGoal.point.y) + ", " +
                    std::to_string(selectedPoseGoal.point.z) +
                    "), added_paths=" +
                    std::to_string(result.paths.size()) +
                    ", snapped_search_paths=" +
                    std::to_string(snappedResultCount) + ".");
            } else if (result.paths.empty()) {
                result.status = snappedResult.status;
                result.error_code = snappedResult.error_code;
                result.message = snappedResult.message;
                logProfileMessage(
                    "Fallback endpoint re-snap failed: status=" +
                    std::to_string(static_cast<int>(snappedResult.status)) +
                    ", error_code=" +
                    std::to_string(static_cast<int>(
                        snappedResult.error_code)) +
                    ", message=" + snappedResult.message);
            }
        }

        logProfileMessage(
            "CoarseAStar fallback added " +
            std::to_string(fallbackAddedPathCount) +
            " unique paths; total=" +
            std::to_string(result.paths.size()) + "/" +
            std::to_string(requestedPathCount));
    }

    if (result.error_code == ::ErrorCode::PATH_NOT_FOUND) {
        return {PlanStatus::NO_PATH, {}};
    }
    if (result.error_code == ::ErrorCode::COMPUTATION_LIMIT_EXCEEDED) {
        if (result.paths.empty()) {
            throw std::runtime_error(result.message.empty()
                ? "SE(3) computation limit exceeded."
                : result.message);
        }
    }
    if (result.error_code == ::ErrorCode::MEMORY_ALLOCATION_FAILED) {
        throw std::bad_alloc();
    }
    if (result.error_code == ::ErrorCode::INVALID_ARGUMENT ||
        result.error_code == ::ErrorCode::START_OUT_OF_BOUNDS ||
        result.error_code == ::ErrorCode::END_OUT_OF_BOUNDS ||
        result.error_code == ::ErrorCode::START_POSE_INVALID ||
        result.error_code == ::ErrorCode::END_POSE_INVALID ||
        result.error_code == ::ErrorCode::START_POINT_BLOCKED ||
        result.error_code == ::ErrorCode::END_POINT_BLOCKED) {
        throw std::invalid_argument(result.message);
    }
    if (result.status == ::PlannerStatus::OK &&
        (result.error_code == ::ErrorCode::NONE ||
         result.error_code ==
            ::ErrorCode::COMPUTATION_LIMIT_EXCEEDED)) {
        logProfileStageBegin("SE3 public result conversion");
        const auto conversionStart = ProfileClock::now();
        std::vector<PathResult> paths;
        paths.reserve(result.paths.size());
        for (::PathResult& internalPath : result.paths) {
            if (internalPath.path.empty() ||
                !(internalPath.path.front() == safeStart) ||
                (!internalPath.goal_tolerance_accepted &&
                 !(internalPath.path.back() == safeGoal))) {
                throw std::runtime_error(
                    "Internal path endpoints do not match shifted endpoints.");
            }
            if (internalPath.goal_tolerance_accepted) {
                const Point3D& acceptedEndpoint = internalPath.path.back();
                const double dx = static_cast<double>(
                    acceptedEndpoint.x - safeGoal.x);
                const double dy = static_cast<double>(
                    acceptedEndpoint.y - safeGoal.y);
                const double dz = static_cast<double>(
                    acceptedEndpoint.z - safeGoal.z);
                const double acceptedDistance = std::sqrt(
                    dx * dx + dy * dy + dz * dz);
                logProfileMessage(
                    "Goal tolerance accepted at distance " +
                    std::to_string(acceptedDistance) +
                    "; returning nearest safe endpoint.");
            }
            PathResult path;
            path.path.reserve(
                internalPath.path.size() +
                static_cast<std::size_t>(
                        std::max({
                        std::abs(safeStart.x - internalStart.x),
                        std::abs(safeStart.y - internalStart.y),
                        std::abs(safeStart.z - internalStart.z)})) +
                static_cast<std::size_t>(
                    std::max({
                        std::abs(internalGoal.x - safeGoal.x),
                        std::abs(internalGoal.y - safeGoal.y),
                        std::abs(internalGoal.z - safeGoal.z)})));
            appendContinuous26(path.path, internalStart);
            appendContinuous26(path.path, safeStart);
            const std::size_t safeStartPathIndex = path.path.size() - 1U;
            std::vector<std::size_t> internalWaypointIndices(
                internalPath.path.size());
            internalWaypointIndices[0] = safeStartPathIndex;
            for (std::size_t index = 1U;
                 index < internalPath.path.size();
                 ++index) {
                appendContinuous26(path.path, internalPath.path[index]);
                internalWaypointIndices[index] = path.path.size() - 1U;
            }
            const std::size_t safeGoalPathIndex = path.path.size() - 1U;
            if (!internalPath.goal_tolerance_accepted) {
                appendContinuous26(path.path, internalGoal);
            }
            path.cost = internalPath.cost +
                pointDistance(internalStart, safeStart) +
                (internalPath.goal_tolerance_accepted
                    ? 0.0
                    : pointDistance(safeGoal, internalGoal));
            path.pose_description.reserve(path.path.size());
            std::vector<bool> described(path.path.size(), false);
            const Vector3D startNormal = unitPoseNormal(
                map.impl_->poses[selectedPoseStart.poseId]);
            const Vector3D startTangent = unitPoseTangent(
                map.impl_->poses[selectedPoseStart.poseId]);
            const Vector3D goalNormal = unitPoseNormal(
                map.impl_->poses[selectedPoseGoal.poseId]);
            const Vector3D goalTangent = unitPoseTangent(
                map.impl_->poses[selectedPoseGoal.poseId]);
            for (std::size_t index = 0U;
                 index <= safeStartPathIndex;
                 ++index) {
                appendPublicConditionalPose(
                    map.impl_->grid,
                    path.path[index],
                    index,
                    startNormal,
                    startTangent,
                    described,
                    path.pose_description);
            }
            for (const ::Pose& internalPose :
                 internalPath.pose_description) {
                if (internalPose.waypointIndex <
                    internalWaypointIndices.size()) {
                    const std::size_t outputWaypointIndex =
                        internalWaypointIndices[internalPose.waypointIndex];
                    appendPublicConditionalPose(
                        map.impl_->grid,
                        path.path[outputWaypointIndex],
                        outputWaypointIndex,
                        internalPose.normal,
                        internalPose.tangent,
                        described,
                        path.pose_description);
                }
            }
            for (std::size_t index = safeGoalPathIndex;
                 index < path.path.size();
                 ++index) {
                appendPublicConditionalPose(
                    map.impl_->grid,
                    path.path[index],
                    index,
                    goalNormal,
                    goalTangent,
                    described,
                    path.pose_description);
            }
            translatePathToPublicCoordinates(
                path.path,
                map.impl_->morphologyOffset);
            paths.push_back(std::move(path));
        }
        std::sort(
            paths.begin(),
            paths.end(),
            [](const PathResult& lhs, const PathResult& rhs) {
                return lhs.cost < rhs.cost;
            });
        const auto conversionEnd = ProfileClock::now();
        logProfileStageEnd(
            "SE3 public result conversion",
            conversionEnd - conversionStart);
        logProfileMessage(
            "returning SE3 paths=" + std::to_string(paths.size()));
        return {PlanStatus::OK, std::move(paths)};
    }
    throw std::runtime_error(result.message);
}

} // namespace voxel_planner
