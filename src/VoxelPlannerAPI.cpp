#include "VoxelPlannerAPI.h"

#include "Module2_Morphology/PoseMaskPopulation.h"
#include "Module2_Morphology/VoxelIO.h"
#include "Module3_AStar/CoarseAStar.h"
#include "Module3_AStar/JumpAStar.h"
#include "Module2_Morphology/PoseGenerator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace {

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

std::uint32_t findNearestAllowedPoseForTangent(
    const VoxelGrid& grid,
    const Point3D& point,
    const Vector3D& requestedTangent,
    const std::vector<voxel_planner::BusbarPose>& poses) {
    if (!grid.isValid(point.x, point.y, point.z)) {
        return std::numeric_limits<std::uint32_t>::max();
    }
    const std::size_t voxelIndex = grid.index(
        static_cast<std::uint32_t>(point.x),
        static_cast<std::uint32_t>(point.y),
        static_cast<std::uint32_t>(point.z));
    double bestAlignment = -std::numeric_limits<double>::infinity();
    std::uint32_t bestPoseId = std::numeric_limits<std::uint32_t>::max();
    for (const voxel_planner::BusbarPose& pose : poses) {
        if (!grid.isPoseAllowed(voxelIndex, pose.poseId)) {
            continue;
        }
        const double tangentLength = std::sqrt(static_cast<double>(
            pose.tx * pose.tx + pose.ty * pose.ty + pose.tz * pose.tz));
        if (tangentLength <= 1e-12) {
            continue;
        }
        const double alignment =
            (requestedTangent.x * pose.tx +
             requestedTangent.y * pose.ty +
             requestedTangent.z * pose.tz) / tangentLength;
        if (alignment > bestAlignment) {
            bestAlignment = alignment;
            bestPoseId = pose.poseId;
        }
    }
    return bestPoseId;
}

void appendJumpPathConditionalPoses(
    const VoxelGrid& grid,
    const std::vector<Point3D>& path,
    const std::vector<voxel_planner::BusbarPose>& poses,
    std::vector<voxel_planner::PoseDescription>& descriptions) {
    descriptions.reserve(path.size());
    for (std::size_t index = 0U; index < path.size(); ++index) {
        const Point3D& point = path[index];
        if (!grid.isValid(point.x, point.y, point.z) ||
            grid.getState(point.x, point.y, point.z) !=
                voxel_planner::VoxelClass::POSE_CONDITIONAL) {
            continue;
        }

        Point3D delta{};
        if (index + 1U < path.size()) {
            delta = {
                path[index + 1U].x - point.x,
                path[index + 1U].y - point.y,
                path[index + 1U].z - point.z};
        } else if (index > 0U) {
            delta = {
                point.x - path[index - 1U].x,
                point.y - path[index - 1U].y,
                point.z - path[index - 1U].z};
        } else {
            delta = {1, 0, 0};
        }
        const double length = std::sqrt(static_cast<double>(
            delta.x * delta.x +
            delta.y * delta.y +
            delta.z * delta.z));
        if (length <= 1e-12) {
            continue;
        }
        const Vector3D requestedTangent{
            delta.x / length,
            delta.y / length,
            delta.z / length};
        const std::uint32_t poseId = findNearestAllowedPoseForTangent(
            grid,
            point,
            requestedTangent,
            poses);
        if (poseId == std::numeric_limits<std::uint32_t>::max()) {
            continue;
        }
        const voxel_planner::BusbarPose& pose = poses[poseId];
        descriptions.push_back({
            unitPoseNormal(pose),
            unitPoseTangent(pose),
            index});
    }
}

float jumpPathCost(const std::vector<Point3D>& path) {
    double cost = 0.0;
    for (std::size_t index = 1U; index < path.size(); ++index) {
        cost += pointDistance(path[index - 1U], path[index]);
    }
    return static_cast<float>(cost);
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

} // namespace

namespace voxel_planner {

struct ProcessedMap::Impl {
    VoxelGrid grid;
    std::vector<BusbarPose> poses;
    PlannerConfig config;
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
    impl->grid = module2_morphology::VoxelIO::loadVoxelMap(filepath);
    impl->poses = generateDiscretePoses(config);
    impl->config = config;
    const module2_morphology::CachedPoseFootprints footprints =
        module2_morphology::precomputePoseFootprints(
            impl->poses,
            impl->config);
    module2_morphology::populateVoxelPoseMasks(
        impl->grid,
        impl->poses,
        footprints);
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

    if (!map.impl_->grid.isValid(start.x, start.y, start.z) ||
        !map.impl_->grid.isValid(goal.x, goal.y, goal.z)) {
        throw std::invalid_argument(
            "Raw start or goal point is outside the voxel map.");
    }

    const std::vector<std::vector<Point3D>> jumpPaths =
        module3_astar::JumpAStar::findPaths(
            map.impl_->grid,
            start,
            goal,
            max_paths);
    if (!jumpPaths.empty()) {
        std::vector<PathResult> paths;
        paths.reserve(jumpPaths.size());
        for (const std::vector<Point3D>& centerline : jumpPaths) {
            if (centerline.empty() ||
                !(centerline.front() == start) ||
                !(centerline.back() == goal)) {
                continue;
            }
            PathResult path;
            path.path = centerline;
            path.cost = jumpPathCost(centerline);
            appendJumpPathConditionalPoses(
                map.impl_->grid,
                path.path,
                map.impl_->poses,
                path.pose_description);
            paths.push_back(std::move(path));
        }
        if (!paths.empty()) {
            std::sort(
                paths.begin(),
                paths.end(),
                [](const PathResult& lhs, const PathResult& rhs) {
                    return lhs.cost < rhs.cost;
                });
            return {PlanStatus::OK, std::move(paths)};
        }
    }

    const double routeDx = static_cast<double>(goal.x - start.x);
    const double routeDy = static_cast<double>(goal.y - start.y);
    const double routeDz = static_cast<double>(goal.z - start.z);
    const double routeLength = std::sqrt(
        routeDx * routeDx + routeDy * routeDy + routeDz * routeDz);
    const Vector3D preferredRouteTangent = routeLength > 1e-9
        ? Vector3D{
            routeDx / routeLength,
            routeDy / routeLength,
            routeDz / routeLength}
        : constraint.start_tangent;
    const ResolvedEndpoint resolvedStart = resolveEndpointOutsideTerminal(
        map.impl_->grid,
        start,
        constraint.start_tangent,
        preferredRouteTangent,
        constraint.min_begin_length,
        1,
        map.impl_->poses,
        "start",
        kEndpointImmunityRadius);
    const ResolvedEndpoint resolvedGoal = resolveEndpointOutsideTerminal(
        map.impl_->grid,
        goal,
        constraint.end_tangent,
        preferredRouteTangent,
        constraint.min_end_length,
        -1,
        map.impl_->poses,
        "end",
        kEndpointImmunityRadius);
    const Point3D& safeStart = resolvedStart.point;
    const Point3D& safeGoal = resolvedGoal.point;

    ::PlanningResult result = module3_astar::CoarseAStar::findPaths(
        map.impl_->grid,
        ::PoseState{safeStart, resolvedStart.poseId},
        ::PoseState{safeGoal, resolvedGoal.poseId},
        max_paths,
        map.impl_->poses,
        map.impl_->config,
        start,
        goal,
        kEndpointImmunityRadius);
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
        std::vector<PathResult> paths;
        paths.reserve(result.paths.size());
        for (::PathResult& internalPath : result.paths) {
            if (internalPath.path.empty() ||
                !(internalPath.path.front() == safeStart) ||
                !(internalPath.path.back() == safeGoal)) {
                throw std::runtime_error(
                    "Internal path endpoints do not match shifted endpoints.");
            }
            PathResult path;
            path.path.reserve(
                internalPath.path.size() +
                static_cast<std::size_t>(
                    std::max({
                        std::abs(safeStart.x - start.x),
                        std::abs(safeStart.y - start.y),
                        std::abs(safeStart.z - start.z)})) +
                static_cast<std::size_t>(
                    std::max({
                        std::abs(goal.x - safeGoal.x),
                        std::abs(goal.y - safeGoal.y),
                        std::abs(goal.z - safeGoal.z)})));
            appendContinuous26(path.path, start);
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
            appendContinuous26(path.path, goal);
            path.cost = internalPath.cost +
                pointDistance(start, safeStart) +
                pointDistance(safeGoal, goal);
            path.pose_description.reserve(path.path.size());
            std::vector<bool> described(path.path.size(), false);
            const Vector3D startNormal = unitPoseNormal(
                map.impl_->poses[resolvedStart.poseId]);
            const Vector3D startTangent = unitPoseTangent(
                map.impl_->poses[resolvedStart.poseId]);
            const Vector3D goalNormal = unitPoseNormal(
                map.impl_->poses[resolvedGoal.poseId]);
            const Vector3D goalTangent = unitPoseTangent(
                map.impl_->poses[resolvedGoal.poseId]);
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
            paths.push_back(std::move(path));
        }
        std::sort(
            paths.begin(),
            paths.end(),
            [](const PathResult& lhs, const PathResult& rhs) {
                return lhs.cost < rhs.cost;
            });
        return {PlanStatus::OK, std::move(paths)};
    }
    throw std::runtime_error(result.message);
}

} // namespace voxel_planner
