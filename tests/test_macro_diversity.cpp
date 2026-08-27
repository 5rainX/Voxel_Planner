#include "Module1_Main/Types.h"
#include "Module2_Morphology/PoseMaskPopulation.h"
#include "Module3_AStar/CoarseAStar.h"
#include "Module2_Morphology/PoseGenerator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "[FAIL] " << message << "\n";
        std::exit(1);
    }
}

double distanceSquared(const Point3D& lhs, const Point3D& rhs) {
    const double dx = static_cast<double>(lhs.x - rhs.x);
    const double dy = static_cast<double>(lhs.y - rhs.y);
    const double dz = static_cast<double>(lhs.z - rhs.z);
    return dx * dx + dy * dy + dz * dz;
}

std::vector<Point3D> interiorPoints(
    const std::vector<Point3D>& path,
    const Point3D& start,
    const Point3D& goal,
    double endpointRadius) {
    const double radiusSquared = endpointRadius * endpointRadius;
    std::vector<Point3D> points;
    for (const Point3D& point : path) {
        if (distanceSquared(point, start) > radiusSquared &&
            distanceSquared(point, goal) > radiusSquared) {
            points.push_back(point);
        }
    }
    return points;
}

double directedAverageSeparation(
    const std::vector<Point3D>& source,
    const std::vector<Point3D>& target) {
    double total = 0.0;
    for (const Point3D& point : source) {
        double nearestSquared = std::numeric_limits<double>::infinity();
        for (const Point3D& other : target) {
            nearestSquared = std::min(
                nearestSquared,
                distanceSquared(point, other));
        }
        total += std::sqrt(nearestSquared);
    }
    return total / static_cast<double>(source.size());
}

double averageSeparation(
    const std::vector<Point3D>& lhs,
    const std::vector<Point3D>& rhs,
    const Point3D& start,
    const Point3D& goal,
    double endpointRadius) {
    const std::vector<Point3D> lhsInterior = interiorPoints(
        lhs, start, goal, endpointRadius);
    const std::vector<Point3D> rhsInterior = interiorPoints(
        rhs, start, goal, endpointRadius);
    require(!lhsInterior.empty() && !rhsInterior.empty(),
            "A candidate path has no interior corridor samples.");
    return 0.5 * (
        directedAverageSeparation(lhsInterior, rhsInterior) +
        directedAverageSeparation(rhsInterior, lhsInterior));
}

double overlapWithPreviousPaths(
    const std::vector<Point3D>& candidate,
    const std::vector<PathResult>& paths,
    std::size_t candidateIndex) {
    const auto key = [](const Point3D& point) {
        return static_cast<std::uint64_t>(point.x) |
               (static_cast<std::uint64_t>(point.y) << 20U) |
               (static_cast<std::uint64_t>(point.z) << 40U);
    };
    std::unordered_set<std::uint64_t> accepted;
    for (std::size_t i = 0U; i < candidateIndex; ++i) {
        for (const Point3D& point : paths[i].path) {
            accepted.insert(key(point));
        }
    }
    std::unordered_set<std::uint64_t> uniqueCandidate;
    std::size_t shared = 0U;
    for (const Point3D& point : candidate) {
        const std::uint64_t voxel = key(point);
        if (uniqueCandidate.insert(voxel).second &&
            accepted.find(voxel) != accepted.end()) {
            ++shared;
        }
    }
    return static_cast<double>(shared) / uniqueCandidate.size();
}

std::uint32_t findPositiveXPose(
    const std::vector<voxel_planner::BusbarPose>& poses) {
    for (const voxel_planner::BusbarPose& pose : poses) {
        if (pose.tx == 1 && pose.ty == 0 && pose.tz == 0 &&
            std::abs(pose.nx) < 1e-5F &&
            std::abs(pose.ny) < 1e-5F &&
            pose.nz > 0.999F) {
            return pose.poseId;
        }
    }
    return std::numeric_limits<std::uint32_t>::max();
}

} // namespace

int main() {
    voxel_planner::PlannerConfig config;
    config.busbar_width = 1.0F;
    config.busbar_thickness = 1.0F;
    config.flat_bend_factor = 1.5F;
    config.vertical_bend_factor = 1.0F;
    // Four roll samples keep this test focused on corridor blocking rather
    // than stress-testing the full 624-state production pose table.
    config.angle_step_deg = 90;
    config.path_blocking_radius = 4.0F;
    config.max_overlap_ratio = 0.3F;

    const std::vector<voxel_planner::BusbarPose> poses =
        voxel_planner::generateDiscretePoses(config);
    const std::uint32_t endpointPoseId = findPositiveXPose(poses);
    require(endpointPoseId < poses.size(),
            "The +X endpoint pose was not generated.");

    VoxelGrid map(40, 40, 13);
    const module2_morphology::CachedPoseFootprints footprints =
        module2_morphology::precomputePoseFootprints(
            poses,
            config);
    module2_morphology::populateVoxelPoseMasks(
        map,
        poses,
        footprints);

    const Point3D start{6, 20, 6};
    const Point3D goal{33, 20, 6};
    const PlanningResult result = module3_astar::CoarseAStar::findPaths(
        map,
        PoseState{start, endpointPoseId},
        PoseState{goal, endpointPoseId},
        3,
        poses,
        config);
    if (result.status != PlannerStatus::OK) {
        std::cerr << "Search error="
                  << static_cast<int>(result.error_code)
                  << " message=" << result.message
                  << " returned_paths=" << result.paths.size() << "\n";
    }
    require(result.status == PlannerStatus::OK,
            "Macro-diversity search returned NO_PATH.");
    if (result.paths.size() != 3U) {
        std::cerr << "Returned paths=" << result.paths.size();
        for (std::size_t i = 0U; i < result.paths.size(); ++i) {
            std::cerr << " path" << (i + 1U)
                      << "_cost=" << result.paths[i].cost
                      << " path" << (i + 1U)
                      << "_points=" << result.paths[i].path.size();
        }
        std::cerr << "\n";
    }
    require(result.paths.size() == 3U,
            "Macro-diversity search did not return three paths.");
    const Point3D originalMidpoint =
        result.paths.front().path[result.paths.front().path.size() / 2U];
    require(map.getState(
                originalMidpoint.x,
                originalMidpoint.y,
                originalMidpoint.z) != VoxelState::BLOCKED,
            "Corridor blocking leaked into the caller's original grid.");
    require(std::is_sorted(
                result.paths.begin(),
                result.paths.end(),
                [](const PathResult& lhs, const PathResult& rhs) {
                    return lhs.cost < rhs.cost;
                }),
            "Returned paths are not sorted by physical cost.");
    double maximumOverlap = 0.0;
    for (std::size_t i = 1U; i < result.paths.size(); ++i) {
        maximumOverlap = std::max(
            maximumOverlap,
            overlapWithPreviousPaths(
                result.paths[i].path,
                result.paths,
                i));
    }
    require(maximumOverlap <= config.max_overlap_ratio,
            "A returned path exceeds max_overlap_ratio.");

    double minimumAverageSeparation =
        std::numeric_limits<double>::infinity();
    for (std::size_t i = 0U; i < result.paths.size(); ++i) {
        for (std::size_t j = i + 1U; j < result.paths.size(); ++j) {
            minimumAverageSeparation = std::min(
                minimumAverageSeparation,
                averageSeparation(
                    result.paths[i].path,
                    result.paths[j].path,
                    start,
                    goal,
                    config.path_blocking_radius));
        }
    }
    require(minimumAverageSeparation >=
                config.path_blocking_radius * 0.75F,
            "Returned paths are only microscopic parallel deviations.");

    std::cout << "[PASS] paths=" << result.paths.size()
              << " minimum_average_separation="
              << minimumAverageSeparation
              << " maximum_overlap=" << maximumOverlap;
    for (std::size_t i = 0U; i < result.paths.size(); ++i) {
        std::cout << " path" << (i + 1U)
                  << "_cost=" << result.paths[i].cost
                  << " path" << (i + 1U)
                  << "_points=" << result.paths[i].path.size();
    }
    std::cout << "\n";
    return 0;
}
