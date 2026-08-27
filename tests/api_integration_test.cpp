#include "VoxelPlannerAPI.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>

namespace {

bool isContinuous26Neighborhood(
    const std::vector<voxel_planner::Point3D>& path) {
    for (std::size_t i = 1U; i < path.size(); ++i) {
        const int dx = std::abs(path[i].x - path[i - 1U].x);
        const int dy = std::abs(path[i].y - path[i - 1U].y);
        const int dz = std::abs(path[i].z - path[i - 1U].z);
        if (dx > 1 || dy > 1 || dz > 1 || dx + dy + dz == 0) {
            return false;
        }
    }
    return true;
}

bool isUnit(const voxel_planner::Vector3D& vector) {
    const double length = std::sqrt(
        vector.x * vector.x +
        vector.y * vector.y +
        vector.z * vector.z);
    return std::abs(length - 1.0) < 1e-6;
}

bool isOrthogonal(
    const voxel_planner::Vector3D& lhs,
    const voxel_planner::Vector3D& rhs) {
    return std::abs(
        lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z) < 1e-6;
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <path_to_map.txt_or_vtk>\n";
        return 2;
    }

    constexpr float width = 16.0F;
    constexpr float thickness = 2.0F;

    const voxel_planner::ProcessedMap map =
        voxel_planner::loadMap(argv[1], width, thickness);
    const voxel_planner::Point3D rawStart{10, 20, 8};
    const voxel_planner::Point3D rawGoal{35, 30, 8};
    const auto result =
        voxel_planner::findPaths(
            map,
            rawStart,
            rawGoal,
            1);
    const std::vector<voxel_planner::PathResult>& paths = result.second;

    if (result.first != voxel_planner::PlanStatus::OK || paths.empty()) {
        std::cerr << "[FAIL] No path returned by the public API.\n";
        return 1;
    }
    if (!std::is_sorted(
            paths.begin(),
            paths.end(),
            [](const voxel_planner::PathResult& lhs,
               const voxel_planner::PathResult& rhs) {
                return lhs.cost < rhs.cost;
            })) {
        std::cerr << "[FAIL] Returned paths are not sorted by cost.\n";
        return 1;
    }

    bool foundConditionalPose = false;
    for (const voxel_planner::PathResult& path : paths) {
        if (path.path.empty() ||
            !(path.path.front() == rawStart) ||
            !(path.path.back() == rawGoal) ||
            !isContinuous26Neighborhood(path.path)) {
            std::cerr << "[FAIL] Invalid 26-neighborhood path.\n";
            return 1;
        }
        for (const voxel_planner::PoseDescription& pose :
             path.pose_description) {
            if (!isUnit(pose.normal) ||
                !isUnit(pose.tangent) ||
                !isOrthogonal(pose.normal, pose.tangent)) {
                std::cerr << "[FAIL] Invalid conditional pose export.\n";
                return 1;
            }
            foundConditionalPose = true;
        }
    }
    if (!foundConditionalPose) {
        std::cerr << "[FAIL] No conditional waypoint pose was exported.\n";
        return 1;
    }

    const auto repeatedResult =
        voxel_planner::findPaths(
            map,
            rawStart,
            rawGoal,
            1);
    const std::vector<voxel_planner::PathResult>& repeatedPaths =
        repeatedResult.second;
    if (repeatedResult.first != voxel_planner::PlanStatus::OK ||
        repeatedPaths.size() != paths.size() ||
        repeatedPaths.front().path != paths.front().path ||
        std::abs(repeatedPaths.front().cost - paths.front().cost) > 1e-9) {
        std::cerr << "[FAIL] ProcessedMap reuse changed the result.\n";
        return 1;
    }

    std::cout << "[PASS] Public API integration test\n"
              << "Map: " << argv[1] << "\n"
              << "Returned paths: " << paths.size() << "\n";
    for (std::size_t i = 0U; i < paths.size(); ++i) {
        std::cout << "Path " << (i + 1U)
                  << ": cost=" << paths[i].cost
                  << ", waypoints=" << paths[i].path.size()
                  << ", poses="
                  << paths[i].pose_description.size()
                  << "\n";
    }
    return 0;
}
