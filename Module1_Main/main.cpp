#include "VoxelPlannerAPI.h"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>

/**
 * @brief CLI wrapper for the public planner facade.
 */
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <voxel_map.txt_or_vtk>\n";
        return 2;
    }

    std::cout << "Busbar width and thickness [default 35 5]: ";
    float width = 35.0F;
    float thickness = 5.0F;
    std::string dimensions;
    std::getline(std::cin, dimensions);
    if (!dimensions.empty()) {
        try {
            const std::size_t separator = dimensions.find(' ');
            width = std::stof(dimensions.substr(0, separator));
            thickness = std::stof(
                dimensions.substr(separator + 1U));
        } catch (const std::exception&) {
            std::cerr << "Invalid physical dimensions.\n";
            return 2;
        }
    }

    voxel_planner::Point3D start;
    voxel_planner::Point3D goal;
    std::cout << "Start point (x y z): ";
    if (!(std::cin >> start.x >> start.y >> start.z)) {
        std::cerr << "Invalid start point.\n";
        return 2;
    }
    std::cout << "Goal point (x y z): ";
    if (!(std::cin >> goal.x >> goal.y >> goal.z)) {
        std::cerr << "Invalid goal point.\n";
        return 2;
    }

    int maxPaths = 0;
    std::cout << "Max paths: ";
    if (!(std::cin >> maxPaths)) {
        std::cerr << "Invalid max_paths.\n";
        return 2;
    }

    try {
        const auto begin = std::chrono::steady_clock::now();
        const voxel_planner::ProcessedMap map =
            voxel_planner::loadMap(argv[1], width, thickness);
        const auto result =
            voxel_planner::findPaths(
                map,
                start,
                goal,
                maxPaths);
        const std::vector<voxel_planner::PathResult>& paths = result.second;
        const auto end = std::chrono::steady_clock::now();
        const double elapsedMs =
            std::chrono::duration<double, std::milli>(end - begin).count();

        std::cout << "Status: "
                  << (result.first == voxel_planner::PlanStatus::OK
                          ? "OK"
                          : "NO_PATH")
                  << "\nElapsedMs: " << elapsedMs
                  << "\nPaths: " << paths.size() << "\n";

        for (std::size_t pathIndex = 0U;
             pathIndex < paths.size();
             ++pathIndex) {
            const voxel_planner::PathResult& path = paths[pathIndex];
            std::cout << "Path " << (pathIndex + 1U)
                      << ": cost=" << path.cost
                      << ", waypoints=" << path.path.size()
                      << ", conditional_poses="
                      << path.pose_description.size() << "\n";
        }
        return result.first == voxel_planner::PlanStatus::OK ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "Planner error: " << error.what() << "\n";
        return 1;
    }
}
