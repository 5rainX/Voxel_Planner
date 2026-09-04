#include "VoxelPlannerAPI.h"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void printUsage(const char* executableName) {
    std::cerr
        << "Usage: " << executableName
        << " <voxel_map.txt_or_vtk>"
        << " <start_x> <start_y> <start_z>"
        << " <goal_x> <goal_y> <goal_z>"
        << " <max_paths>"
        << " [width thickness]"
        << " [flat_bend_factor vertical_bend_factor]\n";
}

int parseIntArgument(const char* value, const char* name) {
    try {
        std::size_t consumed = 0U;
        const int parsed = std::stoi(value, &consumed);
        if (value[consumed] != '\0') {
            throw std::invalid_argument("trailing characters");
        }
        return parsed;
    } catch (const std::exception&) {
        throw std::invalid_argument(
            std::string("Invalid integer argument: ") + name);
    }
}

float parseFloatArgument(const char* value, const char* name) {
    try {
        std::size_t consumed = 0U;
        const float parsed = std::stof(value, &consumed);
        if (value[consumed] != '\0') {
            throw std::invalid_argument("trailing characters");
        }
        return parsed;
    } catch (const std::exception&) {
        throw std::invalid_argument(
            std::string("Invalid numeric argument: ") + name);
    }
}

void printDensePath(
    std::size_t pathIndex,
    const voxel_planner::PathResult& path) {
    std::cout << "Path " << (pathIndex + 1U)
              << ": cost=" << path.cost
              << ", waypoints=" << path.path.size()
              << ", conditional_poses="
              << path.pose_description.size() << "\n";
    for (std::size_t pointIndex = 0U;
         pointIndex < path.path.size();
         ++pointIndex) {
        const voxel_planner::Point3D& point = path.path[pointIndex];
        std::cout << "  voxel " << pointIndex
                  << ": " << point.x
                  << ' ' << point.y
                  << ' ' << point.z << "\n";
    }
}

} // namespace

/**
 * @brief CLI wrapper for the public planner facade.
 */
int main(int argc, char* argv[]) {
    if (argc != 9 && argc != 11 && argc != 13) {
        printUsage(argv[0]);
        return 2;
    }

    try {
        const voxel_planner::Point3D start{
            parseIntArgument(argv[2], "start_x"),
            parseIntArgument(argv[3], "start_y"),
            parseIntArgument(argv[4], "start_z")};
        const voxel_planner::Point3D goal{
            parseIntArgument(argv[5], "goal_x"),
            parseIntArgument(argv[6], "goal_y"),
            parseIntArgument(argv[7], "goal_z")};
        const int maxPaths = parseIntArgument(argv[8], "max_paths");
        const float width = argc >= 11
            ? parseFloatArgument(argv[9], "width")
            : 35.0F;
        const float thickness = argc >= 11
            ? parseFloatArgument(argv[10], "thickness")
            : 5.0F;
        const auto begin = std::chrono::steady_clock::now();
        const voxel_planner::ProcessedMap map = argc == 13
            ? voxel_planner::loadMap(
                  argv[1],
                  width,
                  thickness,
                  parseFloatArgument(argv[11], "flat_bend_factor"),
                  parseFloatArgument(argv[12], "vertical_bend_factor"))
            : voxel_planner::loadMap(argv[1], width, thickness);
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
            printDensePath(pathIndex, paths[pathIndex]);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Planner error: " << error.what() << "\n";
        return 1;
    }
}
