#include "Module1_Main/Types.h"
#include "Module3_AStar/JumpAStar.h"

#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "[FAIL] " << message << "\n";
        std::exit(1);
    }
}

void requireValidPath(
    const VoxelGrid& map,
    const std::vector<Point3D>& path,
    Point3D start,
    Point3D goal) {
    require(!path.empty(), "Jump A* returned an empty path.");
    require(path.front() == start, "Path start does not match the request.");
    require(path.back() == goal, "Path goal does not match the request.");
    for (const Point3D& point : path) {
        require(
            map.isValid(point.x, point.y, point.z),
            "Path contains an out-of-bounds point.");
        require(
            !map.isRawObstacle(map.index(
                static_cast<std::uint32_t>(point.x),
                static_cast<std::uint32_t>(point.y),
                static_cast<std::uint32_t>(point.z))),
            "Path contains a raw obstacle.");
    }
}

} // namespace

int main() {
    VoxelGrid largeOpenMap(1280U, 1103U, 142U);
    const Point3D start{0, 0, 0};
    const Point3D goal{1279, 1102, 141};
    const std::vector<Point3D> directPath =
        module3_astar::JumpAStar::findPath(
            largeOpenMap,
            start,
            goal);
    requireValidPath(largeOpenMap, directPath, start, goal);

    VoxelGrid boundaryMap(4U, 4U, 4U);
    const std::vector<Point3D> boundaryPath =
        module3_astar::JumpAStar::findPath(
            boundaryMap,
            Point3D{0, 0, 0},
            Point3D{3, 0, 0});
    requireValidPath(
        boundaryMap,
        boundaryPath,
        Point3D{0, 0, 0},
        Point3D{3, 0, 0});

    VoxelGrid obstacleMap(9U, 5U, 5U);
    for (int y = 0; y < 5; ++y) {
        for (int z = 0; z < 5; ++z) {
            obstacleMap.setRawObstacle(4, y, z);
        }
    }
    const std::vector<Point3D> noPath =
        module3_astar::JumpAStar::findPath(
            obstacleMap,
            Point3D{1, 2, 2},
            Point3D{7, 2, 2});
    require(noPath.empty(), "A fully separated map must have no path.");

    std::cout << "[PASS] jump_astar large-open, boundary, and obstacle cases\n";
    return 0;
}
