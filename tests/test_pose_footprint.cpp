#include "Module2_Morphology/PoseGenerator.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "[FAIL] " << message << "\n";
        std::exit(1);
    }
}

} // namespace

int main() {
    voxel_planner::PlannerConfig config;
    config.busbar_width = 35.0F;
    config.busbar_thickness = 5.0F;
    const voxel_planner::BusbarPose flatPose{
        1, 0, 0,
        0.0F, 0.0F, 1.0F,
        0U, 0};

    const std::vector<voxel_planner::VoxelOffset> footprint =
        voxel_planner::generatePoseFootprint(flatPose, config);
    require(!footprint.empty(), "Flat-pose footprint is empty.");

    int minX = std::numeric_limits<int>::max();
    int maxX = std::numeric_limits<int>::min();
    int minY = std::numeric_limits<int>::max();
    int maxY = std::numeric_limits<int>::min();
    int minZ = std::numeric_limits<int>::max();
    int maxZ = std::numeric_limits<int>::min();
    for (const voxel_planner::VoxelOffset& offset : footprint) {
        minX = std::min(minX, offset.dx);
        maxX = std::max(maxX, offset.dx);
        minY = std::min(minY, offset.dy);
        maxY = std::max(maxY, offset.dy);
        minZ = std::min(minZ, offset.dz);
        maxZ = std::max(maxZ, offset.dz);
    }

    require(minX == 0 && maxX == 0,
            "Static footprint must not extend along the tangent.");
    require(minY == -17 && maxY == 17,
            "Width does not rasterize to the expected Y extent.");
    require(minZ == -2 && maxZ == 2,
            "Thickness does not rasterize to the expected Z extent.");
    require(footprint.size() == 175U,
            "Flat footprint does not contain the expected voxel offsets.");

    config.busbar_width = 16.0F;
    config.busbar_thickness = 2.0F;
    const std::vector<voxel_planner::VoxelOffset> evenFootprint =
        voxel_planner::generatePoseFootprint(flatPose, config);
    require(evenFootprint.size() == 32U,
            "Even busbar dimensions must rasterize to exactly W * H samples.");
    for (const voxel_planner::VoxelOffset& offset : evenFootprint) {
        require(offset.dx == 0,
                "Static footprint contains a longitudinal extrusion.");
    }

    std::cout << "[PASS] offsets=" << footprint.size()
              << " x=[" << minX << ',' << maxX << ']'
              << " y=[" << minY << ',' << maxY << ']'
              << " z=[" << minZ << ',' << maxZ << "]\n";
    return 0;
}
