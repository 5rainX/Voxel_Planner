#include "Module1_Main/Types.h"
#include "Module2_Morphology/PoseMaskPopulation.h"
#include "Module2_Morphology/PoseGenerator.h"

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

} // namespace

int main() {
    constexpr int gridSize = 8;
    VoxelGrid grid(gridSize, gridSize, gridSize);

    // A complete X=3 wall intersects some width orientations at the adjacent
    // anchor while leaving other cross-section orientations collision-free.
    for (int z = 0; z < gridSize; ++z) {
        for (int y = 0; y < gridSize; ++y) {
            grid.setRawObstacle(3, y, z);
        }
    }

    // Surround one free anchor with all 26 neighboring obstacles.
    const Point3D trapped{5, 3, 3};
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) {
                    continue;
                }
                grid.setRawObstacle(
                    trapped.x + dx,
                    trapped.y + dy,
                    trapped.z + dz);
            }
        }
    }

    voxel_planner::PlannerConfig config;
    config.busbar_width = 3.0F;
    config.busbar_thickness = 1.0F;
    config.angle_step_deg = 15;
    const std::vector<voxel_planner::BusbarPose> poses =
        voxel_planner::generateDiscretePoses(config);
    const module2_morphology::CachedPoseFootprints cachedFootprints =
        module2_morphology::precomputePoseFootprints(
            poses,
            config);

    module2_morphology::populateVoxelPoseMasks(
        grid,
        poses,
        cachedFootprints);

    const Point3D wallAdjacent{2, 1, 1};
    const std::size_t wallAdjacentIndex = grid.index(
        wallAdjacent.x,
        wallAdjacent.y,
        wallAdjacent.z);
    const std::size_t trappedIndex = grid.index(
        trapped.x,
        trapped.y,
        trapped.z);
    const Point3D open{1, 3, 3};
    const std::size_t openIndex = grid.index(open.x, open.y, open.z);

    require(poses.size() == 624U,
            "The 15-degree pose table should contain 624 poses.");
    require(grid.poseMaskWordCount() == 10U,
            "The pose mask should span ten 64-bit words.");
    require(grid.getState(wallAdjacentIndex) ==
                VoxelState::POSE_CONDITIONAL,
            "A voxel adjacent to the wall should be pose-conditional.");
    require(grid.allowedPoseCount(wallAdjacentIndex) > 0U &&
                grid.allowedPoseCount(wallAdjacentIndex) < poses.size(),
            "The wall-adjacent mask should contain only some poses.");
    require(grid.getState(trappedIndex) == VoxelState::BLOCKED,
            "A completely trapped voxel should be blocked.");
    require(grid.allowedPoseCount(trappedIndex) == 0U,
            "A blocked voxel should have no allowed poses.");
    require(!grid.isPoseAllowed(trappedIndex, 623U),
            "A blocked voxel should reject a high poseId bit.");
    require(grid.getState(openIndex) == VoxelState::UNCONDITIONAL,
            "An open voxel should allow every pose.");
    require(grid.allowedPoseCount(openIndex) == poses.size(),
            "An unconditional voxel should expose every pose bit.");
    require(grid.isPoseAllowed(openIndex, 623U),
            "The final pose bit should be available in the tenth mask word.");

    std::cout << "[PASS] poses=" << poses.size()
              << " mask_words=" << grid.poseMaskWordCount()
              << " wall_allowed="
              << grid.allowedPoseCount(wallAdjacentIndex)
              << " trapped_allowed=" << grid.allowedPoseCount(trappedIndex)
              << " open_allowed=" << grid.allowedPoseCount(openIndex)
              << "\n";
    return 0;
}
