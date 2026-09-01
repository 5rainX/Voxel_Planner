#pragma once

#include "Module1_Main/Types.h"
#include "Module2_Morphology/PoseGenerator.h"

#include <vector>

namespace module2_morphology {

using CachedPoseFootprints =
    std::vector<std::vector<voxel_planner::VoxelOffset>>;

/**
 * @brief Precomputes one local static-placement footprint per poseId.
 */
CachedPoseFootprints precomputePoseFootprints(
    const std::vector<voxel_planner::BusbarPose>& poses,
    const voxel_planner::PlannerConfig& config);

/**
 * @brief Precomputes and installs static pose masks for the complete grid.
 *
 * The returned VoxelGrid owns immutable-after-load mask storage. Queries
 * during search only read the packed masks and never allocate or mutate them.
 */
void populateVoxelPoseMasks(
    VoxelGrid& grid,
    const std::vector<voxel_planner::BusbarPose>& poses,
    const CachedPoseFootprints& cachedFootprints);

/**
 * @brief Trusted implementation with access to VoxelGrid's packed storage.
 */
class PoseMaskPopulator {
public:
    static void populate(
        VoxelGrid& grid,
        const std::vector<voxel_planner::BusbarPose>& poses,
        const CachedPoseFootprints& cachedFootprints);
};

} // namespace module2_morphology
