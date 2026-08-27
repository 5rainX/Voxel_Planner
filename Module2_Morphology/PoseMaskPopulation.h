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
 * @brief Installs the footprint table for lazy static-pose evaluation.
 *
 * This function does not traverse the grid. Each voxel is classified when
 * isPoseAllowed() or getState() first requests its mask.
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
