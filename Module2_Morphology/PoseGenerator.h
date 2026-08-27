#pragma once

#include "Module1_Main/Types.h"

#include <vector>

namespace voxel_planner {

std::vector<BusbarPose> generateDiscretePoses(
    const PlannerConfig& config);

std::vector<VoxelOffset> generatePoseFootprint(
    const BusbarPose& pose,
    const PlannerConfig& config);

} // namespace voxel_planner
