#pragma once

#include "Module2_Morphology/PoseGenerator.h"

#include <vector>

namespace voxel_planner {

std::vector<VoxelOffset> generateStraightSweep(
    const BusbarPose& pose,
    const PlannerConfig& config,
    float distanceS);

std::vector<VoxelOffset> generateBendSweep(
    const BusbarPose& poseStart,
    const BusbarPose& poseEnd,
    const PlannerConfig& config,
    float bendRadius,
    float bendAngleRadians);

std::vector<VoxelOffset> generateExplicitBendSweep(
    const BusbarPose& poseStart,
    const PlannerConfig& config,
    float bendRadius,
    float bendAngleRadians,
    float axisX,
    float axisY,
    float axisZ);

std::vector<VoxelOffset> generateTwistSweep(
    const BusbarPose& poseStart,
    const BusbarPose& poseEnd,
    const PlannerConfig& config,
    float twistAngleRadians,
    float twistTravelLength);

} // namespace voxel_planner
