#pragma once

#include "Module1_Main/Types.h"
#include "Module2_Morphology/PoseGenerator.h"

#include <string>
#include <vector>

namespace module3_astar {

/**
 * @brief Non-throwing result used before SE(3) search and edge expansion.
 */
struct ValidationResult {
    ErrorCode error_code = ErrorCode::NONE;
    std::string message;

    bool ok() const noexcept { return error_code == ErrorCode::NONE; }
};

/**
 * @brief Validates endpoint bounds, poseId range, and static pose masks.
 */
ValidationResult validatePoseEndpoints(
    const VoxelGrid& map,
    const PoseState& startPose,
    const PoseState& endPose);

/**
 * @brief Rejects any swept-volume offset that leaves the map or hits a raw
 * obstacle. Bounds are checked before linear indexing.
 */
ValidationResult validateTransitionSweep(
    const VoxelGrid& map,
    const Point3D& anchor,
    const std::vector<voxel_planner::VoxelOffset>& sweepOffsets);

/**
 * @brief Trusted implementation with read-only access to raw occupancy.
 */
class PoseTransitionValidator {
public:
    static ValidationResult validateSweep(
        const VoxelGrid& map,
        const Point3D& anchor,
        const std::vector<voxel_planner::VoxelOffset>& sweepOffsets);
};

} // namespace module3_astar
