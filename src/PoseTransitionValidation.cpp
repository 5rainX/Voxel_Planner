#include "Module3_AStar/PoseTransitionValidation.h"

#include <cstdint>
#include <vector>

namespace module3_astar {
namespace {

ValidationResult failure(ErrorCode errorCode, const char* message) {
    return {errorCode, message};
}

} // namespace

ValidationResult validatePoseEndpoints(
    const VoxelGrid& map,
    const PoseState& startPose,
    const PoseState& endPose) {
    if (!map.isValid(
            startPose.position.x,
            startPose.position.y,
            startPose.position.z)) {
        return failure(
            ErrorCode::START_OUT_OF_BOUNDS,
            "Start pose is outside the voxel map.");
    }
    if (!map.isValid(
            endPose.position.x,
            endPose.position.y,
            endPose.position.z)) {
        return failure(
            ErrorCode::END_OUT_OF_BOUNDS,
            "End pose is outside the voxel map.");
    }
    if (map.poseCount() == 0U || map.poseMaskWordCount() == 0U) {
        return failure(
            ErrorCode::POSE_MASK_UNAVAILABLE,
            "Voxel pose masks have not been populated.");
    }
    if (startPose.poseId >= map.poseCount()) {
        return failure(
            ErrorCode::START_POSE_INVALID,
            "Start poseId is outside the generated pose table.");
    }
    if (endPose.poseId >= map.poseCount()) {
        return failure(
            ErrorCode::END_POSE_INVALID,
            "End poseId is outside the generated pose table.");
    }

    const std::size_t startIndex = map.index(
        static_cast<std::uint32_t>(startPose.position.x),
        static_cast<std::uint32_t>(startPose.position.y),
        static_cast<std::uint32_t>(startPose.position.z));
    if (!map.isPoseAllowed(startIndex, startPose.poseId)) {
        return failure(
            ErrorCode::START_POINT_BLOCKED,
            "Start pose collides with an obstacle or map boundary.");
    }

    const std::size_t endIndex = map.index(
        static_cast<std::uint32_t>(endPose.position.x),
        static_cast<std::uint32_t>(endPose.position.y),
        static_cast<std::uint32_t>(endPose.position.z));
    if (!map.isPoseAllowed(endIndex, endPose.poseId)) {
        return failure(
            ErrorCode::END_POINT_BLOCKED,
            "End pose collides with an obstacle or map boundary.");
    }
    return {};
}

ValidationResult validateTransitionSweep(
    const VoxelGrid& map,
    const Point3D& anchor,
    const std::vector<voxel_planner::VoxelOffset>& sweepOffsets) {
    return PoseTransitionValidator::validateSweep(
        map,
        anchor,
        sweepOffsets);
}

ValidationResult PoseTransitionValidator::validateSweep(
    const VoxelGrid& map,
    const Point3D& anchor,
    const std::vector<voxel_planner::VoxelOffset>& sweepOffsets) {
    if (sweepOffsets.empty()) {
        return failure(
            ErrorCode::INVALID_ARGUMENT,
            "Transition sweep must contain at least one voxel offset.");
    }

    for (const voxel_planner::VoxelOffset& offset : sweepOffsets) {
        const std::int64_t x =
            static_cast<std::int64_t>(anchor.x) + offset.dx;
        const std::int64_t y =
            static_cast<std::int64_t>(anchor.y) + offset.dy;
        const std::int64_t z =
            static_cast<std::int64_t>(anchor.z) + offset.dz;
        if (x < 0 || y < 0 || z < 0 ||
            x >= static_cast<std::int64_t>(map.width_) ||
            y >= static_cast<std::int64_t>(map.height_) ||
            z >= static_cast<std::int64_t>(map.depth_)) {
            return failure(
                ErrorCode::TRANSITION_OUT_OF_BOUNDS,
                "Transition swept volume leaves the voxel map.");
        }

        const std::size_t index = map.index(
            static_cast<std::uint32_t>(x),
            static_cast<std::uint32_t>(y),
            static_cast<std::uint32_t>(z));
        if (map.storage_[index] == VoxelGrid::kRawObstacleValue ||
            map.storage_[index] == VoxelGrid::kBlockedStorageValue) {
            return failure(
                ErrorCode::TRANSITION_COLLISION,
                "Transition swept volume intersects an obstacle or "
                "path keep-out zone.");
        }
    }
    return {};
}

} // namespace module3_astar
