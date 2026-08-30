#include "Module2_Morphology/PoseMaskPopulation.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace module2_morphology {

CachedPoseFootprints precomputePoseFootprints(
    const std::vector<voxel_planner::BusbarPose>& poses,
    const voxel_planner::PlannerConfig& config) {
    CachedPoseFootprints cachedFootprints(poses.size());
    std::vector<std::uint8_t> seen(poses.size(), 0U);
    for (const voxel_planner::BusbarPose& pose : poses) {
        if (pose.poseId >= poses.size() || seen[pose.poseId] != 0U) {
            throw std::invalid_argument(
                "Pose identifiers must be unique and contiguous from zero.");
        }
        seen[pose.poseId] = 1U;
        cachedFootprints[pose.poseId] =
            voxel_planner::generatePoseFootprint(
                pose,
                config);
    }
    return cachedFootprints;
}

void populateVoxelPoseMasks(
    VoxelGrid& grid,
    const std::vector<voxel_planner::BusbarPose>& poses,
    const CachedPoseFootprints& cachedFootprints) {
    PoseMaskPopulator::populate(grid, poses, cachedFootprints);
}

void PoseMaskPopulator::populate(
    VoxelGrid& grid,
    const std::vector<voxel_planner::BusbarPose>& poses,
    const CachedPoseFootprints& cachedFootprints) {
    if (poses.empty() || cachedFootprints.size() != poses.size()) {
        throw std::invalid_argument(
            "Pose table and cached footprints must be non-empty and aligned.");
    }

    std::vector<std::uint8_t> seen(poses.size(), 0U);
    for (const voxel_planner::BusbarPose& pose : poses) {
        if (pose.poseId >= poses.size() || seen[pose.poseId] != 0U ||
            cachedFootprints[pose.poseId].empty()) {
            throw std::invalid_argument(
                "Pose identifiers and cached footprints are not aligned.");
        }
        seen[pose.poseId] = 1U;
    }

    grid.poseCount_ = poses.size();
    grid.poseMaskWordCount_ = (poses.size() + 63U) / 64U;
    grid.allPoseMask_.assign(grid.poseMaskWordCount_, ~std::uint64_t{0});
    const std::size_t remainder = poses.size() % 64U;
    if (remainder != 0U) {
        grid.allPoseMask_.back() =
            (std::uint64_t{1} << remainder) - 1U;
    }

    // No voxel collision checks occur here. The compact bit vectors are the
    // only grid-sized allocations; conditional masks are created on demand.
    grid.mask_evaluated_.assign(grid.storage_.size(), false);
    grid.all_poses_allowed_.assign(grid.storage_.size(), false);
    grid.poseMaskSlots_.clear();
    grid.allowedPoseMask_.clear();
    grid.lazyPoseFootprints_ =
        std::make_shared<const CachedPoseFootprints>(cachedFootprints);
    grid.lazyPoseFootprintBounds_.resize(cachedFootprints.size());
    voxel_planner::VoxelOffsetBounds unconditionalBounds{};
    bool unconditionalBoundsInitialized = false;
    for (std::size_t poseId = 0U;
         poseId < cachedFootprints.size();
         ++poseId) {
        const std::vector<voxel_planner::VoxelOffset>& footprint =
            cachedFootprints[poseId];
        voxel_planner::VoxelOffsetBounds bounds;
        bounds.minDx = bounds.maxDx = footprint.front().dx;
        bounds.minDy = bounds.maxDy = footprint.front().dy;
        bounds.minDz = bounds.maxDz = footprint.front().dz;
        for (const voxel_planner::VoxelOffset& offset : footprint) {
            bounds.minDx = std::min(bounds.minDx, offset.dx);
            bounds.maxDx = std::max(bounds.maxDx, offset.dx);
            bounds.minDy = std::min(bounds.minDy, offset.dy);
            bounds.maxDy = std::max(bounds.maxDy, offset.dy);
            bounds.minDz = std::min(bounds.minDz, offset.dz);
            bounds.maxDz = std::max(bounds.maxDz, offset.dz);
        }
        bounds.valid = true;
        grid.lazyPoseFootprintBounds_[poseId] = bounds;

        if (!unconditionalBoundsInitialized) {
            unconditionalBounds = bounds;
            unconditionalBoundsInitialized = true;
        } else {
            unconditionalBounds.minDx = std::min(
                unconditionalBounds.minDx,
                bounds.minDx);
            unconditionalBounds.maxDx = std::max(
                unconditionalBounds.maxDx,
                bounds.maxDx);
            unconditionalBounds.minDy = std::min(
                unconditionalBounds.minDy,
                bounds.minDy);
            unconditionalBounds.maxDy = std::max(
                unconditionalBounds.maxDy,
                bounds.maxDy);
            unconditionalBounds.minDz = std::min(
                unconditionalBounds.minDz,
                bounds.minDz);
            unconditionalBounds.maxDz = std::max(
                unconditionalBounds.maxDz,
                bounds.maxDz);
        }
    }
    unconditionalBounds.valid = unconditionalBoundsInitialized;
    grid.unconditionalPoseFootprintBounds_ = unconditionalBounds;
}

} // namespace module2_morphology
