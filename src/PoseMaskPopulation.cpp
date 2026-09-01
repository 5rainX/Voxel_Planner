#include "Module2_Morphology/PoseMaskPopulation.h"

#include <array>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

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
    if (grid.poseMaskWordCount_ > 10U) {
        throw std::overflow_error(
            "Pose table exceeds the fixed local mask capacity.");
    }
    grid.allPoseMask_.assign(grid.poseMaskWordCount_, ~std::uint64_t{0});
    const std::size_t remainder = poses.size() % 64U;
    if (remainder != 0U) {
        grid.allPoseMask_.back() =
            (std::uint64_t{1} << remainder) - 1U;
    }

    if (grid.poseMaskWordCount_ != 0U &&
        grid.storage_.size() >
            std::numeric_limits<std::size_t>::max() /
                grid.poseMaskWordCount_) {
        throw std::overflow_error(
            "Expanded conditional pose-mask storage would overflow.");
    }

    std::vector<std::uint8_t> allAllowedFlags(
        grid.storage_.size(),
        0U);
    std::vector<voxel_planner::VoxelOffsetBounds> poseBounds(
        cachedFootprints.size());
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
        poseBounds[poseId] = bounds;

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

    struct ConditionalMaskRecord {
        std::size_t linearIndex = 0U;
        std::array<std::uint64_t, 10U> mask{};
    };

    std::vector<ConditionalMaskRecord> conditionalRecords;

    // Build every voxel's mask while the map is still exclusively owned by
    // the loading thread. Search threads subsequently access only these
    // completed vectors and the sparse conditional-mask table.
    const std::size_t plane =
        static_cast<std::size_t>(grid.width_) * grid.height_;
#ifdef _OPENMP
    const int maxThreads = omp_get_max_threads();
    std::vector<std::vector<ConditionalMaskRecord>> threadRecords(
        static_cast<std::size_t>(maxThreads));
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        std::vector<ConditionalMaskRecord>& localRecords =
            threadRecords[static_cast<std::size_t>(threadId)];
        localRecords.reserve(1024U);
#pragma omp for schedule(dynamic, 1024)
        for (std::size_t linearIndex = 0U;
             linearIndex < grid.storage_.size();
             ++linearIndex) {
            if (grid.isObstacle(linearIndex)) {
                continue;
            }

            const int anchorZ = static_cast<int>(linearIndex / plane);
            const std::size_t planeIndex = linearIndex % plane;
            const int anchorY = static_cast<int>(
                planeIndex / grid.width_);
            const int anchorX = static_cast<int>(
                planeIndex % grid.width_);

            if (grid.isObstacleFreeBox(
                    anchorX,
                    anchorY,
                    anchorZ,
                    unconditionalBounds)) {
                allAllowedFlags[linearIndex] = 1U;
                continue;
            }

            std::array<std::uint64_t, 10U> localMask{};
            std::size_t validPoseCount = 0U;
            for (std::size_t poseId = 0U;
                 poseId < cachedFootprints.size();
                 ++poseId) {
                bool collides = false;
                if (grid.isObstacleFreeBox(
                        anchorX,
                        anchorY,
                        anchorZ,
                        poseBounds[poseId])) {
                    localMask[poseId / 64U] |=
                        std::uint64_t{1} << (poseId % 64U);
                    ++validPoseCount;
                    continue;
                }

                for (const voxel_planner::VoxelOffset& offset :
                     cachedFootprints[poseId]) {
                    const int x = anchorX + offset.dx;
                    const int y = anchorY + offset.dy;
                    const int z = anchorZ + offset.dz;
                    if (!grid.isValid(x, y, z) ||
                        grid.isObstacle(grid.index(
                            static_cast<std::uint32_t>(x),
                            static_cast<std::uint32_t>(y),
                            static_cast<std::uint32_t>(z)))) {
                        collides = true;
                        break;
                    }
                }
                if (!collides) {
                    localMask[poseId / 64U] |=
                        std::uint64_t{1} << (poseId % 64U);
                    ++validPoseCount;
                }
            }

            if (validPoseCount == poses.size()) {
                allAllowedFlags[linearIndex] = 1U;
            } else if (validPoseCount != 0U) {
                ConditionalMaskRecord record;
                record.linearIndex = linearIndex;
                record.mask = localMask;
                localRecords.push_back(record);
            }
        }
    }
    for (std::vector<ConditionalMaskRecord>& localRecords : threadRecords) {
        conditionalRecords.insert(
            conditionalRecords.end(),
            localRecords.begin(),
            localRecords.end());
    }
#else
    for (std::size_t linearIndex = 0U;
         linearIndex < grid.storage_.size();
         ++linearIndex) {
        if (grid.isObstacle(linearIndex)) {
            continue;
        }

        const int anchorZ = static_cast<int>(linearIndex / plane);
        const std::size_t planeIndex = linearIndex % plane;
        const int anchorY = static_cast<int>(
            planeIndex / grid.width_);
        const int anchorX = static_cast<int>(
            planeIndex % grid.width_);

        if (grid.isObstacleFreeBox(
                anchorX,
                anchorY,
                anchorZ,
                unconditionalBounds)) {
            allAllowedFlags[linearIndex] = 1U;
            continue;
        }

        std::array<std::uint64_t, 10U> localMask{};
        std::size_t validPoseCount = 0U;
        for (std::size_t poseId = 0U;
             poseId < cachedFootprints.size();
             ++poseId) {
            bool collides = false;
            if (grid.isObstacleFreeBox(
                    anchorX,
                    anchorY,
                    anchorZ,
                    poseBounds[poseId])) {
                localMask[poseId / 64U] |=
                    std::uint64_t{1} << (poseId % 64U);
                ++validPoseCount;
                continue;
            }

            for (const voxel_planner::VoxelOffset& offset :
                 cachedFootprints[poseId]) {
                const int x = anchorX + offset.dx;
                const int y = anchorY + offset.dy;
                const int z = anchorZ + offset.dz;
                if (!grid.isValid(x, y, z) ||
                    grid.isObstacle(grid.index(
                        static_cast<std::uint32_t>(x),
                        static_cast<std::uint32_t>(y),
                        static_cast<std::uint32_t>(z)))) {
                    collides = true;
                    break;
                }
            }
            if (!collides) {
                localMask[poseId / 64U] |=
                    std::uint64_t{1} << (poseId % 64U);
                ++validPoseCount;
            }
        }

        if (validPoseCount == poses.size()) {
            allAllowedFlags[linearIndex] = 1U;
        } else if (validPoseCount != 0U) {
            ConditionalMaskRecord record;
            record.linearIndex = linearIndex;
            record.mask = localMask;
            conditionalRecords.push_back(record);
        }
    }
#endif

    grid.all_poses_allowed_.assign(
        allAllowedFlags.begin(),
        allAllowedFlags.end());
    grid.poseMaskSlots_.clear();
    grid.allowedPoseMask_.clear();
    if (!conditionalRecords.empty()) {
        std::sort(
            conditionalRecords.begin(),
            conditionalRecords.end(),
            [](const ConditionalMaskRecord& lhs,
               const ConditionalMaskRecord& rhs) {
                return lhs.linearIndex < rhs.linearIndex;
            });
        grid.allowedPoseMask_.reserve(
            conditionalRecords.size() * grid.poseMaskWordCount_);
        for (const ConditionalMaskRecord& record : conditionalRecords) {
            const std::uint32_t slot = static_cast<std::uint32_t>(
                grid.allowedPoseMask_.size() / grid.poseMaskWordCount_);
            grid.poseMaskSlots_.emplace(record.linearIndex, slot);
            const std::size_t base = grid.allowedPoseMask_.size();
            grid.allowedPoseMask_.resize(
                base + grid.poseMaskWordCount_,
                0U);
            for (std::size_t word = 0U;
                 word < grid.poseMaskWordCount_;
                 ++word) {
                grid.allowedPoseMask_[base + word] = record.mask[word];
            }
        }
    }
}

} // namespace module2_morphology
