#pragma once

#include "VoxelPlannerAPI.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace module2_morphology {
class VoxelIO;
class VoxelMorphologyEngine;
class PoseMaskPopulator;
}
namespace module3_astar {
class CoarseAStar;
class PoseTransitionValidator;
}

using VoxelState = voxel_planner::VoxelClass;
using Point3D = voxel_planner::Point3D;
using Vector3D = voxel_planner::Vector3D;

struct Pose {
    Vector3D normal{};
    Vector3D tangent{};
    std::size_t waypointIndex = 0U;
};

struct PathResult {
    std::vector<Point3D> path;
    std::vector<Pose> pose_description;
    double cost = std::numeric_limits<double>::infinity();
    bool goal_tolerance_accepted = false;
};

enum class PlannerStatus : std::uint8_t {
    OK = 0,
    NO_PATH = 1
};

enum class ErrorCode : std::uint8_t {
    NONE = 0,
    INVALID_ARGUMENT,
    POSE_MASK_UNAVAILABLE,
    START_OUT_OF_BOUNDS,
    END_OUT_OF_BOUNDS,
    START_POSE_INVALID,
    END_POSE_INVALID,
    START_POINT_BLOCKED,
    END_POINT_BLOCKED,
    TRANSITION_OUT_OF_BOUNDS,
    TRANSITION_COLLISION,
    PATH_NOT_FOUND,
    MEMORY_ALLOCATION_FAILED,
    COMPUTATION_LIMIT_EXCEEDED
};

struct PlanningResult {
    PlannerStatus status = PlannerStatus::NO_PATH;
    std::vector<PathResult> paths;
    std::string message;
    ErrorCode error_code = ErrorCode::PATH_NOT_FOUND;
};

/**
 * @brief Position and discrete pose identifier for an SE(3) search endpoint.
 */
struct PoseState {
    PoseState() = default;
    PoseState(Point3D point, std::uint32_t id)
        : position(point), poseId(id) {}

    Point3D position{};
    std::uint32_t poseId = 0U;
};

namespace voxel_planner {

/** @brief Private engine configuration derived by the public facade. */
struct PlannerConfig {
    float busbar_width = 35.0F;
    float busbar_thickness = 5.0F;
    float twist_factor = 1.5F;
    float flat_bend_factor = 1.5F;
    float vertical_bend_factor = 0.557143F;
    int angle_step_deg = 15;
    float path_blocking_radius = 50.0F;
    float max_overlap_ratio = 0.3F;
};

/** @brief Private default terminal constraint used by the facade. */
struct EndpointConstraint {
    Vector3D start_normal{0.0, 0.0, 1.0};
    Vector3D start_tangent{1.0, 0.0, 0.0};
    Vector3D end_normal{0.0, 0.0, 1.0};
    Vector3D end_tangent{1.0, 0.0, 0.0};
    float min_begin_length = 0.0F;
    float min_end_length = 0.0F;
};

/** @brief Private discrete SE(3) state used only by implementation modules. */
struct BusbarPose {
    int tx = 0;
    int ty = 0;
    int tz = 0;
    float nx = 0.0F;
    float ny = 0.0F;
    float nz = 0.0F;
    std::uint32_t poseId = 0U;
    int roll_angle_deg = 0;
};

/** @brief Private integer offset in a rasterized swept volume. */
struct VoxelOffset {
    int dx = 0;
    int dy = 0;
    int dz = 0;
};

struct VoxelOffsetBounds {
    int minDx = 0;
    int maxDx = 0;
    int minDy = 0;
    int maxDy = 0;
    int minDz = 0;
    int maxDz = 0;
    bool valid = false;
};

} // namespace voxel_planner

/**
 * @brief Dense raw occupancy with sparse, precomputed pose masks.
 *
 * The storage is flattened with X as the fastest-changing dimension:
 * index = x + y * width + z * width * height.
 */
struct VoxelGrid {
    VoxelGrid() = default;

    VoxelGrid(
        std::uint32_t w,
        std::uint32_t h,
        std::uint32_t d)
        : width_(w),
          height_(h),
          depth_(d),
          storage_(
              static_cast<std::size_t>(w) * h * d,
              0U),
          occupancyBlockWidth_(
              (w + kOccupancyBlockSize - 1U) / kOccupancyBlockSize),
          occupancyBlockHeight_(
              (h + kOccupancyBlockSize - 1U) / kOccupancyBlockSize),
          occupancyBlockDepth_(
              (d + kOccupancyBlockSize - 1U) / kOccupancyBlockSize),
          occupancyBlocks_(
              static_cast<std::size_t>(occupancyBlockWidth_) *
                  occupancyBlockHeight_ * occupancyBlockDepth_,
              0U),
          all_poses_allowed_(storage_.size(), false) {}

    std::uint32_t width() const noexcept { return width_; }
    std::uint32_t height() const noexcept { return height_; }
    std::uint32_t depth() const noexcept { return depth_; }
    std::size_t voxelCount() const noexcept { return storage_.size(); }

    std::size_t poseCount() const noexcept { return poseCount_; }
    std::size_t poseMaskWordCount() const noexcept {
        return poseMaskWordCount_;
    }

    std::uint64_t getAllowedPoseMask(
        std::size_t linearIndex) const {
        return getAllowedPoseMaskWord(linearIndex, 0U);
    }

    std::uint64_t getAllowedPoseMaskWord(
        std::size_t linearIndex,
        std::size_t wordIndex) const {
        if (linearIndex >= storage_.size()) {
            throw std::out_of_range("Voxel index is outside the grid.");
        }
        if (wordIndex >= poseMaskWordCount_) {
            return 0U;
        }

        if (storage_[linearIndex] == kRawObstacleValue ||
            storage_[linearIndex] == kBlockedStorageValue) {
            return 0U;
        }
        if (all_poses_allowed_[linearIndex]) {
            return allPoseMask_[wordIndex];
        }

        const auto slotEntry = poseMaskSlots_.find(linearIndex);
        if (slotEntry == poseMaskSlots_.end()) {
            return 0U;
        }
        const std::uint32_t slot = slotEntry->second;
        return allowedPoseMask_[
            static_cast<std::size_t>(slot) * poseMaskWordCount_ +
            wordIndex];
    }

    bool isPoseAllowed(
        std::size_t linearIndex,
        std::uint32_t poseId) const {
        if (poseId >= poseCount_) {
            return false;
        }
        const std::size_t word = poseId / 64U;
        const std::size_t bit = poseId % 64U;
        return (getAllowedPoseMaskWord(linearIndex, word) &
                (std::uint64_t{1} << bit)) != 0U;
    }

    bool isObstacle(std::size_t linearIndex) const noexcept {
        return linearIndex >= storage_.size() ||
               storage_[linearIndex] == kRawObstacleValue ||
               storage_[linearIndex] == kBlockedStorageValue;
    }

    /**
     * @brief Tests only the source-map occupancy bit.
     *
     * Unlike isObstacle(), this deliberately ignores pose-derived storage
     * classifications so position-only planners can operate on raw input.
     */
    bool isRawObstacle(std::size_t linearIndex) const noexcept {
        return linearIndex >= storage_.size() ||
               storage_[linearIndex] == kRawObstacleValue;
    }

    bool isObstacleFreeBox(
        int anchorX,
        int anchorY,
        int anchorZ,
        const voxel_planner::VoxelOffsetBounds& bounds) const noexcept {
        if (!bounds.valid) {
            return false;
        }
        const std::int64_t minX =
            static_cast<std::int64_t>(anchorX) + bounds.minDx;
        const std::int64_t maxX =
            static_cast<std::int64_t>(anchorX) + bounds.maxDx;
        const std::int64_t minY =
            static_cast<std::int64_t>(anchorY) + bounds.minDy;
        const std::int64_t maxY =
            static_cast<std::int64_t>(anchorY) + bounds.maxDy;
        const std::int64_t minZ =
            static_cast<std::int64_t>(anchorZ) + bounds.minDz;
        const std::int64_t maxZ =
            static_cast<std::int64_t>(anchorZ) + bounds.maxDz;
        if (morphologyPrefix_ && prefixWidth_ > 0U &&
            prefixHeight_ > 0U && prefixDepth_ > 0U) {
            const std::int64_t mapWidth =
                static_cast<std::int64_t>(prefixWidth_) - 1;
            const std::int64_t mapHeight =
                static_cast<std::int64_t>(prefixHeight_) - 1;
            const std::int64_t mapDepth =
                static_cast<std::int64_t>(prefixDepth_) - 1;
            if (minX < 0 || minY < 0 || minZ < 0 ||
                maxX >= mapWidth || maxY >= mapHeight ||
                maxZ >= mapDepth ||
                minX > maxX || minY > maxY || minZ > maxZ) {
                return false;
            }
            const auto at = [
                this](std::uint32_t x,
                      std::uint32_t y,
                      std::uint32_t z) noexcept -> std::uint64_t {
                return (*morphologyPrefix_)[
                    static_cast<std::size_t>(x) +
                    static_cast<std::size_t>(y) * prefixWidth_ +
                    static_cast<std::size_t>(z) * prefixWidth_ *
                        prefixHeight_];
            };
            const std::uint32_t x1 = static_cast<std::uint32_t>(minX);
            const std::uint32_t x2 = static_cast<std::uint32_t>(maxX + 1);
            const std::uint32_t y1 = static_cast<std::uint32_t>(minY);
            const std::uint32_t y2 = static_cast<std::uint32_t>(maxY + 1);
            const std::uint32_t z1 = static_cast<std::uint32_t>(minZ);
            const std::uint32_t z2 = static_cast<std::uint32_t>(maxZ + 1);
            const std::uint64_t sum =
                at(x2, y2, z2) -
                at(x1, y2, z2) -
                at(x2, y1, z2) -
                at(x2, y2, z1) +
                at(x1, y1, z2) +
                at(x1, y2, z1) +
                at(x2, y1, z1) -
                at(x1, y1, z1);
            return sum == 0U;
        }
        if (minX < 0 || minY < 0 || minZ < 0 ||
            maxX >= width_ || maxY >= height_ || maxZ >= depth_) {
            return false;
        }

        const std::uint32_t minBlockX =
            static_cast<std::uint32_t>(minX) / kOccupancyBlockSize;
        const std::uint32_t maxBlockX =
            static_cast<std::uint32_t>(maxX) / kOccupancyBlockSize;
        const std::uint32_t minBlockY =
            static_cast<std::uint32_t>(minY) / kOccupancyBlockSize;
        const std::uint32_t maxBlockY =
            static_cast<std::uint32_t>(maxY) / kOccupancyBlockSize;
        const std::uint32_t minBlockZ =
            static_cast<std::uint32_t>(minZ) / kOccupancyBlockSize;
        const std::uint32_t maxBlockZ =
            static_cast<std::uint32_t>(maxZ) / kOccupancyBlockSize;
        for (std::uint32_t bz = minBlockZ; bz <= maxBlockZ; ++bz) {
            for (std::uint32_t by = minBlockY; by <= maxBlockY; ++by) {
                for (std::uint32_t bx = minBlockX; bx <= maxBlockX; ++bx) {
                    if (occupancyBlocks_[occupancyBlockIndex(bx, by, bz)] !=
                        0U) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     * Read-only single-pose fast path used by parallel A* expansion.
     *
     * Pose masks are populated before a ProcessedMap becomes visible to
     * search, so this query performs no cache allocation or mutation.
     */
    bool isPoseAllowedForSearch(
        std::size_t linearIndex,
        std::uint32_t poseId) const {
        return isPoseAllowed(linearIndex, poseId);
    }

    std::size_t allowedPoseCount(std::size_t linearIndex) const {
        std::size_t count = 0U;
        for (std::size_t word = 0U;
             word < poseMaskWordCount_;
             ++word) {
            std::uint64_t bits = getAllowedPoseMaskWord(
                linearIndex,
                word);
            while (bits != 0U) {
                bits &= bits - 1U;
                ++count;
            }
        }
        return count;
    }

    void setRawObstacle(int x, int y, int z, bool obstacle = true) {
        if (!isValid(x, y, z)) {
            throw std::out_of_range("Voxel coordinate is outside the grid.");
        }
        if (poseCount_ != 0U) {
            throw std::logic_error(
                "Raw obstacles cannot change after pose-mask initialization.");
        }
        const std::size_t linearIndex = index(
            static_cast<std::uint32_t>(x),
            static_cast<std::uint32_t>(y),
            static_cast<std::uint32_t>(z));
        storage_[linearIndex] = obstacle
                ? kRawObstacleValue
                : 0U;
        if (obstacle) {
            markOccupancyBlock(
                static_cast<std::uint32_t>(x),
                static_cast<std::uint32_t>(y),
                static_cast<std::uint32_t>(z));
        } else {
            rebuildOccupancyBlock(
                static_cast<std::uint32_t>(x) / kOccupancyBlockSize,
                static_cast<std::uint32_t>(y) / kOccupancyBlockSize,
                static_cast<std::uint32_t>(z) / kOccupancyBlockSize);
        }
        morphologyPrefix_.reset();
        prefixWidth_ = 0U;
        prefixHeight_ = 0U;
        prefixDepth_ = 0U;
    }

    void attachMorphologyPrefix(
        std::shared_ptr<const std::vector<std::uint64_t>> prefix,
        std::uint32_t prefixWidth,
        std::uint32_t prefixHeight,
        std::uint32_t prefixDepth) {
        if (!prefix) {
            morphologyPrefix_.reset();
            prefixWidth_ = 0U;
            prefixHeight_ = 0U;
            prefixDepth_ = 0U;
            return;
        }
        if (prefixWidth != width_ + 1U ||
            prefixHeight != height_ + 1U ||
            prefixDepth != depth_ + 1U ||
            prefix->size() != static_cast<std::size_t>(prefixWidth) *
                                 prefixHeight * prefixDepth) {
            throw std::invalid_argument(
                "Morphology prefix dimensions do not match its storage.");
        }
        morphologyPrefix_ = std::move(prefix);
        prefixWidth_ = prefixWidth;
        prefixHeight_ = prefixHeight;
        prefixDepth_ = prefixDepth;
    }

    void replaceRawStorage(
        std::uint32_t w,
        std::uint32_t h,
        std::uint32_t d,
        std::vector<std::uint8_t> storage) {
        if (poseCount_ != 0U) {
            throw std::logic_error(
                "Raw storage cannot be replaced after pose-mask initialization.");
        }
        if (storage.size() != static_cast<std::size_t>(w) * h * d) {
            throw std::invalid_argument("Raw voxel storage size mismatch.");
        }

        width_ = w;
        height_ = h;
        depth_ = d;
        storage_ = std::move(storage);

        occupancyBlockWidth_ =
            (w + kOccupancyBlockSize - 1U) / kOccupancyBlockSize;
        occupancyBlockHeight_ =
            (h + kOccupancyBlockSize - 1U) / kOccupancyBlockSize;
        occupancyBlockDepth_ =
            (d + kOccupancyBlockSize - 1U) / kOccupancyBlockSize;
        occupancyBlocks_.assign(
            static_cast<std::size_t>(occupancyBlockWidth_) *
                occupancyBlockHeight_ * occupancyBlockDepth_,
            0U);

        all_poses_allowed_.assign(storage_.size(), false);
        poseMaskSlots_.clear();
        allowedPoseMask_.clear();
        allPoseMask_.clear();
        poseMaskWordCount_ = 0U;
        morphologyPrefix_.reset();
        prefixWidth_ = 0U;
        prefixHeight_ = 0U;
        prefixDepth_ = 0U;
        rebuildOccupancyBlocks();
    }

    std::size_t index(
        std::uint32_t x,
        std::uint32_t y,
        std::uint32_t z) const noexcept {
        return static_cast<std::size_t>(x) +
               static_cast<std::size_t>(y) * width_ +
               static_cast<std::size_t>(z) * width_ * height_;
    }

    bool isValid(int x, int y, int z) const noexcept {
        return x >= 0 &&
               x < static_cast<int>(width_) &&
               y >= 0 &&
               y < static_cast<int>(height_) &&
               z >= 0 &&
               z < static_cast<int>(depth_);
    }

    VoxelState getState(int x, int y, int z) const {
        if (!isValid(x, y, z)) {
            throw std::out_of_range("Voxel coordinate is outside the grid.");
        }
        return getState(index(
            static_cast<std::uint32_t>(x),
            static_cast<std::uint32_t>(y),
            static_cast<std::uint32_t>(z)));
    }

    VoxelState getState(std::size_t linearIndex) const {
        if (linearIndex >= storage_.size()) {
            throw std::out_of_range("Voxel index is outside the grid.");
        }
        const std::uint8_t value = storage_[linearIndex];
        if (value == kBlockedStorageValue ||
            value == kRawObstacleValue) {
            return VoxelState::BLOCKED;
        }
        if (poseCount_ == 0U) {
            return VoxelState::UNCONDITIONAL;
        }

        if (all_poses_allowed_[linearIndex]) {
            return VoxelState::UNCONDITIONAL;
        }
        return poseMaskSlots_.find(linearIndex) != poseMaskSlots_.end()
            ? VoxelState::POSE_CONDITIONAL
            : VoxelState::BLOCKED;
    }

private:
    std::uint32_t width_ = 0;
    std::uint32_t height_ = 0;
    std::uint32_t depth_ = 0;
    std::vector<std::uint8_t> storage_;

    static constexpr std::uint8_t kRawObstacleValue = 255U;
    static constexpr std::uint8_t kBlockedStorageValue = 3U;
    static constexpr std::uint32_t kOccupancyBlockSize = 8U;
    std::uint32_t occupancyBlockWidth_ = 0U;
    std::uint32_t occupancyBlockHeight_ = 0U;
    std::uint32_t occupancyBlockDepth_ = 0U;
    std::vector<std::uint8_t> occupancyBlocks_;
    std::size_t poseCount_ = 0U;
    std::size_t poseMaskWordCount_ = 0U;
    std::vector<std::uint64_t> allPoseMask_;
    std::vector<bool> all_poses_allowed_;
    std::unordered_map<std::size_t, std::uint32_t> poseMaskSlots_;
    std::vector<std::uint64_t> allowedPoseMask_;
    std::shared_ptr<const std::vector<std::uint64_t>> morphologyPrefix_;
    std::uint32_t prefixWidth_ = 0U;
    std::uint32_t prefixHeight_ = 0U;
    std::uint32_t prefixDepth_ = 0U;

    std::size_t occupancyBlockIndex(
        std::uint32_t bx,
        std::uint32_t by,
        std::uint32_t bz) const noexcept {
        return static_cast<std::size_t>(bx) +
            static_cast<std::size_t>(by) * occupancyBlockWidth_ +
            static_cast<std::size_t>(bz) * occupancyBlockWidth_ *
                occupancyBlockHeight_;
    }

    void markOccupancyBlock(
        std::uint32_t x,
        std::uint32_t y,
        std::uint32_t z) noexcept {
        occupancyBlocks_[occupancyBlockIndex(
            x / kOccupancyBlockSize,
            y / kOccupancyBlockSize,
            z / kOccupancyBlockSize)] = 1U;
    }

    void rebuildOccupancyBlock(
        std::uint32_t bx,
        std::uint32_t by,
        std::uint32_t bz) noexcept {
        const std::uint32_t beginX = bx * kOccupancyBlockSize;
        const std::uint32_t beginY = by * kOccupancyBlockSize;
        const std::uint32_t beginZ = bz * kOccupancyBlockSize;
        const std::uint32_t endX = std::min(
            width_, beginX + kOccupancyBlockSize);
        const std::uint32_t endY = std::min(
            height_, beginY + kOccupancyBlockSize);
        const std::uint32_t endZ = std::min(
            depth_, beginZ + kOccupancyBlockSize);
        std::uint8_t occupied = 0U;
        for (std::uint32_t z = beginZ; z < endZ && occupied == 0U; ++z) {
            for (std::uint32_t y = beginY; y < endY && occupied == 0U; ++y) {
                for (std::uint32_t x = beginX; x < endX; ++x) {
                    if (isObstacle(index(x, y, z))) {
                        occupied = 1U;
                        break;
                    }
                }
            }
        }
        occupancyBlocks_[occupancyBlockIndex(bx, by, bz)] = occupied;
    }

    void rebuildOccupancyBlocks() noexcept {
        const std::int64_t blockCount = static_cast<std::int64_t>(
            occupancyBlocks_.size());
#pragma omp parallel for schedule(static) if(blockCount > 1024)
        for (std::int64_t block = 0; block < blockCount; ++block) {
            const std::uint32_t bx = static_cast<std::uint32_t>(
                block % occupancyBlockWidth_);
            const std::uint32_t by = static_cast<std::uint32_t>(
                (block / occupancyBlockWidth_) % occupancyBlockHeight_);
            const std::uint32_t bz = static_cast<std::uint32_t>(
                block /
                (static_cast<std::int64_t>(occupancyBlockWidth_) *
                 occupancyBlockHeight_));
            rebuildOccupancyBlock(bx, by, bz);
        }
    }

    friend class module2_morphology::VoxelIO;
    friend class module2_morphology::VoxelMorphologyEngine;
    friend class module2_morphology::PoseMaskPopulator;
    friend class module3_astar::CoarseAStar;
    friend class module3_astar::PoseTransitionValidator;
};
