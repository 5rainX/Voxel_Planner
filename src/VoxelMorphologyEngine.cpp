#include "Module2_Morphology/VoxelMorphologyEngine.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace module2_morphology {
namespace {

std::size_t checkedVoxelCount(
    std::uint32_t width,
    std::uint32_t height,
    std::uint32_t depth) {
    const std::uint64_t count =
        static_cast<std::uint64_t>(width) * height * depth;
    if (count > std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error("Voxel grid is too large.");
    }
    return static_cast<std::size_t>(count);
}

std::size_t index(
    std::uint32_t x,
    std::uint32_t y,
    std::uint32_t z,
    std::uint32_t width,
    std::uint32_t height) noexcept {
    return static_cast<std::size_t>(x) +
           static_cast<std::size_t>(y) * width +
           static_cast<std::size_t>(z) * width * height;
}

std::vector<std::uint64_t> buildPrefixSum(
    const std::vector<std::uint8_t>& storage,
    std::uint32_t width,
    std::uint32_t height,
    std::uint32_t depth) {
    const std::uint32_t prefixWidth = width + 1U;
    const std::uint32_t prefixHeight = height + 1U;
    const std::uint32_t prefixDepth = depth + 1U;
    std::vector<std::uint64_t> prefix(
        checkedVoxelCount(prefixWidth, prefixHeight, prefixDepth),
        0U);

    for (std::uint32_t z = 1U; z < prefixDepth; ++z) {
        for (std::uint32_t y = 1U; y < prefixHeight; ++y) {
            for (std::uint32_t x = 1U; x < prefixWidth; ++x) {
                const std::uint64_t value =
                    storage[index(
                        x - 1U,
                        y - 1U,
                        z - 1U,
                        width,
                        height)] != 0U
                        ? 1U
                        : 0U;
                prefix[index(x, y, z, prefixWidth, prefixHeight)] =
                    value +
                    prefix[index(x - 1U, y, z, prefixWidth, prefixHeight)] +
                    prefix[index(x, y - 1U, z, prefixWidth, prefixHeight)] +
                    prefix[index(x, y, z - 1U, prefixWidth, prefixHeight)] -
                    prefix[index(
                        x - 1U,
                        y - 1U,
                        z,
                        prefixWidth,
                        prefixHeight)] -
                    prefix[index(
                        x - 1U,
                        y,
                        z - 1U,
                        prefixWidth,
                        prefixHeight)] -
                    prefix[index(
                        x,
                        y - 1U,
                        z - 1U,
                        prefixWidth,
                        prefixHeight)] +
                    prefix[index(
                        x - 1U,
                        y - 1U,
                        z - 1U,
                        prefixWidth,
                        prefixHeight)];
            }
        }
    }

    return prefix;
}

bool prefixDimensionsAreValid(
    const std::vector<std::uint64_t>& prefix,
    std::uint32_t prefixWidth,
    std::uint32_t prefixHeight,
    std::uint32_t prefixDepth) noexcept {
    const std::uint64_t expectedSize =
        static_cast<std::uint64_t>(prefixWidth) *
        prefixHeight *
        prefixDepth;
    return prefixWidth > 0U &&
           prefixHeight > 0U &&
           prefixDepth > 0U &&
           expectedSize <= std::numeric_limits<std::size_t>::max() &&
           prefix.size() == static_cast<std::size_t>(expectedSize);
}

std::uint32_t expandedDimension(std::uint32_t source, int pad) {
    const std::uint64_t expanded =
        static_cast<std::uint64_t>(source) +
        static_cast<std::uint64_t>(pad) * 2U;
    if (expanded > std::numeric_limits<std::uint32_t>::max()) {
        throw std::overflow_error("Expanded voxel dimension overflowed.");
    }
    return static_cast<std::uint32_t>(expanded);
}

} // namespace

MorphologyResult VoxelMorphologyEngine::buildCoarseBlockedMap(
    const VoxelGrid& rawGrid,
    MorphologyKernel kernel) {
    if (kernel.x <= 0 || kernel.y <= 0 || kernel.z <= 0) {
        throw std::invalid_argument("Morphology kernel must be positive.");
    }

    const int padX = kernel.x / 2;
    const int padY = kernel.y / 2;
    const int padZ = kernel.z / 2;
    const std::uint32_t expandedWidth =
        expandedDimension(rawGrid.width(), padX);
    const std::uint32_t expandedHeight =
        expandedDimension(rawGrid.height(), padY);
    const std::uint32_t expandedDepth =
        expandedDimension(rawGrid.depth(), padZ);

    std::vector<std::uint8_t> expanded(
        checkedVoxelCount(expandedWidth, expandedHeight, expandedDepth),
        0U);

    for (std::uint32_t z = 0U; z < rawGrid.depth(); ++z) {
        for (std::uint32_t y = 0U; y < rawGrid.height(); ++y) {
            for (std::uint32_t x = 0U; x < rawGrid.width(); ++x) {
                const std::size_t rawIndex = rawGrid.index(x, y, z);
                if (!rawGrid.isRawObstacle(rawIndex)) {
                    continue;
                }
                expanded[index(
                    x + static_cast<std::uint32_t>(padX),
                    y + static_cast<std::uint32_t>(padY),
                    z + static_cast<std::uint32_t>(padZ),
                    expandedWidth,
                    expandedHeight)] = VoxelGrid::kRawObstacleValue;
            }
        }
    }

    // The morphology engine preserves only raw obstacles in the padded map.
    // PoseMaskPopulator consumes this prefix sum immediately after loading
    // to precompute all static pose classifications.
    std::vector<std::uint64_t> rawPrefix = buildPrefixSum(
        expanded,
        expandedWidth,
        expandedHeight,
        expandedDepth);

    VoxelGrid result;
    result.replaceRawStorage(
        expandedWidth,
        expandedHeight,
        expandedDepth,
        std::move(expanded));

    return {
        std::move(result),
        voxel_planner::Point3D{padX, padY, padZ},
        std::move(rawPrefix),
        expandedWidth + 1U,
        expandedHeight + 1U,
        expandedDepth + 1U};
}

bool VoxelMorphologyEngine::isBoxCollisionFree(
    const std::vector<std::uint64_t>& prefix,
    std::uint32_t prefixWidth,
    std::uint32_t prefixHeight,
    std::uint32_t prefixDepth,
    int anchorX,
    int anchorY,
    int anchorZ,
    const voxel_planner::VoxelOffsetBounds& bounds) noexcept {
    if (!bounds.valid ||
        !prefixDimensionsAreValid(
            prefix,
            prefixWidth,
            prefixHeight,
            prefixDepth)) {
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
    const std::int64_t mapWidth =
        static_cast<std::int64_t>(prefixWidth) - 1;
    const std::int64_t mapHeight =
        static_cast<std::int64_t>(prefixHeight) - 1;
    const std::int64_t mapDepth =
        static_cast<std::int64_t>(prefixDepth) - 1;
    if (minX > maxX || minY > maxY || minZ > maxZ ||
        mapWidth <= 0 || mapHeight <= 0 || mapDepth <= 0) {
        return false;
    }

    const std::int64_t clampedMinX = std::max<std::int64_t>(0, minX);
    const std::int64_t clampedMinY = std::max<std::int64_t>(0, minY);
    const std::int64_t clampedMinZ = std::max<std::int64_t>(0, minZ);
    const std::int64_t clampedMaxX = std::min<std::int64_t>(
        mapWidth - 1,
        maxX);
    const std::int64_t clampedMaxY = std::min<std::int64_t>(
        mapHeight - 1,
        maxY);
    const std::int64_t clampedMaxZ = std::min<std::int64_t>(
        mapDepth - 1,
        maxZ);
    if (clampedMinX > clampedMaxX ||
        clampedMinY > clampedMaxY ||
        clampedMinZ > clampedMaxZ) {
        return false;
    }
    const bool trimmed = clampedMinX != minX || clampedMinY != minY ||
        clampedMinZ != minZ || clampedMaxX != maxX ||
        clampedMaxY != maxY || clampedMaxZ != maxZ;

    const std::uint32_t x1 = static_cast<std::uint32_t>(clampedMinX);
    const std::uint32_t x2 = static_cast<std::uint32_t>(clampedMaxX + 1);
    const std::uint32_t y1 = static_cast<std::uint32_t>(clampedMinY);
    const std::uint32_t y2 = static_cast<std::uint32_t>(clampedMaxY + 1);
    const std::uint32_t z1 = static_cast<std::uint32_t>(clampedMinZ);
    const std::uint32_t z2 = static_cast<std::uint32_t>(clampedMaxZ + 1);

    const auto at = [
        &prefix,
        prefixWidth,
        prefixHeight](std::uint32_t x,
                      std::uint32_t y,
                      std::uint32_t z) noexcept -> std::uint64_t {
        return prefix[
            static_cast<std::size_t>(x) +
            static_cast<std::size_t>(y) * prefixWidth +
            static_cast<std::size_t>(z) * prefixWidth * prefixHeight];
    };

    const std::uint64_t sum =
        at(x2, y2, z2) -
        at(x1, y2, z2) -
        at(x2, y1, z2) -
        at(x2, y2, z1) +
        at(x1, y1, z2) +
        at(x1, y2, z1) +
        at(x2, y1, z1) -
        at(x1, y1, z1);
    return !trimmed && sum == 0U;
}

} // namespace module2_morphology
