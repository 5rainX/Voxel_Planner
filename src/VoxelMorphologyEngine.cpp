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

std::uint64_t boxSum(
    const std::vector<std::uint64_t>& prefix,
    int x1,
    int x2,
    int y1,
    int y2,
    int z1,
    int z2,
    std::uint32_t prefixWidth,
    std::uint32_t prefixHeight) {
    const std::uint32_t minX = static_cast<std::uint32_t>(x1);
    const std::uint32_t maxX = static_cast<std::uint32_t>(x2 + 1);
    const std::uint32_t minY = static_cast<std::uint32_t>(y1);
    const std::uint32_t maxY = static_cast<std::uint32_t>(y2 + 1);
    const std::uint32_t minZ = static_cast<std::uint32_t>(z1);
    const std::uint32_t maxZ = static_cast<std::uint32_t>(z2 + 1);

    return prefix[index(maxX, maxY, maxZ, prefixWidth, prefixHeight)] -
           prefix[index(minX, maxY, maxZ, prefixWidth, prefixHeight)] -
           prefix[index(maxX, minY, maxZ, prefixWidth, prefixHeight)] -
           prefix[index(maxX, maxY, minZ, prefixWidth, prefixHeight)] +
           prefix[index(minX, minY, maxZ, prefixWidth, prefixHeight)] +
           prefix[index(minX, maxY, minZ, prefixWidth, prefixHeight)] +
           prefix[index(maxX, minY, minZ, prefixWidth, prefixHeight)] -
           prefix[index(minX, minY, minZ, prefixWidth, prefixHeight)];
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

struct KernelAxisBounds {
    int minOffset = 0;
    int maxOffset = 0;
};

KernelAxisBounds kernelAxisBounds(int kernelSize) noexcept {
    // For an even-sized discrete kernel, use exactly kernelSize cells:
    // [-kernelSize / 2, kernelSize - kernelSize / 2 - 1].
    const int minOffset = kernelSize / 2;
    return {
        -minOffset,
        kernelSize - minOffset - 1};
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
    const KernelAxisBounds xBounds = kernelAxisBounds(kernel.x);
    const KernelAxisBounds yBounds = kernelAxisBounds(kernel.y);
    const KernelAxisBounds zBounds = kernelAxisBounds(kernel.z);
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

    const std::vector<std::uint64_t> expandedPrefix = buildPrefixSum(
        expanded,
        expandedWidth,
        expandedHeight,
        expandedDepth);

    std::vector<std::uint8_t> coarse = expanded;
    for (std::uint32_t z = 0U; z < expandedDepth; ++z) {
        for (std::uint32_t y = 0U; y < expandedHeight; ++y) {
            for (std::uint32_t x = 0U; x < expandedWidth; ++x) {
                const std::size_t linear =
                    index(x, y, z, expandedWidth, expandedHeight);
                if (coarse[linear] == VoxelGrid::kRawObstacleValue) {
                    continue;
                }

                const int x1 = std::max<int>(
                    0,
                    static_cast<int>(x) + xBounds.minOffset);
                const int x2 = std::min<int>(
                    static_cast<int>(expandedWidth) - 1,
                    static_cast<int>(x) + xBounds.maxOffset);
                const int y1 = std::max<int>(
                    0,
                    static_cast<int>(y) + yBounds.minOffset);
                const int y2 = std::min<int>(
                    static_cast<int>(expandedHeight) - 1,
                    static_cast<int>(y) + yBounds.maxOffset);
                const int z1 = std::max<int>(
                    0,
                    static_cast<int>(z) + zBounds.minOffset);
                const int z2 = std::min<int>(
                    static_cast<int>(expandedDepth) - 1,
                    static_cast<int>(z) + zBounds.maxOffset);

                if (boxSum(
                        expandedPrefix,
                        x1,
                        x2,
                        y1,
                        y2,
                        z1,
                        z2,
                        expandedWidth + 1U,
                        expandedHeight + 1U) > 0U) {
                    coarse[linear] = VoxelGrid::kRawObstacleValue;
                }
            }
        }
    }

    // Dynamic queries must describe the final coarse map, including cells
    // newly blocked by the morphology pass.
    std::vector<std::uint64_t> coarsePrefix = buildPrefixSum(
        coarse,
        expandedWidth,
        expandedHeight,
        expandedDepth);

    VoxelGrid result;
    result.replaceRawStorage(
        expandedWidth,
        expandedHeight,
        expandedDepth,
        std::move(coarse));

    return {
        std::move(result),
        voxel_planner::Point3D{padX, padY, padZ},
        std::move(coarsePrefix),
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
    if (minX < 0 || minY < 0 || minZ < 0 ||
        maxX >= mapWidth || maxY >= mapHeight || maxZ >= mapDepth ||
        minX > maxX || minY > maxY || minZ > maxZ) {
        return false;
    }

    const std::uint32_t x1 = static_cast<std::uint32_t>(minX);
    const std::uint32_t x2 = static_cast<std::uint32_t>(maxX + 1);
    const std::uint32_t y1 = static_cast<std::uint32_t>(minY);
    const std::uint32_t y2 = static_cast<std::uint32_t>(maxY + 1);
    const std::uint32_t z1 = static_cast<std::uint32_t>(minZ);
    const std::uint32_t z2 = static_cast<std::uint32_t>(maxZ + 1);

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
    return sum == 0U;
}

} // namespace module2_morphology
