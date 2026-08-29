#pragma once

#include "Module1_Main/Types.h"

#include <cstdint>
#include <vector>

namespace module2_morphology {

struct MorphologyKernel {
    int x = 1;
    int y = 1;
    int z = 1;
};

struct MorphologyResult {
    VoxelGrid grid;
    voxel_planner::Point3D originOffset{};
    std::vector<std::uint64_t> prefixSum;
    std::uint32_t prefixWidth = 0U;
    std::uint32_t prefixHeight = 0U;
    std::uint32_t prefixDepth = 0U;
};

class VoxelMorphologyEngine {
public:
    static MorphologyResult buildCoarseBlockedMap(
        const VoxelGrid& rawGrid,
        MorphologyKernel kernel);

    /**
     * @brief Tests an inclusive offset box against a 3D prefix-sum field.
     *
     * The query is O(1), uses the eight inclusion-exclusion corners, and
     * returns false when the requested box is invalid or outside the field.
     */
    static bool isBoxCollisionFree(
        const std::vector<std::uint64_t>& prefix,
        std::uint32_t prefixWidth,
        std::uint32_t prefixHeight,
        std::uint32_t prefixDepth,
        int anchorX,
        int anchorY,
        int anchorZ,
        const voxel_planner::VoxelOffsetBounds& bounds) noexcept;
};

} // namespace module2_morphology
