#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

struct VoxelGrid;

namespace environment {

/**
 * @brief Three-dimensional physical-space vector.
 *
 * VoxelGrid currently uses an implicit physical frame: origin (0, 0, 0)
 * and unit voxel spacing. GlobalSdf stores that frame explicitly so the
 * continuous representation can later accept physical map metadata without
 * changing the existing VoxelGrid contract.
 */
struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/**
 * @brief Signed Euclidean distance field built from a VoxelGrid.
 *
 * Free-space samples are positive, occupied samples are negative. The
 * distance samples are stored in a flattened X-fastest array:
 *
 *   index = z * width * height + y * width + x
 *
 * The constructor taking only a VoxelGrid uses the coordinate convention
 * already implied by VoxelGrid: origin (0, 0, 0) and voxel size (1, 1, 1).
 */
class GlobalSdf {
public:
    explicit GlobalSdf(const VoxelGrid& grid);

    /**
     * @brief Builds an SDF with explicit physical frame metadata.
     *
     * This overload does not alter VoxelGrid. It is provided for callers
     * whose external map format carries non-unit spacing or a non-zero
     * physical origin.
     */
    GlobalSdf(
        const VoxelGrid& grid,
        Vec3 origin,
        Vec3 voxel_size);

    const Vec3& origin() const noexcept;
    const Vec3& voxelSize() const noexcept;
    const std::array<std::size_t, 3>& dimensions() const noexcept;

    /**
     * @brief Returns a trilinearly interpolated signed distance in physical
     * coordinates.
     *
     * Positions outside the sampled domain return a large negative soft-wall
     * distance so optimizers can reject or recover without exception-based
     * control flow.
     */
    double getDistance(const Vec3& pos) const noexcept;

    /**
     * @brief Returns the physical-space SDF gradient using finite differences.
     */
    Vec3 getGradient(const Vec3& pos) const noexcept;

    /**
     * @brief Returns true when a physical-space query lies inside the sampled
     * SDF domain.
     */
    bool contains(const Vec3& pos) const noexcept;

    /**
     * @brief Provides read-only access to the flattened float field.
     */
    const std::vector<float>& data() const noexcept;

private:
    static constexpr float kInfinity = 1.0e30F;

    std::size_t index(
        std::size_t x,
        std::size_t y,
        std::size_t z) const noexcept;

    static bool finiteVec(const Vec3& value) noexcept;
    static void validateVoxelSize(const Vec3& voxel_size);
    static std::size_t checkedVoxelCount(
        const std::array<std::size_t, 3>& dimensions);

    std::vector<float> buildSquaredDistanceField(
        const VoxelGrid& grid,
        bool obstacleIsFeature) const;

    float sample(
        std::size_t x,
        std::size_t y,
        std::size_t z) const noexcept;

    float sampleClamped(
        std::ptrdiff_t x,
        std::ptrdiff_t y,
        std::ptrdiff_t z) const noexcept;

    Vec3 toGridCoordinate(const Vec3& pos) const noexcept;
    Vec3 domainCenter() const noexcept;
    Vec3 directionToCenter(const Vec3& pos) const noexcept;
    double surfaceOffset() const noexcept;

    Vec3 origin_{};
    Vec3 voxel_size_{1.0, 1.0, 1.0};
    std::array<std::size_t, 3> dimensions_{0U, 0U, 0U};
    std::vector<float> distances_;
};

/**
 * @brief Reserved multi-scale SDF container.
 *
 * The current MVP uses a global volume for all queries. Local high-resolution
 * patches can be registered without changing downstream code that consumes a
 * MultiResolutionSdf.
 */
class MultiResolutionSdf {
public:
    explicit MultiResolutionSdf(std::shared_ptr<GlobalSdf> global_volume);
    explicit MultiResolutionSdf(const VoxelGrid& grid);

    const std::shared_ptr<GlobalSdf>& globalVolume() const noexcept;

    void addLocalPatch(std::shared_ptr<GlobalSdf> patch);
    const std::vector<std::shared_ptr<GlobalSdf>>& localPatches()
        const noexcept;

    bool contains(const Vec3& pos) const noexcept;
    double getDistance(const Vec3& pos) const noexcept;
    Vec3 getGradient(const Vec3& pos) const noexcept;

private:
    std::shared_ptr<GlobalSdf> global_volume;
    std::vector<std::shared_ptr<GlobalSdf>> local_patches;
};

} // namespace environment
