#include "environment/SdfVolume.h"

#include "Module1_Main/Types.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace environment {
namespace {

constexpr float kTransformInfinity = 1.0e30F;
constexpr double kSoftWallDistance = -1.0e6;

double safeSquare(double value) noexcept {
    return value * value;
}

/**
 * @brief One-dimensional lower-envelope squared distance transform.
 *
 * The input contains squared distances accumulated along the already
 * processed axes. The coordinate spacing is included in the parabola
 * intersections, so the transform also supports anisotropic voxels.
 */
std::vector<float> transformLine(
    const std::vector<float>& input,
    double spacing) {
    const std::size_t count = input.size();
    std::vector<float> output(count, kTransformInfinity);
    if (count == 0U) {
        return output;
    }

    std::vector<std::size_t> featureSites;
    featureSites.reserve(count);
    for (std::size_t site = 0U; site < count; ++site) {
        if (std::isfinite(input[site]) &&
            input[site] < kTransformInfinity * 0.5F) {
            featureSites.push_back(site);
        }
    }
    if (featureSites.empty()) {
        return output;
    }

    std::vector<std::size_t> parabolaSites(featureSites.size());
    std::vector<double> intersections(featureSites.size() + 1U);
    std::size_t envelopeSize = 0U;
    intersections[0] = -std::numeric_limits<double>::infinity();
    intersections[1] = std::numeric_limits<double>::infinity();

    for (const std::size_t inputSite : featureSites) {
        if (envelopeSize == 0U) {
            parabolaSites[0] = inputSite;
            envelopeSize = 1U;
            intersections[0] = -std::numeric_limits<double>::infinity();
            intersections[1] = std::numeric_limits<double>::infinity();
            continue;
        }

        while (true) {
            const std::size_t previousSite =
                parabolaSites[envelopeSize - 1U];
            const double inputCoordinate =
                static_cast<double>(inputSite) * spacing;
            const double previousCoordinate =
                static_cast<double>(previousSite) * spacing;
            const double numerator =
                static_cast<double>(input[inputSite]) +
                safeSquare(inputCoordinate) -
                static_cast<double>(input[previousSite]) -
                safeSquare(previousCoordinate);
            const double denominator =
                2.0 * spacing *
                static_cast<double>(inputSite - previousSite);
            const double intersection = numerator / denominator;

            if (intersection > intersections[envelopeSize - 1U]) {
                parabolaSites[envelopeSize] = inputSite;
                intersections[envelopeSize] = intersection;
                intersections[envelopeSize + 1U] =
                    std::numeric_limits<double>::infinity();
                ++envelopeSize;
                break;
            }

            --envelopeSize;
            if (envelopeSize == 0U) {
                parabolaSites[0] = inputSite;
                envelopeSize = 1U;
                intersections[0] =
                    -std::numeric_limits<double>::infinity();
                intersections[1] =
                    std::numeric_limits<double>::infinity();
                break;
            }
        }
    }

    std::size_t activeParabola = 0U;
    for (std::size_t coordinateIndex = 0U;
         coordinateIndex < count;
         ++coordinateIndex) {
        const double coordinate =
            static_cast<double>(coordinateIndex) * spacing;
        while (activeParabola + 1U < envelopeSize &&
               intersections[activeParabola + 1U] < coordinate) {
            ++activeParabola;
        }
        const std::size_t site = parabolaSites[activeParabola];
        const double delta =
            coordinate - static_cast<double>(site) * spacing;
        output[coordinateIndex] = static_cast<float>(
            safeSquare(delta) + static_cast<double>(input[site]));
    }
    return output;
}

} // namespace

GlobalSdf::GlobalSdf(const VoxelGrid& grid)
    : GlobalSdf(grid, Vec3{}, Vec3{1.0, 1.0, 1.0}) {}

GlobalSdf::GlobalSdf(
    const VoxelGrid& grid,
    Vec3 origin,
    Vec3 voxel_size)
    : origin_(origin),
      voxel_size_(voxel_size),
      dimensions_{
          static_cast<std::size_t>(grid.width()),
          static_cast<std::size_t>(grid.height()),
          static_cast<std::size_t>(grid.depth())} {
    if (dimensions_[0] == 0U ||
        dimensions_[1] == 0U ||
        dimensions_[2] == 0U) {
        throw std::invalid_argument(
            "GlobalSdf requires a non-empty VoxelGrid.");
    }
    if (!finiteVec(origin_)) {
        throw std::invalid_argument(
            "SDF origin must contain finite coordinates.");
    }
    validateVoxelSize(voxel_size_);

    const std::size_t voxelCount = checkedVoxelCount(dimensions_);
    if (voxelCount != grid.voxelCount()) {
        throw std::invalid_argument(
            "VoxelGrid dimensions do not match its storage.");
    }

    std::size_t obstacleCount = 0U;
    for (std::size_t linearIndex = 0U;
         linearIndex < voxelCount;
         ++linearIndex) {
        if (grid.isRawObstacle(linearIndex)) {
            ++obstacleCount;
        }
    }
    if (obstacleCount == 0U) {
        distances_.assign(voxelCount, kInfinity);
        return;
    }
    if (obstacleCount == voxelCount) {
        distances_.assign(voxelCount, -kInfinity);
        return;
    }

    const std::vector<float> distanceToObstacles =
        buildSquaredDistanceField(grid, true);
    const std::vector<float> distanceToFreeSpace =
        buildSquaredDistanceField(grid, false);

    distances_.resize(voxelCount, 0.0F);
    const float offset = static_cast<float>(surfaceOffset());
    for (std::size_t z = 0U; z < dimensions_[2]; ++z) {
        for (std::size_t y = 0U; y < dimensions_[1]; ++y) {
            for (std::size_t x = 0U; x < dimensions_[0]; ++x) {
                const std::size_t linearIndex = index(x, y, z);
                const bool occupied = grid.isRawObstacle(linearIndex);
                const float squaredDistance = occupied
                    ? distanceToFreeSpace[linearIndex]
                    : distanceToObstacles[linearIndex];
                if (!std::isfinite(squaredDistance) ||
                    squaredDistance >= kInfinity * 0.5F) {
                    distances_[linearIndex] = occupied
                        ? -kInfinity
                        : kInfinity;
                    continue;
                }
                const float magnitude = std::sqrt(
                    std::max(0.0F, squaredDistance));
                const float surfaceMagnitude = std::max(
                    0.0F,
                    magnitude - offset);
                distances_[linearIndex] = occupied
                    ? -surfaceMagnitude
                    : surfaceMagnitude;
            }
        }
    }
}

const Vec3& GlobalSdf::origin() const noexcept {
    return origin_;
}

const Vec3& GlobalSdf::voxelSize() const noexcept {
    return voxel_size_;
}

const std::array<std::size_t, 3>& GlobalSdf::dimensions() const noexcept {
    return dimensions_;
}

double GlobalSdf::getDistance(const Vec3& pos) const noexcept {
    if (!contains(pos)) {
        return kSoftWallDistance;
    }

    const Vec3 coordinate = toGridCoordinate(pos);
    const double x0 = std::floor(coordinate.x);
    const double y0 = std::floor(coordinate.y);
    const double z0 = std::floor(coordinate.z);
    const std::size_t x = static_cast<std::size_t>(x0);
    const std::size_t y = static_cast<std::size_t>(y0);
    const std::size_t z = static_cast<std::size_t>(z0);
    const double fx = coordinate.x - x0;
    const double fy = coordinate.y - y0;
    const double fz = coordinate.z - z0;

    const double c000 = sample(x, y, z);
    const double c100 = sampleClamped(
        static_cast<std::ptrdiff_t>(x) + 1,
        static_cast<std::ptrdiff_t>(y),
        static_cast<std::ptrdiff_t>(z));
    const double c010 = sampleClamped(
        static_cast<std::ptrdiff_t>(x),
        static_cast<std::ptrdiff_t>(y) + 1,
        static_cast<std::ptrdiff_t>(z));
    const double c110 = sampleClamped(
        static_cast<std::ptrdiff_t>(x) + 1,
        static_cast<std::ptrdiff_t>(y) + 1,
        static_cast<std::ptrdiff_t>(z));
    const double c001 = sampleClamped(
        static_cast<std::ptrdiff_t>(x),
        static_cast<std::ptrdiff_t>(y),
        static_cast<std::ptrdiff_t>(z) + 1);
    const double c101 = sampleClamped(
        static_cast<std::ptrdiff_t>(x) + 1,
        static_cast<std::ptrdiff_t>(y),
        static_cast<std::ptrdiff_t>(z) + 1);
    const double c011 = sampleClamped(
        static_cast<std::ptrdiff_t>(x),
        static_cast<std::ptrdiff_t>(y) + 1,
        static_cast<std::ptrdiff_t>(z) + 1);
    const double c111 = sampleClamped(
        static_cast<std::ptrdiff_t>(x) + 1,
        static_cast<std::ptrdiff_t>(y) + 1,
        static_cast<std::ptrdiff_t>(z) + 1);

    const double c00 = c000 * (1.0 - fx) + c100 * fx;
    const double c10 = c010 * (1.0 - fx) + c110 * fx;
    const double c01 = c001 * (1.0 - fx) + c101 * fx;
    const double c11 = c011 * (1.0 - fx) + c111 * fx;
    const double c0 = c00 * (1.0 - fy) + c10 * fy;
    const double c1 = c01 * (1.0 - fy) + c11 * fy;
    return c0 * (1.0 - fz) + c1 * fz;
}

Vec3 GlobalSdf::getGradient(const Vec3& pos) const noexcept {
    if (!contains(pos)) {
        return directionToCenter(pos);
    }

    const Vec3 coordinate = toGridCoordinate(pos);
    const double maxX = static_cast<double>(dimensions_[0] - 1U);
    const double maxY = static_cast<double>(dimensions_[1] - 1U);
    const double maxZ = static_cast<double>(dimensions_[2] - 1U);

    const auto difference = [this](
                                const Vec3& point,
                                int axis,
                                double gridCoordinate,
                                double maxCoordinate,
                                double spacing) {
        Vec3 lower = point;
        Vec3 upper = point;
        if (axis == 0) {
            lower.x -= spacing;
            upper.x += spacing;
        } else if (axis == 1) {
            lower.y -= spacing;
            upper.y += spacing;
        } else {
            lower.z -= spacing;
            upper.z += spacing;
        }

        if (maxCoordinate <= 0.0) {
            return 0.0;
        }
        if (gridCoordinate <= 0.0) {
            const double center = getDistance(point);
            return (getDistance(upper) - center) / spacing;
        }
        if (gridCoordinate >= maxCoordinate) {
            const double center = getDistance(point);
            return (center - getDistance(lower)) / spacing;
        }
        return (getDistance(upper) - getDistance(lower)) /
            (2.0 * spacing);
    };

    return {
        difference(
            pos,
            0,
            coordinate.x,
            maxX,
            voxel_size_.x),
        difference(
            pos,
            1,
            coordinate.y,
            maxY,
            voxel_size_.y),
        difference(
            pos,
            2,
            coordinate.z,
            maxZ,
            voxel_size_.z)};
}

bool GlobalSdf::contains(const Vec3& pos) const noexcept {
    if (!finiteVec(pos)) {
        return false;
    }
    const Vec3 maximum{
        origin_.x + voxel_size_.x *
            static_cast<double>(dimensions_[0] - 1U),
        origin_.y + voxel_size_.y *
            static_cast<double>(dimensions_[1] - 1U),
        origin_.z + voxel_size_.z *
            static_cast<double>(dimensions_[2] - 1U)};
    return pos.x >= origin_.x &&
        pos.y >= origin_.y &&
        pos.z >= origin_.z &&
        pos.x <= maximum.x &&
        pos.y <= maximum.y &&
        pos.z <= maximum.z;
}

const std::vector<float>& GlobalSdf::data() const noexcept {
    return distances_;
}

std::size_t GlobalSdf::index(
    std::size_t x,
    std::size_t y,
    std::size_t z) const noexcept {
    return z * dimensions_[0] * dimensions_[1] +
        y * dimensions_[0] +
        x;
}

bool GlobalSdf::finiteVec(const Vec3& value) noexcept {
    return std::isfinite(value.x) &&
        std::isfinite(value.y) &&
        std::isfinite(value.z);
}

void GlobalSdf::validateVoxelSize(const Vec3& voxel_size) {
    if (!finiteVec(voxel_size) ||
        voxel_size.x <= 0.0 ||
        voxel_size.y <= 0.0 ||
        voxel_size.z <= 0.0) {
        throw std::invalid_argument(
            "SDF voxel size must be finite and strictly positive.");
    }
}

std::size_t GlobalSdf::checkedVoxelCount(
    const std::array<std::size_t, 3>& dimensions) {
    const auto checkedMultiply = [](
                                     std::size_t lhs,
                                     std::size_t rhs) {
        if (lhs != 0U &&
            rhs > std::numeric_limits<std::size_t>::max() / lhs) {
            throw std::overflow_error(
                "SDF dimensions overflow size_t storage.");
        }
        return lhs * rhs;
    };

    return checkedMultiply(
        checkedMultiply(dimensions[0], dimensions[1]),
        dimensions[2]);
}

std::vector<float> GlobalSdf::buildSquaredDistanceField(
    const VoxelGrid& grid,
    bool obstacleIsFeature) const {
    const std::size_t width = dimensions_[0];
    const std::size_t height = dimensions_[1];
    const std::size_t depth = dimensions_[2];
    const std::size_t voxelCount = checkedVoxelCount(dimensions_);

    std::vector<float> current(voxelCount, kInfinity);
    for (std::size_t z = 0U; z < depth; ++z) {
        for (std::size_t y = 0U; y < height; ++y) {
            for (std::size_t x = 0U; x < width; ++x) {
                const std::size_t linearIndex = index(x, y, z);
                const bool isObstacle = grid.isRawObstacle(linearIndex);
                const bool isFeature = obstacleIsFeature
                    ? isObstacle
                    : !isObstacle;
                current[linearIndex] = isFeature ? 0.0F : kInfinity;
            }
        }
    }

    std::vector<float> lineInput;
    std::vector<float> lineOutput;
    lineInput.reserve(std::max({width, height, depth}));

    std::vector<float> afterX(voxelCount, kInfinity);
    lineInput.resize(width);
    for (std::size_t z = 0U; z < depth; ++z) {
        for (std::size_t y = 0U; y < height; ++y) {
            for (std::size_t x = 0U; x < width; ++x) {
                lineInput[x] = current[index(x, y, z)];
            }
            lineOutput = transformLine(lineInput, voxel_size_.x);
            for (std::size_t x = 0U; x < width; ++x) {
                afterX[index(x, y, z)] = lineOutput[x];
            }
        }
    }

    std::vector<float> afterY(voxelCount, kInfinity);
    lineInput.resize(height);
    for (std::size_t z = 0U; z < depth; ++z) {
        for (std::size_t x = 0U; x < width; ++x) {
            for (std::size_t y = 0U; y < height; ++y) {
                lineInput[y] = afterX[index(x, y, z)];
            }
            lineOutput = transformLine(lineInput, voxel_size_.y);
            for (std::size_t y = 0U; y < height; ++y) {
                afterY[index(x, y, z)] = lineOutput[y];
            }
        }
    }

    std::vector<float> result(voxelCount, kInfinity);
    lineInput.resize(depth);
    for (std::size_t y = 0U; y < height; ++y) {
        for (std::size_t x = 0U; x < width; ++x) {
            for (std::size_t z = 0U; z < depth; ++z) {
                lineInput[z] = afterY[index(x, y, z)];
            }
            lineOutput = transformLine(lineInput, voxel_size_.z);
            for (std::size_t z = 0U; z < depth; ++z) {
                result[index(x, y, z)] = lineOutput[z];
            }
        }
    }
    return result;
}

float GlobalSdf::sample(
    std::size_t x,
    std::size_t y,
    std::size_t z) const noexcept {
    return distances_[index(x, y, z)];
}

float GlobalSdf::sampleClamped(
    std::ptrdiff_t x,
    std::ptrdiff_t y,
    std::ptrdiff_t z) const noexcept {
    const std::ptrdiff_t maxX =
        static_cast<std::ptrdiff_t>(dimensions_[0] - 1U);
    const std::ptrdiff_t maxY =
        static_cast<std::ptrdiff_t>(dimensions_[1] - 1U);
    const std::ptrdiff_t maxZ =
        static_cast<std::ptrdiff_t>(dimensions_[2] - 1U);
    const std::size_t clampedX = static_cast<std::size_t>(
        std::clamp(x, std::ptrdiff_t{0}, maxX));
    const std::size_t clampedY = static_cast<std::size_t>(
        std::clamp(y, std::ptrdiff_t{0}, maxY));
    const std::size_t clampedZ = static_cast<std::size_t>(
        std::clamp(z, std::ptrdiff_t{0}, maxZ));
    return sample(clampedX, clampedY, clampedZ);
}

Vec3 GlobalSdf::toGridCoordinate(const Vec3& pos) const noexcept {
    return {
        (pos.x - origin_.x) / voxel_size_.x,
        (pos.y - origin_.y) / voxel_size_.y,
        (pos.z - origin_.z) / voxel_size_.z};
}

Vec3 GlobalSdf::domainCenter() const noexcept {
    return {
        origin_.x + 0.5 * voxel_size_.x *
            static_cast<double>(dimensions_[0] - 1U),
        origin_.y + 0.5 * voxel_size_.y *
            static_cast<double>(dimensions_[1] - 1U),
        origin_.z + 0.5 * voxel_size_.z *
            static_cast<double>(dimensions_[2] - 1U)};
}

Vec3 GlobalSdf::directionToCenter(const Vec3& pos) const noexcept {
    if (!finiteVec(pos)) {
        return {1.0, 0.0, 0.0};
    }
    const Vec3 center = domainCenter();
    const Vec3 direction{
        center.x - pos.x,
        center.y - pos.y,
        center.z - pos.z};
    const double length = std::sqrt(
        direction.x * direction.x +
        direction.y * direction.y +
        direction.z * direction.z);
    if (!std::isfinite(length) || length <= 1.0e-12) {
        return {1.0, 0.0, 0.0};
    }
    return {
        direction.x / length,
        direction.y / length,
        direction.z / length};
}

double GlobalSdf::surfaceOffset() const noexcept {
    return 0.5 * std::max({
        voxel_size_.x,
        voxel_size_.y,
        voxel_size_.z});
}

MultiResolutionSdf::MultiResolutionSdf(
    std::shared_ptr<GlobalSdf> global_volume_in)
    : global_volume(std::move(global_volume_in)) {
    if (!global_volume) {
        throw std::invalid_argument(
            "MultiResolutionSdf requires a global SDF volume.");
    }
}

MultiResolutionSdf::MultiResolutionSdf(const VoxelGrid& grid)
    : global_volume(std::make_shared<GlobalSdf>(grid)) {}

const std::shared_ptr<GlobalSdf>&
MultiResolutionSdf::globalVolume() const noexcept {
    return global_volume;
}

void MultiResolutionSdf::addLocalPatch(std::shared_ptr<GlobalSdf> patch) {
    if (!patch) {
        throw std::invalid_argument("Local SDF patch must not be null.");
    }
    local_patches.push_back(std::move(patch));
}

const std::vector<std::shared_ptr<GlobalSdf>>&
MultiResolutionSdf::localPatches() const noexcept {
    return local_patches;
}

bool MultiResolutionSdf::contains(const Vec3& pos) const noexcept {
    for (const std::shared_ptr<GlobalSdf>& patch : local_patches) {
        if (patch && patch->contains(pos)) {
            return true;
        }
    }
    return global_volume && global_volume->contains(pos);
}

double MultiResolutionSdf::getDistance(const Vec3& pos) const noexcept {
    for (const std::shared_ptr<GlobalSdf>& patch : local_patches) {
        if (patch && patch->contains(pos)) {
            return patch->getDistance(pos);
        }
    }
    return global_volume
        ? global_volume->getDistance(pos)
        : kSoftWallDistance;
}

Vec3 MultiResolutionSdf::getGradient(const Vec3& pos) const noexcept {
    for (const std::shared_ptr<GlobalSdf>& patch : local_patches) {
        if (patch && patch->contains(pos)) {
            return patch->getGradient(pos);
        }
    }
    return global_volume
        ? global_volume->getGradient(pos)
        : Vec3{1.0, 0.0, 0.0};
}

} // namespace environment
