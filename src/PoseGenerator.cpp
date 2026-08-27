#include "Module2_Morphology/PoseGenerator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace voxel_planner {
namespace {

struct UnitVector {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

UnitVector normalize(int x, int y, int z) {
    const double length = std::sqrt(static_cast<double>(
        x * x + y * y + z * z));
    return {
        static_cast<double>(x) / length,
        static_cast<double>(y) / length,
        static_cast<double>(z) / length};
}

UnitVector normalize(float x, float y, float z) {
    const double length = std::sqrt(
        static_cast<double>(x) * x +
        static_cast<double>(y) * y +
        static_cast<double>(z) * z);
    if (length <= 1e-9) {
        throw std::invalid_argument("Pose vector must not be zero.");
    }
    return {
        static_cast<double>(x) / length,
        static_cast<double>(y) / length,
        static_cast<double>(z) / length};
}

UnitVector cross(const UnitVector& lhs, const UnitVector& rhs) {
    return {
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x};
}

double dot(const UnitVector& lhs, const UnitVector& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

UnitVector normalized(const UnitVector& vector) {
    const double length = std::sqrt(
        vector.x * vector.x +
        vector.y * vector.y +
        vector.z * vector.z);
    return {
        vector.x / length,
        vector.y / length,
        vector.z / length};
}

UnitVector referenceAxis(const UnitVector& tangent) {
    const double absX = std::abs(tangent.x);
    const double absY = std::abs(tangent.y);
    const double absZ = std::abs(tangent.z);
    if (absX <= absY && absX <= absZ) {
        return {1.0, 0.0, 0.0};
    }
    if (absY <= absZ) {
        return {0.0, 1.0, 0.0};
    }
    return {0.0, 0.0, 1.0};
}

int crossSectionSampleCount(float physicalSize) {
    return std::max(
        1,
        static_cast<int>(std::ceil(
            static_cast<double>(physicalSize) - 1e-6)));
}

double centeredSample(int sample, int sampleCount) {
    return static_cast<double>(sample) -
        static_cast<double>(sampleCount - 1) / 2.0;
}

int roundVoxel(double value) {
    return static_cast<int>(std::floor(value + 0.5));
}

} // namespace

std::vector<BusbarPose> generateDiscretePoses(
    const PlannerConfig& config) {
    if (config.angle_step_deg <= 0 ||
        360 % config.angle_step_deg != 0) {
        throw std::invalid_argument(
            "angle_step_deg must be a positive divisor of 360.");
    }

    const int rollsPerTangent = 360 / config.angle_step_deg;
    std::vector<BusbarPose> poses;
    poses.reserve(static_cast<std::size_t>(26 * rollsPerTangent));
    std::uint32_t poseId = 0U;
    constexpr double pi = 3.14159265358979323846;

    for (int tz = -1; tz <= 1; ++tz) {
        for (int ty = -1; ty <= 1; ++ty) {
            for (int tx = -1; tx <= 1; ++tx) {
                if (tx == 0 && ty == 0 && tz == 0) {
                    continue;
                }

                const UnitVector tangent = normalize(tx, ty, tz);
                const UnitVector baseNormal = normalized(cross(
                    tangent,
                    referenceAxis(tangent)));
                const UnitVector binormal = cross(tangent, baseNormal);

                for (int roll = 0;
                     roll < 360;
                     roll += config.angle_step_deg) {
                    const double radians =
                        static_cast<double>(roll) * pi / 180.0;
                    const double cosine = std::cos(radians);
                    const double sine = std::sin(radians);
                    const UnitVector normal{
                        baseNormal.x * cosine + binormal.x * sine,
                        baseNormal.y * cosine + binormal.y * sine,
                        baseNormal.z * cosine + binormal.z * sine};

                    poses.push_back({
                        tx,
                        ty,
                        tz,
                        static_cast<float>(normal.x),
                        static_cast<float>(normal.y),
                        static_cast<float>(normal.z),
                        poseId++,
                        roll});
                }
            }
        }
    }

    return poses;
}

std::vector<VoxelOffset> generatePoseFootprint(
    const BusbarPose& pose,
    const PlannerConfig& config) {
    if (config.busbar_width <= 0.0F ||
        config.busbar_thickness <= 0.0F) {
        throw std::invalid_argument(
            "Busbar dimensions must be positive.");
    }

    const UnitVector tangent = normalize(pose.tx, pose.ty, pose.tz);
    const UnitVector normal = normalize(pose.nx, pose.ny, pose.nz);
    if (std::abs(dot(tangent, normal)) > 1e-5) {
        throw std::invalid_argument(
            "Pose normal must be orthogonal to its tangent.");
    }
    const UnitVector binormal = normalized(cross(tangent, normal));

    const int widthSamples = crossSectionSampleCount(config.busbar_width);
    const int thicknessSamples = crossSectionSampleCount(
        config.busbar_thickness);

    std::set<std::tuple<int, int, int>> uniqueOffsets;
    for (int widthSample = 0;
         widthSample < widthSamples;
         ++widthSample) {
        const double width = centeredSample(widthSample, widthSamples);
        for (int thicknessSample = 0;
             thicknessSample < thicknessSamples;
             ++thicknessSample) {
            const double thickness = centeredSample(
                thicknessSample,
                thicknessSamples);
            uniqueOffsets.emplace(
                roundVoxel(
                    binormal.x * width + normal.x * thickness),
                roundVoxel(
                    binormal.y * width + normal.y * thickness),
                roundVoxel(
                    binormal.z * width + normal.z * thickness));
        }
    }

    std::vector<VoxelOffset> footprint;
    footprint.reserve(uniqueOffsets.size());
    for (const auto& offset : uniqueOffsets) {
        footprint.push_back({
            std::get<0>(offset),
            std::get<1>(offset),
            std::get<2>(offset)});
    }
    return footprint;
}

} // namespace voxel_planner
