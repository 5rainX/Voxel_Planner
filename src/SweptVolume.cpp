#include "Module2_Morphology/SweptVolume.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace voxel_planner {
namespace {

constexpr double kMaximumSurfaceStep = 0.49;
constexpr double kVectorTolerance = 1e-5;

struct Vector {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vector makeUnit(int x, int y, int z) {
    const double length = std::sqrt(static_cast<double>(
        x * x + y * y + z * z));
    if (length <= kVectorTolerance) {
        throw std::invalid_argument("Pose tangent must not be zero.");
    }
    return {
        static_cast<double>(x) / length,
        static_cast<double>(y) / length,
        static_cast<double>(z) / length};
}

Vector makeUnit(float x, float y, float z) {
    const double length = std::sqrt(
        static_cast<double>(x) * x +
        static_cast<double>(y) * y +
        static_cast<double>(z) * z);
    if (length <= kVectorTolerance) {
        throw std::invalid_argument("Pose normal must not be zero.");
    }
    return {
        static_cast<double>(x) / length,
        static_cast<double>(y) / length,
        static_cast<double>(z) / length};
}

double dot(const Vector& lhs, const Vector& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Vector cross(const Vector& lhs, const Vector& rhs) {
    return {
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x};
}

Vector normalized(const Vector& vector) {
    const double length = std::sqrt(dot(vector, vector));
    if (length <= kVectorTolerance) {
        throw std::invalid_argument("Sweep frame is degenerate.");
    }
    return {vector.x / length, vector.y / length, vector.z / length};
}

Vector rotateAroundAxis(
    const Vector& vector,
    const Vector& axis,
    double radians) {
    const double cosine = std::cos(radians);
    const double sine = std::sin(radians);
    const double parallel = dot(axis, vector) * (1.0 - cosine);
    const Vector axisCrossVector = cross(axis, vector);
    return {
        vector.x * cosine + axisCrossVector.x * sine + axis.x * parallel,
        vector.y * cosine + axisCrossVector.y * sine + axis.y * parallel,
        vector.z * cosine + axisCrossVector.z * sine + axis.z * parallel};
}

void validateDimensions(const PlannerConfig& config) {
    if (config.busbar_width <= 0.0F ||
        config.busbar_thickness <= 0.0F) {
        throw std::invalid_argument(
            "Busbar width and thickness must be positive.");
    }
}

void validateFrame(const Vector& tangent, const Vector& normal) {
    if (std::abs(dot(tangent, normal)) > kVectorTolerance) {
        throw std::invalid_argument(
            "Pose normal must be orthogonal to its tangent.");
    }
}

int roundVoxel(double value) {
    return static_cast<int>(std::floor(value + 0.5));
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

void addCrossSection(
    const Vector& center,
    const Vector& tangent,
    const Vector& normal,
    const PlannerConfig& config,
    std::set<std::tuple<int, int, int>>& offsets) {
    const Vector binormal = normalized(cross(tangent, normal));
    const int widthSamples = crossSectionSampleCount(config.busbar_width);
    const int thicknessSamples = crossSectionSampleCount(
        config.busbar_thickness);

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
            offsets.emplace(
                roundVoxel(
                    center.x + binormal.x * width + normal.x * thickness),
                roundVoxel(
                    center.y + binormal.y * width + normal.y * thickness),
                roundVoxel(
                    center.z + binormal.z * width + normal.z * thickness));
        }
    }
}

std::vector<VoxelOffset> toOffsets(
    const std::set<std::tuple<int, int, int>>& uniqueOffsets) {
    std::vector<VoxelOffset> offsets;
    offsets.reserve(uniqueOffsets.size());
    for (const auto& offset : uniqueOffsets) {
        offsets.push_back({
            std::get<0>(offset),
            std::get<1>(offset),
            std::get<2>(offset)});
    }
    return offsets;
}

} // namespace

std::vector<VoxelOffset> generateStraightSweep(
    const BusbarPose& pose,
    const PlannerConfig& config,
    float distance_s) {
    validateDimensions(config);
    if (distance_s <= 0.0F) {
        throw std::invalid_argument("Straight sweep distance must be positive.");
    }

    const Vector tangent = makeUnit(pose.tx, pose.ty, pose.tz);
    const Vector normal = makeUnit(pose.nx, pose.ny, pose.nz);
    validateFrame(tangent, normal);

    const int sampleCount = std::max(
        1,
        static_cast<int>(std::ceil(
            static_cast<double>(distance_s) / kMaximumSurfaceStep)));
    std::set<std::tuple<int, int, int>> uniqueOffsets;
    for (int sample = 0; sample <= sampleCount; ++sample) {
        const double distance =
            static_cast<double>(distance_s) * sample / sampleCount;
        addCrossSection(
            {tangent.x * distance, tangent.y * distance, tangent.z * distance},
            tangent,
            normal,
            config,
            uniqueOffsets);
    }
    return toOffsets(uniqueOffsets);
}

std::vector<VoxelOffset> generateBendSweep(
    const BusbarPose& pose_start,
    const BusbarPose& pose_end,
    const PlannerConfig& config,
    float bend_radius,
    float bend_angle_radians) {
    validateDimensions(config);
    if (bend_radius <= 0.0F || bend_angle_radians <= 0.0F ||
        bend_angle_radians >= 3.14159265358979323846F) {
        throw std::invalid_argument(
            "Bend radius must be positive and angle must be in (0, pi).");
    }

    const Vector startTangent = makeUnit(
        pose_start.tx,
        pose_start.ty,
        pose_start.tz);
    const Vector endTangent = makeUnit(
        pose_end.tx,
        pose_end.ty,
        pose_end.tz);
    const Vector startNormal = makeUnit(
        pose_start.nx,
        pose_start.ny,
        pose_start.nz);
    const Vector endNormal = makeUnit(
        pose_end.nx,
        pose_end.ny,
        pose_end.nz);
    validateFrame(startTangent, startNormal);
    validateFrame(endTangent, endNormal);

    const Vector axis = normalized(cross(startTangent, endTangent));
    const Vector expectedEndTangent = rotateAroundAxis(
        startTangent,
        axis,
        bend_angle_radians);
    if (dot(expectedEndTangent, endTangent) < 1.0 - kVectorTolerance) {
        throw std::invalid_argument(
            "Bend angle and end tangent are inconsistent.");
    }

    const Vector expectedEndNormal = rotateAroundAxis(
        startNormal,
        axis,
        bend_angle_radians);
    const double endNormalAlignment = dot(expectedEndNormal, endNormal);
    if (endNormalAlignment < 1.0 - kVectorTolerance) {
        throw std::invalid_argument(
            "Bend requires a parallel-transported end normal.");
    }

    const Vector radialStart = cross(startTangent, axis);
    const double crossSectionRadius = std::sqrt(
        std::pow(static_cast<double>(config.busbar_width) / 2.0, 2.0) +
        std::pow(static_cast<double>(config.busbar_thickness) / 2.0, 2.0));
    const double maximumSweepDistance =
        static_cast<double>(bend_radius) +
        crossSectionRadius;
    const int sampleCount = std::max(
        1,
        static_cast<int>(std::ceil(
            static_cast<double>(bend_angle_radians) *
                maximumSweepDistance /
            kMaximumSurfaceStep)));

    std::set<std::tuple<int, int, int>> uniqueOffsets;
    for (int sample = 0; sample <= sampleCount; ++sample) {
        const double angle =
            static_cast<double>(bend_angle_radians) * sample / sampleCount;
        const Vector radial = rotateAroundAxis(radialStart, axis, angle);
        const Vector tangent = rotateAroundAxis(startTangent, axis, angle);
        const Vector bendNormal = normalized(rotateAroundAxis(
            startNormal,
            axis,
            angle));
        const Vector normal = bendNormal;
        const Vector center{
            static_cast<double>(bend_radius) * (radial.x - radialStart.x),
            static_cast<double>(bend_radius) * (radial.y - radialStart.y),
            static_cast<double>(bend_radius) * (radial.z - radialStart.z)};
        addCrossSection(center, tangent, normal, config, uniqueOffsets);
    }
    return toOffsets(uniqueOffsets);
}

std::vector<VoxelOffset> generateExplicitBendSweep(
    const BusbarPose& pose_start,
    const PlannerConfig& config,
    float bend_radius,
    float bend_angle_radians,
    float axis_x,
    float axis_y,
    float axis_z) {
    validateDimensions(config);
    if (bend_radius <= 0.0F || bend_angle_radians <= 0.0F ||
        bend_angle_radians > 3.14159265358979323846F / 2.0F) {
        throw std::invalid_argument(
            "Explicit bend radius must be positive and angle must be in "
            "(0, pi/2].");
    }

    const Vector startTangent = makeUnit(
        pose_start.tx,
        pose_start.ty,
        pose_start.tz);
    const Vector startNormal = makeUnit(
        pose_start.nx,
        pose_start.ny,
        pose_start.nz);
    const Vector axis = makeUnit(axis_x, axis_y, axis_z);
    validateFrame(startTangent, startNormal);
    if (std::abs(dot(startTangent, axis)) > kVectorTolerance) {
        throw std::invalid_argument(
            "Explicit bend axis must be orthogonal to the start tangent.");
    }

    const Vector radialStart = cross(startTangent, axis);
    const double crossSectionRadius = std::sqrt(
        std::pow(static_cast<double>(config.busbar_width) / 2.0, 2.0) +
        std::pow(static_cast<double>(config.busbar_thickness) / 2.0, 2.0));
    const double maximumSweepDistance =
        static_cast<double>(bend_radius) + crossSectionRadius;
    const int sampleCount = std::max(
        1,
        static_cast<int>(std::ceil(
            static_cast<double>(bend_angle_radians) *
                maximumSweepDistance /
            kMaximumSurfaceStep)));

    std::set<std::tuple<int, int, int>> uniqueOffsets;
    for (int sample = 0; sample <= sampleCount; ++sample) {
        const double angle =
            static_cast<double>(bend_angle_radians) * sample / sampleCount;
        const Vector radial = rotateAroundAxis(radialStart, axis, angle);
        const Vector tangent = rotateAroundAxis(startTangent, axis, angle);
        const Vector normal = normalized(rotateAroundAxis(
            startNormal,
            axis,
            angle));
        const Vector center{
            static_cast<double>(bend_radius) * (radial.x - radialStart.x),
            static_cast<double>(bend_radius) * (radial.y - radialStart.y),
            static_cast<double>(bend_radius) * (radial.z - radialStart.z)};
        addCrossSection(
            center,
            tangent,
            normal,
            config,
            uniqueOffsets);
    }
    return toOffsets(uniqueOffsets);
}

std::vector<VoxelOffset> generateTwistSweep(
    const BusbarPose& pose_start,
    const BusbarPose& pose_end,
    const PlannerConfig& config,
    float twist_angle_radians) {
    validateDimensions(config);
    if (std::abs(twist_angle_radians) <= kVectorTolerance ||
        std::abs(twist_angle_radians) >
            3.14159265358979323846F) {
        throw std::invalid_argument(
            "Twist angle magnitude must be in (0, pi].");
    }

    const Vector tangent = makeUnit(
        pose_start.tx,
        pose_start.ty,
        pose_start.tz);
    const Vector endTangent = makeUnit(
        pose_end.tx,
        pose_end.ty,
        pose_end.tz);
    const Vector startNormal = makeUnit(
        pose_start.nx,
        pose_start.ny,
        pose_start.nz);
    const Vector endNormal = makeUnit(
        pose_end.nx,
        pose_end.ny,
        pose_end.nz);
    validateFrame(tangent, startNormal);
    validateFrame(endTangent, endNormal);
    if (dot(tangent, endTangent) < 1.0 - kVectorTolerance) {
        throw std::invalid_argument(
            "Twist must preserve the tangent direction.");
    }
    const Vector expectedEndNormal = rotateAroundAxis(
        startNormal,
        tangent,
        twist_angle_radians);
    if (dot(expectedEndNormal, endNormal) < 1.0 - kVectorTolerance) {
        throw std::invalid_argument(
            "Twist angle and end normal are inconsistent.");
    }

    const double crossSectionRadius = std::sqrt(
        std::pow(static_cast<double>(config.busbar_width) / 2.0, 2.0) +
        std::pow(static_cast<double>(config.busbar_thickness) / 2.0, 2.0));
    const int sampleCount = std::max(
        1,
        static_cast<int>(std::ceil(
            std::abs(static_cast<double>(twist_angle_radians)) *
            crossSectionRadius / kMaximumSurfaceStep)));

    std::set<std::tuple<int, int, int>> uniqueOffsets;
    for (int sample = 0; sample <= sampleCount; ++sample) {
        const double angle =
            static_cast<double>(twist_angle_radians) * sample / sampleCount;
        const Vector normal = normalized(rotateAroundAxis(
            startNormal,
            tangent,
            angle));
        addCrossSection(
            {0.0, 0.0, 0.0},
            tangent,
            normal,
            config,
            uniqueOffsets);
    }
    return toOffsets(uniqueOffsets);
}

} // namespace voxel_planner
