#include "Module2_Morphology/PoseGenerator.h"
#include "Module2_Morphology/SweptVolume.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "[FAIL] " << message << "\n";
        std::exit(1);
    }
}

} // namespace

int main() {
    voxel_planner::PlannerConfig config;
    config.busbar_width = 3.0F;
    config.busbar_thickness = 1.0F;

    const voxel_planner::BusbarPose start{
        1, 0, 0,
        0.0F, 0.0F, 1.0F,
        0U, 0};
    const voxel_planner::BusbarPose end{
        0, 1, 0,
        0.0F, 0.0F, 1.0F,
        1U, 0};

    const std::vector<voxel_planner::VoxelOffset> staticFootprint =
        voxel_planner::generatePoseFootprint(start, config);
    const std::vector<voxel_planner::VoxelOffset> straightSweep =
        voxel_planner::generateStraightSweep(start, config, 3.0F);
    require(straightSweep.size() > staticFootprint.size(),
            "Straight sweep should exceed the static segment footprint.");

    constexpr float halfPi = 1.57079632679F;
    const float flatRadius =
        config.flat_bend_factor * config.busbar_thickness;
    const std::vector<voxel_planner::VoxelOffset> bendSweep =
        voxel_planner::generateBendSweep(
            start,
            end,
            config,
            flatRadius,
            halfPi);
    require(!bendSweep.empty(), "90-degree flat bend sweep is empty.");
    require(bendSweep.size() > staticFootprint.size(),
            "Bend sweep should exceed a static footprint.");

    const voxel_planner::BusbarPose axisStarts[3] = {
        {0, 1, 0, 0.0F, 0.0F, 1.0F, 0U, 0},
        {1, 0, 0, 0.0F, 0.0F, 1.0F, 0U, 0},
        {1, 0, 0, 0.0F, 1.0F, 0.0F, 0U, 0}};
    const std::tuple<float, float, float> axes[3] = {
        {1.0F, 0.0F, 0.0F},
        {0.0F, 1.0F, 0.0F},
        {0.0F, 0.0F, 1.0F}};
    std::size_t explicitSweepCount = 0U;
    constexpr float pi = 3.14159265358979323846F;
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
        for (int degrees = 15; degrees <= 90; degrees += 15) {
            const auto [axisX, axisY, axisZ] = axes[axis];
            const std::vector<voxel_planner::VoxelOffset> explicitSweep =
                voxel_planner::generateExplicitBendSweep(
                    axisStarts[axis],
                    config,
                    4.0F,
                    degrees * pi / 180.0F,
                    axisX,
                    axisY,
                    axisZ);
            require(!explicitSweep.empty(),
                    "An explicit 15-degree-axis bend sweep is empty.");
            ++explicitSweepCount;
        }
    }
    require(explicitSweepCount == 18U,
            "Explicit bend coverage must contain six angles on three axes.");

    constexpr float fifteenDegrees = pi / 12.0F;
    const voxel_planner::BusbarPose twisted{
        1, 0, 0,
        0.0F,
        -std::sin(fifteenDegrees),
        std::cos(fifteenDegrees),
        2U,
        15};
    const std::vector<voxel_planner::VoxelOffset> twistSweep =
        voxel_planner::generateTwistSweep(
            start,
            twisted,
            config,
            fifteenDegrees);
    require(!twistSweep.empty(), "15-degree twist sweep is empty.");
    require(twistSweep.size() >= staticFootprint.size(),
            "Twist sweep does not cover its static endpoint footprints.");

    std::cout << "[PASS] static=" << staticFootprint.size()
              << " straight=" << straightSweep.size()
              << " bend=" << bendSweep.size()
              << " flat_radius=" << flatRadius
              << " twist=" << twistSweep.size()
              << " explicit_sweeps=" << explicitSweepCount << "\n";
    return 0;
}
