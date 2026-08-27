#include "Module2_Morphology/PoseGenerator.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace {

using TangentKey = std::tuple<int, int, int>;

double dot(const voxel_planner::BusbarPose& pose) {
    return static_cast<double>(pose.tx) * pose.nx +
           static_cast<double>(pose.ty) * pose.ny +
           static_cast<double>(pose.tz) * pose.nz;
}

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "[FAIL] " << message << "\n";
        std::exit(1);
    }
}

} // namespace

int main() {
    voxel_planner::PlannerConfig config;
    config.angle_step_deg = 15;

    const std::vector<voxel_planner::BusbarPose> poses =
        voxel_planner::generateDiscretePoses(config);
    constexpr int expectedTangentCount = 26;
    constexpr int expectedRollCount = 360 / 15;
    constexpr double orthogonalityTolerance = 1e-6;

    require(
        poses.size() == static_cast<std::size_t>(
                            expectedTangentCount * expectedRollCount),
        "Unexpected total pose count.");

    std::set<TangentKey> expectedTangents;
    std::set<TangentKey> generatedTangents;
    std::set<std::uint32_t> poseIds;
    std::map<TangentKey, std::set<int>> rollsByTangent;

    int tangentIndex = 0;
    for (int tz = -1; tz <= 1; ++tz) {
        for (int ty = -1; ty <= 1; ++ty) {
            for (int tx = -1; tx <= 1; ++tx) {
                if (tx == 0 && ty == 0 && tz == 0) {
                    continue;
                }
                expectedTangents.emplace(tx, ty, tz);
                ++tangentIndex;
            }
        }
    }
    require(
        tangentIndex == expectedTangentCount,
        "The expected 26-neighborhood tangent table is incomplete.");

    for (const voxel_planner::BusbarPose& pose : poses) {
        const TangentKey tangent{pose.tx, pose.ty, pose.tz};
        require(
            expectedTangents.count(tangent) == 1U,
            "Generated an invalid tangent direction.");
        require(
            pose.roll_angle_deg >= 0 && pose.roll_angle_deg < 360,
            "Generated a roll angle outside [0, 360).");
        require(
            pose.roll_angle_deg % config.angle_step_deg == 0,
            "Generated a roll angle not aligned with angle_step_deg.");
        require(
            std::abs(dot(pose)) <= orthogonalityTolerance,
            "Generated normal is not orthogonal to tangent.");
        generatedTangents.insert(tangent);
        poseIds.insert(pose.poseId);
        rollsByTangent[tangent].insert(pose.roll_angle_deg);
    }

    for (const TangentKey& tangent : expectedTangents) {
        const std::set<int>& rolls = rollsByTangent.at(tangent);
        require(
            rolls.size() == expectedRollCount,
            "A tangent does not have the expected number of roll angles.");
        for (int roll = 0; roll < 360; roll += config.angle_step_deg) {
            require(
                rolls.count(roll) == 1U,
                "A tangent is missing a configured roll angle.");
        }
    }
    require(
        generatedTangents == expectedTangents,
        "Generated tangent directions do not cover the full 26-neighborhood.");
    require(
        poseIds.size() == poses.size(),
        "Pose identifiers are not unique.");

    std::cout << "[PASS] generated " << poses.size()
              << " orthogonal discrete poses\n";
    return 0;
}
