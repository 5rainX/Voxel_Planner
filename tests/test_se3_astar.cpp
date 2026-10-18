#include "Module1_Main/Types.h"
#include "Module2_Morphology/PoseMaskPopulation.h"
#include "Module3_AStar/CoarseAStar.h"
#include "Module3_AStar/PoseTransitionValidation.h"
#include "Module2_Morphology/PoseGenerator.h"
#include "Module2_Morphology/SweptVolume.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "[FAIL] " << message << "\n";
        std::exit(1);
    }
}

std::uint32_t findPose(
    const std::vector<voxel_planner::BusbarPose>& poses,
    int tx,
    int ty,
    int tz) {
    for (const voxel_planner::BusbarPose& pose : poses) {
        if (pose.tx == tx && pose.ty == ty && pose.tz == tz &&
            std::abs(pose.nx) < 1e-5F &&
            std::abs(pose.ny) < 1e-5F &&
            pose.nz > 0.999F) {
            return pose.poseId;
        }
    }
    return std::numeric_limits<std::uint32_t>::max();
}

} // namespace

int main() {
    voxel_planner::PlannerConfig config;
    config.busbar_width = 3.0F;
    config.busbar_thickness = 1.0F;
    config.flat_bend_factor = 1.5F;
    config.vertical_bend_factor = 1.0F;
    config.angle_step_deg = 15;

    const std::vector<voxel_planner::BusbarPose> poses =
        voxel_planner::generateDiscretePoses(config);
    const std::uint32_t startPoseId = findPose(poses, 1, 0, 0);
    const std::uint32_t endPoseId = findPose(poses, 0, 1, 0);
    require(startPoseId < poses.size() && endPoseId < poses.size(),
            "Required orthogonal test poses were not generated.");

    VoxelGrid map(12, 12, 12);
    // The selected bend remains in Z=4. This obstacle is outside that sweep,
    // but other width orientations at (5,4,4) reach it and make the waypoint
    // pose-conditional for export verification.
    map.setRawObstacle(5, 4, 5);
    const module2_morphology::CachedPoseFootprints footprints =
        module2_morphology::precomputePoseFootprints(
            poses,
            config);
    module2_morphology::populateVoxelPoseMasks(
        map,
        poses,
        footprints);

    const PlanningResult result = module3_astar::CoarseAStar::findPaths(
        map,
        PoseState{Point3D{4, 4, 4}, startPoseId},
        PoseState{Point3D{6, 6, 4}, endPoseId},
        1,
        poses,
        config);

    constexpr float halfPi = 1.57079632679F;
    const std::vector<voxel_planner::VoxelOffset> directSweep =
        voxel_planner::generateBendSweep(
            poses[startPoseId],
            poses[endPoseId],
            config,
            1.5F,
            halfPi);
    const std::size_t directTarget = map.index(6U, 6U, 4U);
    const module3_astar::ValidationResult directValidation =
        module3_astar::validateTransitionSweep(
            map,
            Point3D{4, 4, 4},
            directSweep);
    if (!map.isPoseAllowed(directTarget, endPoseId) ||
        !directValidation.ok()) {
        std::cerr << "Direct gate static="
                  << map.isPoseAllowed(directTarget, endPoseId)
                  << " dynamic=" << directValidation.ok()
                  << " error="
                  << static_cast<int>(directValidation.error_code)
                  << "\n";
    }
    require(result.status == PlannerStatus::OK,
            "Pose-aware A* did not find the 90-degree bend.");
    require(result.paths.size() == 1U,
            "Pose-aware A* returned an unexpected path count.");
    // CoarseAStar searches a 26-neighbor voxel lattice. It does not construct
    // a manufacturing arc; the selected two diagonal steps each cost sqrt(2).
    const double expectedLatticeLength = 2.0 * std::sqrt(2.0);
    if (std::abs(result.paths.front().cost - expectedLatticeLength) >= 1e-5) {
        std::cerr << "Actual cost=" << result.paths.front().cost
                  << " expected=" << expectedLatticeLength
                  << " states=" << result.paths.front().path.size()
                  << "\n";
        for (const Point3D& point : result.paths.front().path) {
            std::cerr << "  (" << point.x << ',' << point.y << ','
                      << point.z << ")\n";
        }
    }
    require(std::abs(result.paths.front().cost - expectedLatticeLength) < 1e-5,
            "Coarse lattice cost differs from the selected step lengths.");
    // A separate straight request crosses the pose-conditional voxel.
    // The shortest diagonal route above is not required to visit it.
    const PlanningResult conditionalResult = module3_astar::CoarseAStar::findPaths(
        map,
        PoseState{Point3D{4, 4, 4}, startPoseId},
        PoseState{Point3D{6, 4, 4}, startPoseId},
        1,
        poses,
        config);
    require(conditionalResult.status == PlannerStatus::OK &&
                conditionalResult.paths.size() == 1U,
            "Conditional waypoint request failed.");
    const Point3D conditionalPoint{5, 4, 4};
    require(map.getState(
                conditionalPoint.x,
                conditionalPoint.y,
                conditionalPoint.z) == VoxelState::POSE_CONDITIONAL,
            "The fixture did not create a conditional waypoint.");
    std::size_t conditionalWaypoint = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = 1U; i < conditionalResult.paths.front().path.size(); ++i) {
        const Point3D& previous =
            conditionalResult.paths.front().path[i - 1U];
        const Point3D& current =
            conditionalResult.paths.front().path[i];
        require(std::abs(current.x - previous.x) <= 1 &&
                    std::abs(current.y - previous.y) <= 1 &&
                    std::abs(current.z - previous.z) <= 1,
                "Conditional route is not 26-neighborhood continuous.");
        if (current == conditionalPoint) {
            conditionalWaypoint = i;
        }
    }
    require(conditionalWaypoint != std::numeric_limits<std::size_t>::max(),
            "The straight route omitted its conditional waypoint.");
    bool exportedIntermediatePose = false;
    for (const Pose& pose : conditionalResult.paths.front().pose_description) {
        if (pose.waypointIndex != conditionalWaypoint) {
            continue;
        }
        const double tangentLength = std::sqrt(
            pose.tangent.x * pose.tangent.x +
            pose.tangent.y * pose.tangent.y +
            pose.tangent.z * pose.tangent.z);
        const double normalLength = std::sqrt(
            pose.normal.x * pose.normal.x +
            pose.normal.y * pose.normal.y +
            pose.normal.z * pose.normal.z);
        const double orthogonality =
            pose.normal.x * pose.tangent.x +
            pose.normal.y * pose.tangent.y +
            pose.normal.z * pose.tangent.z;
        exportedIntermediatePose =
            std::abs(tangentLength - 1.0) < 1e-6 &&
            std::abs(normalLength - 1.0) < 1e-6 &&
            std::abs(orthogonality) < 1e-6;
        break;
    }
    require(exportedIntermediatePose,
            "A conditional intermediate waypoint lacks an exported pose.");

    std::cout << "[PASS] cost=" << result.paths.front().cost
              << " expected_lattice=" << expectedLatticeLength
              << " states=" << result.paths.front().path.size()
              << " intermediate_pose_index=" << conditionalWaypoint
              << "\n";
    return 0;
}
