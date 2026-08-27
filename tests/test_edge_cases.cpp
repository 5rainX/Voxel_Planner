#include "Module1_Main/Types.h"
#include "Module2_Morphology/PoseMaskPopulation.h"
#include "Module3_AStar/CoarseAStar.h"
#include "Module3_AStar/PoseTransitionValidation.h"
#include "Module2_Morphology/PoseGenerator.h"
#include "Module2_Morphology/SweptVolume.h"

#include <cstdlib>
#include <iostream>
#include <string>
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
        if (pose.tx == tx && pose.ty == ty && pose.tz == tz) {
            return pose.poseId;
        }
    }
    return std::numeric_limits<std::uint32_t>::max();
}

} // namespace

int main() {
    voxel_planner::PlannerConfig config;
    config.busbar_width = 1.0F;
    config.busbar_thickness = 1.0F;
    config.angle_step_deg = 15;

    const std::vector<voxel_planner::BusbarPose> poses =
        voxel_planner::generateDiscretePoses(config);
    const module2_morphology::CachedPoseFootprints footprints =
        module2_morphology::precomputePoseFootprints(
            poses,
            config);

    VoxelGrid openMap(9, 9, 9);
    module2_morphology::populateVoxelPoseMasks(
        openMap,
        poses,
        footprints);
    const PoseState stationaryState{Point3D{4, 4, 4}, 0U};
    require(openMap.isPoseAllowed(
                openMap.index(4U, 4U, 4U),
                stationaryState.poseId),
            "The zero-distance fixture pose must be collision-free.");

    const PlanningResult zeroDistance =
        module3_astar::CoarseAStar::findPaths(
            openMap,
            stationaryState,
            stationaryState,
            1,
            poses,
            config);
    require(zeroDistance.status == PlannerStatus::OK &&
                zeroDistance.error_code == ErrorCode::NONE,
            "An identical start/end state must succeed immediately.");
    require(zeroDistance.paths.size() == 1U &&
                zeroDistance.paths.front().path.size() == 1U &&
                zeroDistance.paths.front().cost == 0.0,
            "A zero-distance route must contain one zero-cost waypoint.");

    const PlanningResult clampedRequest =
        module3_astar::CoarseAStar::findPaths(
            openMap,
            stationaryState,
            stationaryState,
            500,
            poses,
            config);
    require(clampedRequest.status == PlannerStatus::OK &&
                clampedRequest.paths.size() <= 10U,
            "An excessive maxPaths request must be bounded safely.");
    require(clampedRequest.message.find("clamped to 10") !=
                std::string::npos,
            "The result must report that maxPaths was clamped.");

    VoxelGrid poseMap(9, 9, 9);
    const Point3D buriedPoint{4, 4, 4};
    poseMap.setRawObstacle(
        buriedPoint.x,
        buriedPoint.y,
        buriedPoint.z);
    module2_morphology::populateVoxelPoseMasks(
        poseMap,
        poses,
        footprints);

    const PoseState buriedStart{buriedPoint, 0U};
    const PoseState validEnd{{2, 2, 2}, 0U};
    const module3_astar::ValidationResult endpointValidation =
        module3_astar::validatePoseEndpoints(
            poseMap,
            buriedStart,
            validEnd);
    require(
        endpointValidation.error_code == ErrorCode::START_POINT_BLOCKED,
        "A start pose buried in a wall must return START_POINT_BLOCKED.");

    const voxel_planner::BusbarPose outwardPose{
        -1, 0, 0,
        0.0F, 0.0F, 1.0F,
        0U, 0};
    const std::vector<voxel_planner::VoxelOffset> outwardSweep =
        voxel_planner::generateStraightSweep(
            outwardPose,
            config,
            1.0F);
    const module3_astar::ValidationResult sweepValidation =
        module3_astar::validateTransitionSweep(
            poseMap,
            Point3D{0, 0, 0},
            outwardSweep);
    require(
        sweepValidation.error_code ==
            ErrorCode::TRANSITION_OUT_OF_BOUNDS,
        "An outward edge sweep must be rejected before array indexing.");

    VoxelGrid disconnectedMap(5, 5, 5);
    for (int z = 0; z < 5; ++z) {
        for (int y = 0; y < 5; ++y) {
            disconnectedMap.setRawObstacle(2, y, z);
        }
    }
    module2_morphology::populateVoxelPoseMasks(
        disconnectedMap,
        poses,
        footprints);
    const std::uint32_t negativeXPose = findPose(poses, -1, 0, 0);
    const std::uint32_t positiveXPose = findPose(poses, 1, 0, 0);
    require(negativeXPose < poses.size() && positiveXPose < poses.size(),
            "Required endpoint poses were not generated.");
    const PlanningResult noPath = module3_astar::CoarseAStar::findPaths(
        disconnectedMap,
        PoseState{Point3D{1, 2, 2}, negativeXPose},
        PoseState{Point3D{3, 2, 2}, positiveXPose},
        1,
        poses,
        config);
    require(noPath.status == PlannerStatus::NO_PATH,
            "A disconnected map must return NO_PATH.");
    require(noPath.error_code == ErrorCode::PATH_NOT_FOUND,
            "An exhausted open set must return PATH_NOT_FOUND.");
    require(noPath.paths.empty(),
            "A failed search must return an empty path collection.");

    std::cout << "[PASS] buried_start="
              << static_cast<int>(endpointValidation.error_code)
              << " outward_sweep="
              << static_cast<int>(sweepValidation.error_code)
              << " disconnected_search="
              << static_cast<int>(noPath.error_code)
              << " zero_cost="
              << zeroDistance.paths.front().cost
              << " clamped_paths="
              << clampedRequest.paths.size()
              << "\n";
    return 0;
}
