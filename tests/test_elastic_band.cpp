#include "Module1_Main/Types.h"
#include "environment/SdfVolume.h"
#include "optimizer/SqpProjector.h"
#include "optimizer/Trajectory.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "[FAIL] " << message << "\n";
        std::exit(1);
    }
}

double length(const environment::Vec3& value) {
    return std::sqrt(
        value.x * value.x +
        value.y * value.y +
        value.z * value.z);
}

} // namespace

int main() {
    VoxelGrid grid(20U, 20U, 5U);
    for (int z = 0; z < 5; ++z) {
        grid.setRawObstacle(10, 10, z);
    }
    const environment::GlobalSdf sdf(grid);

    const std::vector<voxel_planner::Point3D> zigzag{
        {2, 8, 2},
        {4, 9, 2},
        {6, 8, 2},
        {8, 9, 2},
        {10, 8, 2},
        {12, 9, 2},
        {14, 8, 2},
        {16, 8, 2}};

    const optimizer::ConfigurationSpaceRefiner refiner;
    const optimizer::PoseTrajectory initial =
        refiner.fromVoxelPath(zigzag, false);
    optimizer::SqpProjectorOptions options;
    options.max_iterations = 40;
    options.safe_margin = 1.0;
    options.initial_step = 0.05;
    const optimizer::SqpProjector projector(options);
    const optimizer::SqpProjectorResult result =
        projector.project(initial, sdf);

    require(result.success, "Elastic band projection did not converge.");
    require(result.trajectory.size() == initial.size(),
            "Projection changed the node count.");
    require(result.trajectory[0].position.x == initial[0].position.x &&
            result.trajectory[0].position.y == initial[0].position.y,
            "Projection moved the start endpoint.");
    require(result.trajectory[result.trajectory.size() - 1U].position.x ==
                initial[initial.size() - 1U].position.x &&
            result.trajectory[result.trajectory.size() - 1U].position.y ==
                initial[initial.size() - 1U].position.y,
            "Projection moved the goal endpoint.");

    for (const optimizer::TrajectoryNode& node : result.trajectory.nodes()) {
        require(std::isfinite(node.position.x) &&
                std::isfinite(node.position.y) &&
                std::isfinite(node.position.z),
                "Projection produced a non-finite position.");
        require(std::abs(length(node.tangent) - 1.0) < 1.0e-6,
                "Recomputed tangent is not unit length.");
        require(std::abs(length(node.normal) - 1.0) < 1.0e-6,
                "Recomputed normal is not unit length.");
    }

    VoxelGrid wallGrid(20U, 20U, 5U);
    for (int y = 0; y < 20; ++y) {
        for (int z = 0; z < 5; ++z) {
            wallGrid.setRawObstacle(10, y, z);
        }
    }
    const environment::GlobalSdf wallSdf(wallGrid);
    const optimizer::PoseTrajectory wallTrajectory =
        refiner.fromVoxelPath(
            {{2, 10, 2}, {18, 10, 2}},
            false);
    optimizer::SqpProjectorOptions wallOptions;
    wallOptions.safe_margin = 1.0;
    wallOptions.max_iterations = 4;
    const optimizer::SqpProjector wallProjector(wallOptions);
    const optimizer::SqpProjectorResult wallResult =
        wallProjector.project(wallTrajectory, wallSdf);
    require(!wallResult.success,
            "A colliding segment was accepted despite safe endpoints.");
    require(wallResult.trajectory.size() == wallTrajectory.size() &&
            wallResult.trajectory[0].position.x ==
                wallTrajectory[0].position.x &&
            wallResult.trajectory[1].position.x ==
                wallTrajectory[1].position.x,
            "Failed projection did not return the original trajectory.");

    std::cout << "[PASS] Elastic band iterations=" << result.iterations
              << " energy=" << result.final_energy
              << " gradient_norm=" << result.gradient_norm << "\n";
    return 0;
}
