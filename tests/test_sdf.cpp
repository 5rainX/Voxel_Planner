#include "environment/SdfVolume.h"
#include "Module1_Main/Types.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>

int main() {
    VoxelGrid grid(9U, 9U, 9U);
    for (int z = 3; z <= 5; ++z) {
        for (int y = 3; y <= 5; ++y) {
            for (int x = 3; x <= 5; ++x) {
                grid.setRawObstacle(x, y, z);
            }
        }
    }

    const environment::GlobalSdf sdf(grid);
    const environment::Vec3 edge{5.5, 4.0, 4.0};
    const environment::Vec3 inside{4.0, 4.0, 4.0};
    const environment::Vec3 outside{7.0, 4.0, 4.0};
    const environment::Vec3 outOfBounds{-3.0, 4.0, 4.0};

    const double edgeDistance = sdf.getDistance(edge);
    const double insideDistance = sdf.getDistance(inside);
    const double outsideDistance = sdf.getDistance(outside);
    const double outOfBoundsDistance = sdf.getDistance(outOfBounds);
    const environment::Vec3 edgeGradient = sdf.getGradient(edge);
    const environment::Vec3 outsideGradient = sdf.getGradient(outside);
    const environment::Vec3 outOfBoundsGradient =
        sdf.getGradient(outOfBounds);

    std::cout << std::fixed << std::setprecision(6)
              << "edge distance=" << edgeDistance
              << ", gradient=("
              << edgeGradient.x << ", "
              << edgeGradient.y << ", "
              << edgeGradient.z << ")\n"
              << "inside distance=" << insideDistance << "\n"
              << "outside distance=" << outsideDistance
              << ", gradient=("
              << outsideGradient.x << ", "
              << outsideGradient.y << ", "
              << outsideGradient.z << ")\n";

    if (!(insideDistance < 0.0) ||
        !(outsideDistance > 0.0) ||
        std::abs(edgeDistance) > 1.0e-6 ||
        outsideGradient.x <= 0.0) {
        std::cerr << "[FAIL] Unexpected signed distance or gradient.\n";
        return 1;
    }

    if (sdf.contains(outOfBounds) ||
        outOfBoundsDistance > -999999.0 ||
        outOfBoundsGradient.x <= 0.0) {
        std::cerr << "[FAIL] SDF soft-wall query behavior is invalid.\n";
        return 1;
    }

    if (sdf.dimensions() != std::array<std::size_t, 3>{9U, 9U, 9U}) {
        std::cerr << "[FAIL] SDF dimensions do not match VoxelGrid.\n";
        return 1;
    }
    if (sdf.data().size() != grid.voxelCount()) {
        std::cerr << "[FAIL] SDF storage is not voxel-contiguous.\n";
        return 1;
    }

    VoxelGrid freeGrid(3U, 3U, 3U);
    const environment::GlobalSdf freeSdf(freeGrid);
    if (!std::isfinite(freeSdf.getDistance({1.0, 1.0, 1.0})) ||
        freeSdf.getDistance({1.0, 1.0, 1.0}) <= 1.0e20) {
        std::cerr << "[FAIL] All-free SDF did not return positive INF.\n";
        return 1;
    }

    VoxelGrid blockedGrid(3U, 3U, 3U);
    for (int z = 0; z < 3; ++z) {
        for (int y = 0; y < 3; ++y) {
            for (int x = 0; x < 3; ++x) {
                blockedGrid.setRawObstacle(x, y, z);
            }
        }
    }
    const environment::GlobalSdf blockedSdf(blockedGrid);
    if (!std::isfinite(blockedSdf.getDistance({1.0, 1.0, 1.0})) ||
        blockedSdf.getDistance({1.0, 1.0, 1.0}) >= -1.0e20) {
        std::cerr << "[FAIL] All-obstacle SDF did not return negative INF.\n";
        return 1;
    }

    const environment::MultiResolutionSdf multi(
        std::make_shared<environment::GlobalSdf>(sdf));
    if (!multi.contains(outside) ||
        multi.getDistance(outOfBounds) > -999999.0) {
        std::cerr << "[FAIL] MultiResolutionSdf query behavior is invalid.\n";
        return 1;
    }

    std::cout << "[PASS] GlobalSdf single-obstacle test\n";
    return 0;
}
