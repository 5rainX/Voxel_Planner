# Voxel Planner

## Project Overview

Voxel Planner is a C++17 library for voxel-map pre-processing and coarse
SE(3) planning. It loads an ASCII voxel map (`.txt`, or an ASCII structured
volume in `.vtk` format), evaluates pose-aware swept-volume feasibility, and
returns cost-sorted candidate centerlines between two voxel coordinates.

The public, stable interface is declared in `include/VoxelPlannerAPI.h`.

## Prerequisites

- CMake 3.16 or newer
- A C++17-capable compiler (GCC, Clang, or MSVC)
- OpenMP development/runtime support

## Build Instructions

From the project root, run:

```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

The static library is written to `lib/` and the demo executable to `bin/`.
On a multi-configuration generator such as Visual Studio, select Release
when building (for example, `cmake --build . --config Release`).

## How to Run the Executable

Run the batch/demo caller with a voxel map path:

```sh
./bin/simulate_caller path/to/map.txt
```

On Windows, use:

```powershell
.\bin\simulate_caller.exe path\to\map.txt
```

The program prompts for busbar width and thickness (default `35 5`), start
and goal coordinates, and the maximum number of paths. It also accepts ASCII
structured `.vtk` voxel maps.

## Integration Guide (C++ API)

Link your application with the `voxel_planner` target and add the repository's
`include/` directory to your include path. The complete public workflow is:

```cpp
#include "VoxelPlannerAPI.h"

#include <iostream>

int main() {
    // Planner configuration supplied when the map context is created.
    const float busbarWidth = 35.0F;
    const float busbarThickness = 5.0F;
    const voxel_planner::ProcessedMap map =
        voxel_planner::loadMap("map.txt", busbarWidth, busbarThickness);

    const voxel_planner::Point3D start{10, 20, 8};
    const voxel_planner::Point3D goal{35, 30, 8};
    const int maxPaths = 3;
    const auto [status, paths] = voxel_planner::findPaths(
        map, start, goal, maxPaths);

    if (status == voxel_planner::PlanStatus::NO_PATH) {
        std::cout << "No path found.\n";
        return 0;
    }

    for (const voxel_planner::PathResult& result : paths) {
        std::cout << "cost=" << result.cost
                  << ", waypoints=" << result.path.size() << '\n';
        for (const voxel_planner::Point3D& point : result.path) {
            std::cout << "  voxel (" << point.x << ',' << point.y
                      << ',' << point.z << ")\n";
        }
        for (const voxel_planner::PoseDescription& pose :
             result.pose_description) {
            std::cout << "  pose normal=(" << pose.normal.x << ','
                      << pose.normal.y << ',' << pose.normal.z
                      << "), tangent=(" << pose.tangent.x << ','
                      << pose.tangent.y << ',' << pose.tangent.z << ")\n";
        }
    }
}
```

`loadMap` throws `std::invalid_argument` for invalid dimensions and
`std::runtime_error` for unreadable or malformed maps. `findPaths` returns
`PlanStatus::NO_PATH` with an empty vector when no route is available; other
invalid requests are reported as exceptions.
