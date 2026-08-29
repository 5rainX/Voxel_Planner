cmake_minimum_required(VERSION 3.16)

project(VoxelPlanner
    VERSION 1.0.0
    DESCRIPTION "Voxel pre-processing and SE(3) coarse planner"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Keep generated artifacts in predictable, source-tree locations for handoff.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
foreach(config DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${config}
        "${CMAKE_SOURCE_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${config}
        "${CMAKE_SOURCE_DIR}/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${config}
        "${CMAKE_SOURCE_DIR}/bin")
endforeach()

find_package(OpenMP REQUIRED COMPONENTS CXX)

add_library(voxel_planner STATIC
    src/PoseGenerator.cpp
    src/PoseMaskPopulation.cpp
    src/VoxelMorphologyEngine.cpp
    src/PoseTransitionValidation.cpp
    src/SweptVolume.cpp
    src/JumpAStar.cpp
    src/VoxelPlannerAPI.cpp
)

target_compile_features(voxel_planner PUBLIC cxx_std_17)
target_include_directories(voxel_planner
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_link_libraries(voxel_planner PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(voxel_planner PRIVATE
    $<$<CONFIG:Release>:$<$<CXX_COMPILER_ID:MSVC>:/O2>>
    $<$<CONFIG:Release>:$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3>>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall;-Wextra;-Wpedantic>
)

add_executable(simulate_caller Module1_Main/main.cpp)
target_link_libraries(simulate_caller PRIVATE voxel_planner)
target_compile_options(simulate_caller PRIVATE
    $<$<CONFIG:Release>:$<$<CXX_COMPILER_ID:MSVC>:/O2>>
    $<$<CONFIG:Release>:$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3>>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall;-Wextra;-Wpedantic>
)
