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
    src/environment/SdfVolume.cpp
    src/optimizer/SqpProjector.cpp
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

function(add_voxel_planner_cli target_name)
    add_executable(${target_name} Module1_Main/main.cpp)
    target_link_libraries(${target_name} PRIVATE voxel_planner)
    target_compile_options(${target_name} PRIVATE
        $<$<CONFIG:Release>:$<$<CXX_COMPILER_ID:MSVC>:/O2>>
        $<$<CONFIG:Release>:$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3>>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall;-Wextra;-Wpedantic>
    )
endfunction()

add_voxel_planner_cli(simulate_caller)
add_voxel_planner_cli(benchmark_runner)
add_voxel_planner_cli(busbar_example)

enable_testing()

add_executable(test_sdf tests/test_sdf.cpp)
target_link_libraries(test_sdf PRIVATE voxel_planner)
target_include_directories(test_sdf
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_compile_features(test_sdf PRIVATE cxx_std_17)
target_compile_options(test_sdf PRIVATE
    $<$<CONFIG:Release>:$<$<CXX_COMPILER_ID:MSVC>:/O2>>
    $<$<CONFIG:Release>:$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3>>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall;-Wextra;-Wpedantic>
)
add_test(NAME test_sdf COMMAND test_sdf)

add_executable(test_elastic_band tests/test_elastic_band.cpp)
target_link_libraries(test_elastic_band PRIVATE voxel_planner)
target_include_directories(test_elastic_band
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_compile_features(test_elastic_band PRIVATE cxx_std_17)
target_compile_options(test_elastic_band PRIVATE
    $<$<CONFIG:Release>:$<$<CXX_COMPILER_ID:MSVC>:/O2>>
    $<$<CONFIG:Release>:$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3>>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall;-Wextra;-Wpedantic>
)
add_test(NAME test_elastic_band COMMAND test_elastic_band)

add_executable(generate_dataset src/tools/DatasetGenerator.cpp)
target_link_libraries(generate_dataset PRIVATE voxel_planner)
target_include_directories(generate_dataset
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
)
