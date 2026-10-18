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

# Maintained legacy source tests live in tests/; historical run records are separate.
foreach(test_name IN ITEMS
    test_jump_astar test_pose_generation test_pose_footprint
    test_macro_diversity test_swept_volume test_mask_population
    test_se3_astar test_edge_cases)
    add_executable(${test_name} "tests/${test_name}.cpp")
    target_link_libraries(${test_name} PRIVATE voxel_planner)
    target_include_directories(${test_name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

add_executable(api_integration_test tests/api_integration_test.cpp)
target_link_libraries(api_integration_test PRIVATE voxel_planner)
add_test(NAME api_integration_test COMMAND api_integration_test
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/api_open_map.txt")

add_executable(verify_missing_flag tests/verify_missing_flag.cpp)
target_link_libraries(verify_missing_flag PRIVATE voxel_planner)
target_include_directories(verify_missing_flag PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME verify_missing_flag COMMAND verify_missing_flag
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/missing_flag_map.txt")

add_executable(generate_dataset src/tools/DatasetGenerator.cpp)
target_link_libraries(generate_dataset PRIVATE voxel_planner)
target_include_directories(generate_dataset
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Fine geometry is a separate module; the coarse planner target and API stay intact.
add_library(fine_planner STATIC
    src/fine_planner/FinePlanner.cpp
    src/fine_planner/GeometryValidator.cpp
    src/fine_planner/ReferenceExport.cpp
)
target_include_directories(fine_planner PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(fine_planner PUBLIC voxel_planner)
target_compile_features(fine_planner PUBLIC cxx_std_17)

add_executable(test_fine_planner tests/fine_planner/test_fine_planner.cpp)
target_link_libraries(test_fine_planner PRIVATE fine_planner)
target_include_directories(test_fine_planner PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME test_fine_planner COMMAND test_fine_planner)

find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME test_fine_planner_reference_integration
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/fine_planner/test_reference_integration.py"
            "$<TARGET_FILE:test_fine_planner>")
    set_tests_properties(test_fine_planner_reference_integration PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME test_twist_reference_integration
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/fine_planner/test_twist_reference_integration.py"
            "$<TARGET_FILE:test_fine_planner>")
    set_tests_properties(test_twist_reference_integration PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME test_twist_cross_validation
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/fine_planner/test_twist_cross_validation.py"
            "$<TARGET_FILE:test_fine_planner>")
    set_tests_properties(test_twist_cross_validation PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endif()
