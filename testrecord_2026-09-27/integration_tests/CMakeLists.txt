cmake_minimum_required(VERSION 3.16)
project(BusbarIntegrationTests LANGUAGES CXX)
enable_testing()

set(BUSBAR_SOURCE_DIR "" CACHE PATH "busbar_hybrid_state directory containing src/main.cpp")
set(VOXEL_PLANNER_SOURCE_DIR "" CACHE PATH "Voxel_Planner source directory")
foreach(root BUSBAR_SOURCE_DIR VOXEL_PLANNER_SOURCE_DIR)
    if(NOT EXISTS "${${root}}/CMakeLists.txt")
        message(FATAL_ERROR "Set ${root} to the actual source directory")
    endif()
endforeach()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(MSVC)
    add_compile_options(/utf-8)
    # Voxel_Planner uses unsigned induction variables in OpenMP loops.
    set(OpenMP_CXX_FLAGS "/openmp:llvm" CACHE STRING "MSVC OpenMP mode" FORCE)
endif()
find_package(OpenMP REQUIRED COMPONENTS CXX)
add_subdirectory("${VOXEL_PLANNER_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/voxel_planner")
foreach(target test_sdf test_elastic_band)
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
    foreach(config DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
        set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${config} "${CMAKE_BINARY_DIR}/bin")
    endforeach()
endforeach()
target_compile_definitions(voxel_planner PUBLIC NOMINMAX)
set_target_properties(voxel_planner PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
foreach(config DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set_target_properties(voxel_planner PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY_${config} "${CMAKE_BINARY_DIR}/lib")
endforeach()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(voxel_planner PUBLIC stdc++fs)
endif()

# Compile the collaborator's original main and sources. Its root CMake defines
# busbar_example again, which conflicts with the dependency's demo target.
set(busbar_sources
    astar_jump_planner.cpp vector3.cpp busbar_segment.cpp segment_generator.cpp
    space_analyzer.cpp busbar_planner.cpp coarse_planner.cpp voxel_preprocessor.cpp
    busbar_planner_api.cpp debug_export.cpp debug_logger.cpp
    voxel_distance_calculator.cpp voxel_distance_api.cpp geometry_utils.cpp)
list(TRANSFORM busbar_sources PREPEND "${BUSBAR_SOURCE_DIR}/src/")
add_library(integration_busbar STATIC ${busbar_sources})
target_include_directories(integration_busbar PUBLIC "${BUSBAR_SOURCE_DIR}/include")
target_compile_definitions(integration_busbar PUBLIC BUSBAR_ROUTER_STATIC NOMINMAX)
target_link_libraries(integration_busbar PUBLIC OpenMP::OpenMP_CXX)

add_executable(hybrid_validation_cli "${BUSBAR_SOURCE_DIR}/src/main.cpp")
target_compile_definitions(hybrid_validation_cli PRIVATE BUSBAR_HAS_VOXEL_PLANNER_RELEASE=1)
target_link_libraries(hybrid_validation_cli PRIVATE integration_busbar voxel_planner)
set_target_properties(hybrid_validation_cli PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
foreach(config DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set_target_properties(hybrid_validation_cli PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY_${config} "${CMAKE_BINARY_DIR}/bin")
endforeach()
