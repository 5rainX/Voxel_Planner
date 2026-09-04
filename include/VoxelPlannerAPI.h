#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace voxel_planner {

/**
 * @brief Public voxel classification used by the planner contract.
 */
enum class VoxelClass {
    BLOCKED = 0,
    UNCONDITIONAL = 1,
    POSE_CONDITIONAL = 2
};

/** @brief Coarse A* completion status. */
enum class PlanStatus { OK, NO_PATH };

/**
 * @brief Discrete voxel-center coordinate.
 */
struct Point3D {
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const Point3D& other) const noexcept {
        return x == other.x && y == other.y && z == other.z;
    }
};

/**
 * @brief Public three-dimensional vector.
 */
struct Vector3D {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/**
 * @brief Requested terminal orientation for a planning endpoint.
 */
struct EndpointPose {
    Vector3D normal{0.0, 0.0, 1.0};
    Vector3D tangent{1.0, 0.0, 0.0};
};

/**
 * @brief Orientation information attached to a conditional waypoint.
 */
struct PoseDescription {
    Vector3D normal{};
    Vector3D tangent{};
    std::size_t path_index = 0U;
};

/**
 * @brief Ordered waypoint constraint requested by an external caller.
 */
struct WaypointConstraint {
    Point3D point{};
    Vector3D tangent{1.0, 0.0, 0.0};
    int tolerance = 0;
};

/**
 * @brief Metadata describing where an ordered waypoint was hit.
 */
struct WaypointHit {
    std::size_t constraint_index = 0U;
    std::size_t path_index = 0U;
    Point3D point{};
    Vector3D normal{};
    Vector3D tangent{};
};

/**
 * @brief Requested terminal orientation and straight lengths.
 */
struct TerminalConstraints {
    Vector3D start_normal{0.0, 0.0, 1.0};
    Vector3D start_tangent{1.0, 0.0, 0.0};
    Vector3D end_normal{0.0, 0.0, 1.0};
    Vector3D end_tangent{1.0, 0.0, 0.0};
    float min_begin_length = 0.0F;
    float min_end_length = 0.0F;
};

/**
 * @brief Optional computation limits for coarse search stages.
 */
struct CoarseSearchLimits {
    std::size_t max_expansions_per_stage = 0U;
};

/**
 * @brief One costed centerline returned by the multi-path planner.
 */
struct PathResult {
    std::vector<Point3D> path;
    std::vector<PoseDescription> pose_description;
    std::vector<PoseDescription> centerline_poses;
    std::vector<WaypointHit> waypoint_hits;
    double cost = 0.0;
};

/**
 * @brief Loaded map and planner configuration used by routing requests.
 *
 * The implementation hides raw occupancy, pose tables, masks, and search
 * internals. Copies share the loaded map context, while route-specific
 * diversity penalties remain private to each search request.
 */
class ProcessedMap {
public:
    ProcessedMap(const ProcessedMap&) noexcept = default;
    ProcessedMap(ProcessedMap&&) noexcept = default;
    ProcessedMap& operator=(const ProcessedMap&) noexcept = default;
    ProcessedMap& operator=(ProcessedMap&&) noexcept = default;
    ~ProcessedMap();

    /**
     * @brief Exports the pose-aware voxel classification as an ASCII VTK file.
     *
     * The exported VoxelClass scalar uses 0 for BLOCKED, 1 for
     * UNCONDITIONAL, and 2 for POSE_CONDITIONAL.
     */
    void exportVoxelClassificationToVtk(
        const std::string& filename) const;

    /**
     * @brief Exports voxel classes and allowed-pose counts as ASCII VTK.
     *
     * In addition to VoxelClass, the file contains AllowedPoseCount,
     * which is zero for blocked voxels and otherwise reports the number
     * of collision-free discrete poses at each voxel.
     */
    void exportVoxelClassificationAndPosesToVtk(
        const std::string& filename) const;

private:
    struct Impl;
    std::shared_ptr<Impl> impl_;

    explicit ProcessedMap(std::shared_ptr<Impl> impl) noexcept;

    friend ProcessedMap loadMap(
        const std::string& filepath,
        float width,
        float thickness);
    friend ProcessedMap loadMap(
        const std::string& filepath,
        float width,
        float thickness,
        float flat_bend_factor,
        float vertical_bend_factor);
    friend std::pair<PlanStatus, std::vector<PathResult>> findPaths(
        const ProcessedMap& map,
        const Point3D& start,
        const Point3D& goal,
        int max_paths);
    friend std::pair<PlanStatus, std::vector<PathResult>> findPaths(
        const ProcessedMap& map,
        const Point3D& start,
        const Point3D& goal,
        int max_paths,
        const EndpointPose& start_pose,
        const EndpointPose& end_pose);
    friend std::pair<PlanStatus, std::vector<PathResult>> findPaths(
        const ProcessedMap& map,
        const Point3D& start,
        const Point3D& goal,
        int max_paths,
        const TerminalConstraints& terminal_constraints,
        const std::vector<WaypointConstraint>& waypoint_constraints,
        const CoarseSearchLimits& coarse_limits);
    friend bool isPoseAllowed(
        const ProcessedMap& map,
        const Point3D& point,
        const Vector3D& normal,
        const Vector3D& tangent);
};

/**
 * @brief Loads a voxel map and prepares its opaque SE(3) context.
 *
 * @param filepath TXT voxel-map path.
 * @param width Busbar width W in voxels.
 * @param thickness Busbar thickness H in voxels.
 * @return Reusable, opaque map context.
 * @throws std::invalid_argument for an empty path or invalid configuration.
 * @throws std::runtime_error when the map cannot be parsed.
 */
ProcessedMap loadMap(
    const std::string& filepath,
    float width,
    float thickness);

/**
 * @brief Loads a voxel map using explicit flat/vertical bend factors.
 */
ProcessedMap loadMap(
    const std::string& filepath,
    float width,
    float thickness,
    float flat_bend_factor,
    float vertical_bend_factor);

/**
 * @brief Finds cost-sorted coarse SE(3) routes in a loaded map.
 *
 * NO_PATH is returned with an empty path vector. Invalid requests,
 * computation limits, and resource failures are reported as exceptions.
 *
 * @param map Context returned by loadMap().
 * @param start Start voxel coordinate.
 * @param goal Goal voxel coordinate.
 * @param max_paths Maximum number of candidate routes requested.
 * @return Explicit status and cost-sorted candidate paths.
 */
std::pair<PlanStatus, std::vector<PathResult>> findPaths(
    const ProcessedMap& map,
    const Point3D& start,
    const Point3D& goal,
    int max_paths);

/**
 * @brief Finds routes using explicit start and goal endpoint orientations.
 */
std::pair<PlanStatus, std::vector<PathResult>> findPaths(
    const ProcessedMap& map,
    const Point3D& start,
    const Point3D& goal,
    int max_paths,
    const EndpointPose& start_pose,
    const EndpointPose& end_pose);

/**
 * @brief Finds routes with terminal, ordered-waypoint, and coarse limits.
 */
std::pair<PlanStatus, std::vector<PathResult>> findPaths(
    const ProcessedMap& map,
    const Point3D& start,
    const Point3D& goal,
    int max_paths,
    const TerminalConstraints& terminal_constraints,
    const std::vector<WaypointConstraint>& waypoint_constraints,
    const CoarseSearchLimits& coarse_limits);

/**
 * @brief Tests whether a continuous normal/tangent pose is allowed at a point.
 */
bool isPoseAllowed(
    const ProcessedMap& map,
    const Point3D& point,
    const Vector3D& normal,
    const Vector3D& tangent);

} // namespace voxel_planner
