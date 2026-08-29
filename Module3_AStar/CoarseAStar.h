#pragma once

#include "Module1_Main/Types.h"
#include "Module2_Morphology/VoxelMorphologyEngine.h"
#include "Module3_AStar/PoseTransitionValidation.h"
#include "Module2_Morphology/SweptVolume.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <new>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace module3_astar {

/**
 * @brief Pose-aware state-lattice A* over (voxelIndex, poseId).
 *
 * Every edge passes the target voxel's static pose mask and the raw-obstacle
 * swept-volume validator. Search cost is physical centerline length plus
 * route-specific soft penalties for macro-path diversity:
 * linear distance for straight actions and radius x angle for bends.
 */
class CoarseAStar {
public:
    static PlanningResult findPaths(
        const VoxelGrid& map,
        PoseState startPose,
        PoseState endPose,
        int maxPaths,
        const std::vector<voxel_planner::BusbarPose>& poses,
        const voxel_planner::PlannerConfig& config,
        Point3D rawStart,
        Point3D rawGoal,
        float endpointImmunityRadius,
        const std::vector<Point3D>* topologyHint = nullptr) {
        PlanningResult result;
        const int boundedMaxPaths = std::min(
            maxPaths,
            kMaxPathsPerRequest);
        const bool pathLimitClamped = boundedMaxPaths != maxPaths;
        if (maxPaths <= 0 || poses.empty() ||
            config.busbar_width <= 0.0F ||
            config.busbar_thickness <= 0.0F ||
            config.angle_step_deg <= 0 ||
            360 % config.angle_step_deg != 0 ||
            !std::isfinite(endpointImmunityRadius) ||
            endpointImmunityRadius < 0.0F ||
            !std::isfinite(config.max_overlap_ratio) ||
            config.max_overlap_ratio < 0.0F ||
            config.max_overlap_ratio > 1.0F ||
            (boundedMaxPaths > 1 &&
             (!std::isfinite(config.path_blocking_radius) ||
              config.path_blocking_radius <= 0.0F))) {
            result.error_code = ErrorCode::INVALID_ARGUMENT;
            result.message = "Invalid SE(3) planner arguments.";
            return result;
        }
        if (!poseTableMatchesMap(map, poses)) {
            result.error_code = ErrorCode::POSE_MASK_UNAVAILABLE;
            result.message =
                "Pose table does not match the populated voxel pose masks.";
            return result;
        }

        const ValidationResult endpoints = validatePoseEndpointsWithImmunity(
            map,
            startPose,
            endPose,
            rawStart,
            rawGoal,
            endpointImmunityRadius);
        if (!endpoints.ok()) {
            result.error_code = endpoints.error_code;
            result.message = endpoints.message;
            return result;
        }

        if (startPose.position == endPose.position &&
            startPose.poseId == endPose.poseId) {
            try {
                PathResult trivialPath;
                trivialPath.path.push_back(startPose.position);
                trivialPath.cost = 0.0;
                const std::size_t voxelIndex = map.index(
                    static_cast<std::uint32_t>(startPose.position.x),
                    static_cast<std::uint32_t>(startPose.position.y),
                    static_cast<std::uint32_t>(startPose.position.z));
                if (map.getState(voxelIndex) ==
                    VoxelState::POSE_CONDITIONAL) {
                    const voxel_planner::BusbarPose& pose =
                        poses[startPose.poseId];
                    const Vector normal = unitNormal(pose);
                    const Vector tangent = unitTangent(pose);
                    trivialPath.pose_description.push_back({
                        {normal.x, normal.y, normal.z},
                        {tangent.x, tangent.y, tangent.z},
                        0U});
                }
                result.paths.push_back(std::move(trivialPath));
            } catch (const std::bad_alloc&) {
                result.error_code = ErrorCode::MEMORY_ALLOCATION_FAILED;
                result.message =
                    "Memory allocation failed while returning the trivial "
                    "path.";
                return result;
            } catch (const std::exception& error) {
                result.error_code = ErrorCode::INVALID_ARGUMENT;
                result.message = error.what();
                return result;
            }
            result.status = PlannerStatus::OK;
            result.error_code = ErrorCode::NONE;
            result.message = "Start and end states are identical; returned "
                             "the zero-distance path.";
            if (pathLimitClamped) {
                result.message += " maxPaths was clamped to " +
                    std::to_string(kMaxPathsPerRequest) + ".";
            }
            return result;
        }

        ActionCatalog actions;
        try {
            actions = buildActionCatalog(poses, config);
        } catch (const std::bad_alloc&) {
            result.error_code = ErrorCode::MEMORY_ALLOCATION_FAILED;
            result.message =
                "Memory allocation failed while building SE(3) actions.";
            return result;
        } catch (const std::exception& error) {
            result.error_code = ErrorCode::INVALID_ARGUMENT;
            result.message = error.what();
            return result;
        }

        const std::uint64_t startKey = stateKey(
            map.index(
                static_cast<std::uint32_t>(startPose.position.x),
                static_cast<std::uint32_t>(startPose.position.y),
                static_cast<std::uint32_t>(startPose.position.z)),
            startPose.poseId,
            poses.size());
        const std::uint64_t goalKey = stateKey(
            map.index(
                static_cast<std::uint32_t>(endPose.position.x),
                static_cast<std::uint32_t>(endPose.position.y),
                static_cast<std::uint32_t>(endPose.position.z)),
            endPose.poseId,
            poses.size());

        BackwardDistanceField backwardDistance;
        try {
            backwardDistance = buildBackwardDistanceField(
                map,
                startKey / poses.size(),
                goalKey / poses.size(),
                rawStart,
                rawGoal,
                endpointImmunityRadius,
                topologyHint,
                std::max(config.busbar_width * 2.0F, 30.0F));
        } catch (const std::bad_alloc&) {
            result.error_code = ErrorCode::MEMORY_ALLOCATION_FAILED;
            result.message =
                "Memory allocation failed while building the backward "
                "distance field.";
            return result;
        }
        CostPenaltyMap penalties;

        int penaltyRetries = 0;
        int rejectedOverlapCandidates = 0;
        bool diversityRetryLimitReached = false;
        while (static_cast<int>(result.paths.size()) < boundedMaxPaths) {
            SearchOutcome outcome;
            try {
                outcome = runSingleSearch(
                    map,
                    startKey,
                    goalKey,
                    poses,
                    config,
                    actions,
                    backwardDistance,
                    penalties,
                    rawStart,
                    rawGoal,
                    endpointImmunityRadius);
            } catch (const std::bad_alloc&) {
                std::sort(
                    result.paths.begin(),
                    result.paths.end(),
                    [](const PathResult& lhs, const PathResult& rhs) {
                        return lhs.cost < rhs.cost;
                    });
                result.status = result.paths.empty()
                    ? PlannerStatus::NO_PATH
                    : PlannerStatus::OK;
                result.error_code = ErrorCode::MEMORY_ALLOCATION_FAILED;
                result.message =
                    "Memory allocation failed during SE(3) search; "
                    "returning paths completed so far.";
                return result;
            }
            if (outcome.path.path.empty()) {
                if (outcome.error_code ==
                    ErrorCode::COMPUTATION_LIMIT_EXCEEDED) {
                    std::sort(
                        result.paths.begin(),
                        result.paths.end(),
                        [](const PathResult& lhs, const PathResult& rhs) {
                            return lhs.cost < rhs.cost;
                        });
                    result.status = result.paths.empty()
                        ? PlannerStatus::NO_PATH
                        : PlannerStatus::OK;
                    result.error_code =
                        ErrorCode::COMPUTATION_LIMIT_EXCEEDED;
                    result.message =
                        "SE(3) expansion limit exceeded after " +
                        std::to_string(outcome.expansions) +
                        " expansions; closest goal distance=" +
                        std::to_string(outcome.closestGoalDistance) +
                        ", generated states=" +
                        std::to_string(outcome.generatedStates) +
                        ", open states=" +
                        std::to_string(outcome.remainingOpenStates) +
                        ", backward BFS reachable voxels=" +
                        std::to_string(
                            outcome.backwardReachableVoxels) +
                        ", backward start distance=" +
                        std::to_string(
                            outcome.backwardStartDistance) +
                        ", closest backward distance=" +
                        std::to_string(
                            outcome.closestBackwardDistance) +
                        ", best generated goal alignment=" +
                        std::to_string(
                            outcome.bestGeneratedGoalAlignment) +
                        "; returning paths completed so far.";
                    return result;
                }
                break;
            }

            try {
                if (!result.paths.empty() &&
                    pathOverlapRatio(
                        map,
                        outcome.path.path,
                        result.paths) > config.max_overlap_ratio) {
                    ++penaltyRetries;
                    ++rejectedOverlapCandidates;
                    addSoftPathPenalty(
                        map,
                        outcome.path.path,
                        rawStart,
                        rawGoal,
                        endpointImmunityRadius,
                        penalties);
                    if (penaltyRetries < kMaxPenaltyRetriesPerPath) {
                        continue;
                    }
                    diversityRetryLimitReached = true;
                    break;
                }
                penaltyRetries = 0;
                result.paths.push_back(std::move(outcome.path));
                if (static_cast<int>(result.paths.size()) <
                    boundedMaxPaths) {
                    addSoftPathPenalty(
                        map,
                        result.paths.back().path,
                        rawStart,
                        rawGoal,
                        endpointImmunityRadius,
                        penalties);
                }
            } catch (const std::bad_alloc&) {
                std::sort(
                    result.paths.begin(),
                    result.paths.end(),
                    [](const PathResult& lhs, const PathResult& rhs) {
                        return lhs.cost < rhs.cost;
                    });
                result.status = result.paths.empty()
                    ? PlannerStatus::NO_PATH
                    : PlannerStatus::OK;
                result.error_code = ErrorCode::MEMORY_ALLOCATION_FAILED;
                result.message =
                    "Memory allocation failed while storing or blocking a "
                    "path; returning paths completed so far.";
                return result;
            }
        }

        if (result.paths.empty()) {
            result.error_code = ErrorCode::PATH_NOT_FOUND;
            result.message = "No feasible pose-aware path exists.";
            return result;
        }

        std::sort(
            result.paths.begin(),
            result.paths.end(),
            [](const PathResult& lhs, const PathResult& rhs) {
                return lhs.cost < rhs.cost;
            });
        result.status = PlannerStatus::OK;
        result.error_code = ErrorCode::NONE;
        result.message =
            "Generated " + std::to_string(result.paths.size()) +
            " pose-aware path(s).";
        if (rejectedOverlapCandidates != 0) {
            result.message += " Rejected " +
                std::to_string(rejectedOverlapCandidates) +
                " candidate(s) above the overlap threshold.";
        }
        if (diversityRetryLimitReached) {
            result.message +=
                " Stopped after exhausting diversity retries without "
                "accepting an over-threshold path.";
        }
        if (pathLimitClamped) {
            result.message += " maxPaths was clamped to " +
                std::to_string(kMaxPathsPerRequest) + ".";
        }
        return result;
    }

    static PlanningResult findPaths(
        const VoxelGrid& map,
        PoseState startPose,
        PoseState endPose,
        int maxPaths,
        const std::vector<voxel_planner::BusbarPose>& poses,
        const voxel_planner::PlannerConfig& config) {
        return findPaths(
            map,
            startPose,
            endPose,
            maxPaths,
            poses,
            config,
            startPose.position,
            endPose.position,
            0.0F);
    }

    /**
     * @brief Compatibility entry point that selects the first valid endpoint
     * poses and delegates to the same SE(3) search loop.
     */
    static PlanningResult findPaths(
        const VoxelGrid& map,
        Point3D start,
        Point3D goal,
        int maxPaths,
        const std::vector<voxel_planner::BusbarPose>& poses,
        const voxel_planner::PlannerConfig& config) {
        PlanningResult result;
        if (!map.isValid(start.x, start.y, start.z)) {
            result.error_code = ErrorCode::START_OUT_OF_BOUNDS;
            result.message = "Start point is outside the voxel map.";
            return result;
        }
        if (!map.isValid(goal.x, goal.y, goal.z)) {
            result.error_code = ErrorCode::END_OUT_OF_BOUNDS;
            result.message = "Goal point is outside the voxel map.";
            return result;
        }
        if (!poseTableMatchesMap(map, poses)) {
            result.error_code = ErrorCode::POSE_MASK_UNAVAILABLE;
            result.message = "Voxel pose masks have not been populated.";
            return result;
        }

        const std::size_t startIndex = map.index(
            static_cast<std::uint32_t>(start.x),
            static_cast<std::uint32_t>(start.y),
            static_cast<std::uint32_t>(start.z));
        const std::size_t goalIndex = map.index(
            static_cast<std::uint32_t>(goal.x),
            static_cast<std::uint32_t>(goal.y),
            static_cast<std::uint32_t>(goal.z));
        const std::uint32_t sharedPoseId = bestSharedEndpointPose(
            map,
            startIndex,
            goalIndex,
            start,
            goal,
            poses);
        const std::uint32_t startPoseId = sharedPoseId != kInvalidPoseId
            ? sharedPoseId
            : firstAllowedPose(map, startIndex);
        const std::uint32_t goalPoseId = sharedPoseId != kInvalidPoseId
            ? sharedPoseId
            : firstAllowedPose(map, goalIndex);
        if (startPoseId == kInvalidPoseId) {
            result.error_code = ErrorCode::START_POINT_BLOCKED;
            result.message = "No busbar pose is valid at the start point.";
            return result;
        }
        if (goalPoseId == kInvalidPoseId) {
            result.error_code = ErrorCode::END_POINT_BLOCKED;
            result.message = "No busbar pose is valid at the goal point.";
            return result;
        }
        return findPaths(
            map,
            PoseState{start, startPoseId},
            PoseState{goal, goalPoseId},
            maxPaths,
            poses,
            config);
    }

    /**
     * @brief Obsolete point-mass signature retained only for a clear migration
     * error; it never runs the removed point-mass loop.
     */
    static PlanningResult findPaths(
        const VoxelGrid&,
        Point3D,
        Point3D,
        int) {
        PlanningResult result;
        result.error_code = ErrorCode::POSE_MASK_UNAVAILABLE;
        result.message =
            "SE(3) search requires populated pose masks, a pose table, and "
            "PlannerConfig.";
        return result;
    }

private:
    struct Vector {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
    };

    struct BendAction {
        std::uint32_t toPoseId = 0U;
        int dx = 0;
        int dy = 0;
        int dz = 0;
        float radius = 0.0F;
        float angleRadians = 0.0F;
        float physicalLength = 0.0F;
        float axisX = 0.0F;
        float axisY = 0.0F;
        float axisZ = 0.0F;
        int requestedDegrees = 0;
    };

    struct TwistAction {
        std::uint32_t toPoseId = 0U;
        float angleRadians = 0.0F;
    };

    struct ActionCatalog {
        std::vector<std::vector<voxel_planner::VoxelOffset>> straightSweeps;
        std::vector<voxel_planner::VoxelOffsetBounds> straightSweepBounds;
        std::vector<float> straightLengths;
        std::vector<std::vector<BendAction>> bends;
        std::vector<std::vector<std::vector<voxel_planner::VoxelOffset>>>
            bendSweeps;
        std::vector<std::vector<voxel_planner::VoxelOffsetBounds>>
            bendSweepBounds;
        std::vector<std::vector<TwistAction>> twists;
        std::vector<std::vector<std::vector<voxel_planner::VoxelOffset>>>
            twistSweeps;
        std::vector<std::vector<voxel_planner::VoxelOffsetBounds>>
            twistSweepBounds;
        std::vector<std::uint32_t> poseEquivalenceClass;
        std::size_t poseEquivalenceClassCount = 0U;
        float heuristicScale = 1.0F;
    };

    struct NeighborCandidate {
        std::size_t targetVoxel = 0U;
        std::uint32_t targetPoseId = 0U;
        std::uint32_t transitionAction =
            std::numeric_limits<std::uint32_t>::max() - 1U;
        float physicalLength = 0.0F;
        Point3D targetPoint{};
        bool evaluate = false;
        bool endpointImmune = false;
        bool valid = false;
    };

    struct Node {
        std::uint64_t stateKey = 0U;
        std::uint32_t voxelIndex = 0U;
        std::uint32_t poseId = 0U;
        float g = std::numeric_limits<float>::infinity();
        float f = std::numeric_limits<float>::infinity();
        float h = std::numeric_limits<float>::infinity();
    };

    struct NodeCompare {
        bool operator()(const Node& lhs, const Node& rhs) const {
            if (lhs.f != rhs.f) {
                return lhs.f > rhs.f;
            }
            if (lhs.h != rhs.h) {
                return lhs.h > rhs.h;
            }
            if (lhs.g != rhs.g) {
                return lhs.g < rhs.g;
            }
            return lhs.stateKey > rhs.stateKey;
        }
    };

    static constexpr std::uint32_t kInvalidPoseId =
        std::numeric_limits<std::uint32_t>::max();
    static constexpr std::uint64_t kNoParent =
        std::numeric_limits<std::uint64_t>::max();
    static constexpr std::uint32_t kStraightAction =
        std::numeric_limits<std::uint32_t>::max();
    static constexpr std::uint32_t kNoAction =
        std::numeric_limits<std::uint32_t>::max() - 1U;
    static constexpr std::uint32_t kTwistActionBase = 0x80000000U;

    struct StateRecord {
        float g = std::numeric_limits<float>::infinity();
        std::uint64_t parent = kNoParent;
        std::uint32_t parentAction = kNoAction;
        std::uint8_t status = 0U;
    };

    struct DominanceRecord {
        float g = std::numeric_limits<float>::infinity();
        std::uint64_t stateKey = kNoParent;
    };

    struct SearchOutcome {
        PathResult path;
        std::vector<std::uint64_t> statePath;
        ErrorCode error_code = ErrorCode::PATH_NOT_FOUND;
        std::size_t expansions = 0U;
        std::size_t generatedStates = 0U;
        std::size_t remainingOpenStates = 0U;
        float closestGoalDistance =
            std::numeric_limits<float>::infinity();
        std::size_t backwardReachableVoxels = 0U;
        std::uint32_t backwardStartDistance =
            std::numeric_limits<std::uint32_t>::max();
        std::uint32_t closestBackwardDistance =
            std::numeric_limits<std::uint32_t>::max();
        float bestGeneratedGoalAlignment = -1.0F;
    };

    struct BackwardDistanceField {
        static constexpr std::uint16_t kCompactUnreachable =
            std::numeric_limits<std::uint16_t>::max();
        static constexpr std::uint32_t kWideUnreachable =
            std::numeric_limits<std::uint32_t>::max();

        std::uint32_t width = 0U;
        std::uint32_t height = 0U;
        std::uint32_t depth = 0U;
        std::vector<std::uint16_t> compactDistances;
        std::vector<std::uint32_t> wideDistances;
        std::size_t reachableVoxels = 0U;

        std::uint32_t distance(std::size_t voxelIndex) const noexcept {
            if (!compactDistances.empty()) {
                const std::uint16_t value = compactDistances[voxelIndex];
                return value == kCompactUnreachable
                    ? kWideUnreachable
                    : static_cast<std::uint32_t>(value);
            }
            return wideDistances.empty()
                ? kWideUnreachable
                : wideDistances[voxelIndex];
        }
    };

    struct CostPenaltyMap {
        std::unordered_map<std::size_t, float> values;

        float get(std::size_t voxelIndex) const noexcept {
            const auto entry = values.find(voxelIndex);
            return entry == values.end() ? 0.0F : entry->second;
        }
    };

    static constexpr std::uint8_t kOpen = 1U;
    static constexpr std::uint8_t kClosed = 2U;
    static constexpr int kMaxPathsPerRequest = 10;
    static constexpr int kNeighborThreadCount = 6;
    static constexpr std::size_t kNeighborParallelThreshold = 12U;
    static constexpr float kHeuristicWeight = 10.0F;
    static constexpr float kGoalOrientationPenaltyDistance = 5.0F;
    static constexpr float kTravelHeadingPenaltyDistance = 20.0F;
    static constexpr int kMaxPenaltyRetriesPerPath = 10;
    static constexpr std::size_t kMaxExpansionsPerPath = 50000U;
    static constexpr std::size_t kMaxExpansionsPerPenalizedPath = 100000U;
    static constexpr double kPi = 3.14159265358979323846;
    static constexpr double kGeometryTolerance = 1e-4;
    static constexpr float kDominanceCostTolerance = 1e-5F;
    static constexpr float kCenterlineReusePenalty = 100000.0F;
    static constexpr float kNeighborReusePenalty = 1000.0F;

    static const std::vector<Point3D>& emptyTopologyHint() {
        static const std::vector<Point3D> empty;
        return empty;
    }

    static bool isTwistAction(std::uint32_t action) noexcept {
        return action != kStraightAction &&
            action != kNoAction &&
            (action & kTwistActionBase) != 0U;
    }

    static std::size_t twistActionIndex(std::uint32_t action) noexcept {
        return static_cast<std::size_t>(action & ~kTwistActionBase);
    }

    static bool isNearEndpoint(
        const Point3D& point,
        const Point3D& rawStart,
        const Point3D& rawGoal,
        float radius) {
        if (radius <= 0.0F) {
            return false;
        }
        const auto withinRadius = [&point, radius](const Point3D& endpoint) {
            const double dx = static_cast<double>(point.x - endpoint.x);
            const double dy = static_cast<double>(point.y - endpoint.y);
            const double dz = static_cast<double>(point.z - endpoint.z);
            return dx * dx + dy * dy + dz * dz <
                static_cast<double>(radius) * radius;
        };
        return withinRadius(rawStart) || withinRadius(rawGoal);
    }

    static ValidationResult validatePoseEndpointsWithImmunity(
        const VoxelGrid& map,
        const PoseState& startPose,
        const PoseState& endPose,
        const Point3D& rawStart,
        const Point3D& rawGoal,
        float immunityRadius) {
        if (!map.isValid(
                startPose.position.x,
                startPose.position.y,
                startPose.position.z)) {
            return {ErrorCode::START_OUT_OF_BOUNDS,
                    "Start pose is outside the voxel map."};
        }
        if (!map.isValid(
                endPose.position.x,
                endPose.position.y,
                endPose.position.z)) {
            return {ErrorCode::END_OUT_OF_BOUNDS,
                    "End pose is outside the voxel map."};
        }
        if (map.poseCount() == 0U || map.poseMaskWordCount() == 0U) {
            return {ErrorCode::POSE_MASK_UNAVAILABLE,
                    "Voxel pose masks have not been populated."};
        }
        if (startPose.poseId >= map.poseCount()) {
            return {ErrorCode::START_POSE_INVALID,
                    "Start poseId is outside the generated pose table."};
        }
        if (endPose.poseId >= map.poseCount()) {
            return {ErrorCode::END_POSE_INVALID,
                    "End poseId is outside the generated pose table."};
        }

        const std::size_t startIndex = map.index(
            static_cast<std::uint32_t>(startPose.position.x),
            static_cast<std::uint32_t>(startPose.position.y),
            static_cast<std::uint32_t>(startPose.position.z));
        if (!isNearEndpoint(
                startPose.position,
                rawStart,
                rawGoal,
                immunityRadius) &&
            !map.isPoseAllowedForSearch(startIndex, startPose.poseId)) {
            return {ErrorCode::START_POINT_BLOCKED,
                    "Start pose collides with an obstacle or map boundary."};
        }

        const std::size_t endIndex = map.index(
            static_cast<std::uint32_t>(endPose.position.x),
            static_cast<std::uint32_t>(endPose.position.y),
            static_cast<std::uint32_t>(endPose.position.z));
        if (!isNearEndpoint(
                endPose.position,
                rawStart,
                rawGoal,
                immunityRadius) &&
            !map.isPoseAllowedForSearch(endIndex, endPose.poseId)) {
            return {ErrorCode::END_POINT_BLOCKED,
                    "End pose collides with an obstacle or map boundary."};
        }
        return {};
    }

    static bool isSweepCollisionFreeForSearch(
        const VoxelGrid& map,
        const Point3D& anchor,
        const std::vector<voxel_planner::VoxelOffset>& sweep,
        const voxel_planner::VoxelOffsetBounds& bounds) noexcept {
        if (sweep.empty()) {
            return false;
        }
        if (map.morphologyPrefix_) {
            const bool prefixFree =
                module2_morphology::VoxelMorphologyEngine::
                    isBoxCollisionFree(
                        *map.morphologyPrefix_,
                        map.prefixWidth_,
                        map.prefixHeight_,
                        map.prefixDepth_,
                        anchor.x,
                        anchor.y,
                        anchor.z,
                        bounds);
            if (prefixFree) {
                return true;
            }
        } else if (map.isObstacleFreeBox(
                       anchor.x,
                       anchor.y,
                       anchor.z,
                       bounds)) {
            return true;
        }

        // A rotated footprint can make its enclosing box much larger than
        // the actual sweep. Split it into three slabs and use the same
        // O(1) query on each local box before falling back to exact offsets.
        if (map.morphologyPrefix_ && bounds.valid) {
            const int spanX = bounds.maxDx - bounds.minDx;
            const int spanY = bounds.maxDy - bounds.minDy;
            const int spanZ = bounds.maxDz - bounds.minDz;
            const int axis = spanX >= spanY && spanX >= spanZ
                ? 0
                : (spanY >= spanZ ? 1 : 2);
            voxel_planner::VoxelOffsetBounds slabs[3];
            std::uint8_t used[3] = {0U, 0U, 0U};
            const int span = axis == 0
                ? spanX
                : (axis == 1 ? spanY : spanZ);
            for (const voxel_planner::VoxelOffset& offset : sweep) {
                const int coordinate = axis == 0
                    ? offset.dx
                    : (axis == 1 ? offset.dy : offset.dz);
                const int relative = coordinate -
                    (axis == 0
                        ? bounds.minDx
                        : (axis == 1 ? bounds.minDy : bounds.minDz));
                const int slab = span == 0
                    ? 0
                    : std::min(2, (relative * 3) / (span + 1));
                if (used[slab] == 0U) {
                    slabs[slab].minDx = slabs[slab].maxDx = offset.dx;
                    slabs[slab].minDy = slabs[slab].maxDy = offset.dy;
                    slabs[slab].minDz = slabs[slab].maxDz = offset.dz;
                    slabs[slab].valid = true;
                    used[slab] = 1U;
                } else {
                    slabs[slab].minDx = std::min(
                        slabs[slab].minDx,
                        offset.dx);
                    slabs[slab].maxDx = std::max(
                        slabs[slab].maxDx,
                        offset.dx);
                    slabs[slab].minDy = std::min(
                        slabs[slab].minDy,
                        offset.dy);
                    slabs[slab].maxDy = std::max(
                        slabs[slab].maxDy,
                        offset.dy);
                    slabs[slab].minDz = std::min(
                        slabs[slab].minDz,
                        offset.dz);
                    slabs[slab].maxDz = std::max(
                        slabs[slab].maxDz,
                        offset.dz);
                }
            }
            bool allSlabsFree = true;
            for (int slab = 0; slab < 3; ++slab) {
                if (used[slab] != 0U &&
                    !module2_morphology::VoxelMorphologyEngine::
                        isBoxCollisionFree(
                        *map.morphologyPrefix_,
                        map.prefixWidth_,
                        map.prefixHeight_,
                        map.prefixDepth_,
                        anchor.x,
                        anchor.y,
                        anchor.z,
                        slabs[slab])) {
                    allSlabsFree = false;
                    break;
                }
            }
            if (allSlabsFree) {
                return true;
            }
        } else if (map.isObstacleFreeBox(
                       anchor.x,
                       anchor.y,
                       anchor.z,
                       bounds)) {
            return true;
        }
        for (const voxel_planner::VoxelOffset& offset : sweep) {
            const int x = anchor.x + offset.dx;
            const int y = anchor.y + offset.dy;
            const int z = anchor.z + offset.dz;
            if (!map.isValid(x, y, z) ||
                map.isObstacle(map.index(
                    static_cast<std::uint32_t>(x),
                    static_cast<std::uint32_t>(y),
                    static_cast<std::uint32_t>(z)))) {
                return false;
            }
        }
        return true;
    }

    static voxel_planner::VoxelOffsetBounds offsetBounds(
        const std::vector<voxel_planner::VoxelOffset>& offsets) {
        voxel_planner::VoxelOffsetBounds bounds;
        if (offsets.empty()) {
            return bounds;
        }
        bounds.minDx = bounds.maxDx = offsets.front().dx;
        bounds.minDy = bounds.maxDy = offsets.front().dy;
        bounds.minDz = bounds.maxDz = offsets.front().dz;
        for (const voxel_planner::VoxelOffset& offset : offsets) {
            bounds.minDx = std::min(bounds.minDx, offset.dx);
            bounds.maxDx = std::max(bounds.maxDx, offset.dx);
            bounds.minDy = std::min(bounds.minDy, offset.dy);
            bounds.maxDy = std::max(bounds.maxDy, offset.dy);
            bounds.minDz = std::min(bounds.minDz, offset.dz);
            bounds.maxDz = std::max(bounds.maxDz, offset.dz);
        }
        bounds.valid = true;
        return bounds;
    }

    static bool poseTableMatchesMap(
        const VoxelGrid& map,
        const std::vector<voxel_planner::BusbarPose>& poses) {
        if (poses.empty() || map.poseCount() != poses.size()) {
            return false;
        }
        std::vector<std::uint8_t> seen(poses.size(), 0U);
        for (const voxel_planner::BusbarPose& pose : poses) {
            if (pose.poseId >= poses.size() || seen[pose.poseId] != 0U) {
                return false;
            }
            seen[pose.poseId] = 1U;
        }
        return true;
    }

    static std::uint32_t firstAllowedPose(
        const VoxelGrid& map,
        std::size_t voxelIndex) {
        for (std::uint32_t poseId = 0U;
             poseId < map.poseCount();
             ++poseId) {
            if (map.isPoseAllowed(voxelIndex, poseId)) {
                return poseId;
            }
        }
        return kInvalidPoseId;
    }

    static std::uint32_t bestSharedEndpointPose(
        const VoxelGrid& map,
        std::size_t startIndex,
        std::size_t goalIndex,
        const Point3D& start,
        const Point3D& goal,
        const std::vector<voxel_planner::BusbarPose>& poses) {
        const double dx = static_cast<double>(goal.x - start.x);
        const double dy = static_cast<double>(goal.y - start.y);
        const double dz = static_cast<double>(goal.z - start.z);
        double bestAlignment = -std::numeric_limits<double>::infinity();
        std::uint32_t selected = kInvalidPoseId;
        for (const voxel_planner::BusbarPose& pose : poses) {
            if (!map.isPoseAllowed(startIndex, pose.poseId) ||
                !map.isPoseAllowed(goalIndex, pose.poseId)) {
                continue;
            }
            const Vector tangent = unitTangent(pose);
            const double alignment =
                tangent.x * dx + tangent.y * dy + tangent.z * dz;
            if (alignment > bestAlignment) {
                bestAlignment = alignment;
                selected = pose.poseId;
            }
        }
        return selected;
    }

    static Vector unitTangent(const voxel_planner::BusbarPose& pose) {
        const double length = std::sqrt(static_cast<double>(
            pose.tx * pose.tx + pose.ty * pose.ty + pose.tz * pose.tz));
        if (length <= kGeometryTolerance) {
            throw std::invalid_argument("Pose tangent must not be zero.");
        }
        return {pose.tx / length, pose.ty / length, pose.tz / length};
    }

    static Vector unitNormal(const voxel_planner::BusbarPose& pose) {
        const double length = std::sqrt(
            static_cast<double>(pose.nx) * pose.nx +
            static_cast<double>(pose.ny) * pose.ny +
            static_cast<double>(pose.nz) * pose.nz);
        if (length <= kGeometryTolerance) {
            throw std::invalid_argument("Pose normal must not be zero.");
        }
        return {pose.nx / length, pose.ny / length, pose.nz / length};
    }

    static double dot(const Vector& lhs, const Vector& rhs) {
        return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
    }

    static Vector cross(const Vector& lhs, const Vector& rhs) {
        return {
            lhs.y * rhs.z - lhs.z * rhs.y,
            lhs.z * rhs.x - lhs.x * rhs.z,
            lhs.x * rhs.y - lhs.y * rhs.x};
    }

    static Vector normalized(const Vector& value) {
        const double length = std::sqrt(dot(value, value));
        if (length <= kGeometryTolerance) {
            throw std::invalid_argument("Degenerate bend transition axis.");
        }
        return {value.x / length, value.y / length, value.z / length};
    }

    static Vector rotate(
        const Vector& value,
        const Vector& axis,
        double radians) {
        const double cosine = std::cos(radians);
        const double sine = std::sin(radians);
        const double parallel = dot(axis, value) * (1.0 - cosine);
        const Vector crossed = cross(axis, value);
        return {
            value.x * cosine + crossed.x * sine + axis.x * parallel,
            value.y * cosine + crossed.y * sine + axis.y * parallel,
            value.z * cosine + crossed.z * sine + axis.z * parallel};
    }

    static int roundVoxel(double value) {
        return static_cast<int>(std::lround(
            value + std::copysign(1e-7, value)));
    }

    static float euclideanDistance(
        const Point3D& lhs,
        const Point3D& rhs) noexcept {
        const double dx = static_cast<double>(rhs.x - lhs.x);
        const double dy = static_cast<double>(rhs.y - lhs.y);
        const double dz = static_cast<double>(rhs.z - lhs.z);
        return static_cast<float>(std::sqrt(
            dx * dx + dy * dy + dz * dz));
    }

    static int straightStepCount(
        const Point3D& from,
        const Point3D& goal,
        const voxel_planner::BusbarPose& pose) noexcept {
        const int delta[3] = {
            goal.x - from.x,
            goal.y - from.y,
            goal.z - from.z};
        const int step[3] = {pose.tx, pose.ty, pose.tz};
        int count = -1;
        for (int axis = 0; axis < 3; ++axis) {
            if (step[axis] == 0) {
                if (delta[axis] != 0) {
                    return -1;
                }
                continue;
            }
            if (delta[axis] % step[axis] != 0) {
                return -1;
            }
            const int axisCount = delta[axis] / step[axis];
            if (axisCount <= 0 ||
                (count >= 0 && axisCount != count)) {
                return -1;
            }
            count = axisCount;
        }
        return count;
    }

    static BackwardDistanceField buildBackwardDistanceField(
        const VoxelGrid& map,
        std::size_t startVoxel,
        std::size_t goalVoxel,
        const Point3D& rawStart,
        const Point3D& rawGoal,
        float endpointImmunityRadius,
        const std::vector<Point3D>* topologyHint,
        float topologyHintTolerance) {
        if (map.voxelCount() >
            static_cast<std::size_t>(
                std::numeric_limits<std::uint32_t>::max())) {
            throw std::overflow_error(
                "Backward BFS requires voxel indices to fit in 32 bits.");
        }

        BackwardDistanceField field;
        field.width = map.width();
        field.height = map.height();
        field.depth = map.depth();
        const std::uint64_t maximumDistance =
            static_cast<std::uint64_t>(field.width - 1U) +
            static_cast<std::uint64_t>(field.height - 1U) +
            static_cast<std::uint64_t>(field.depth - 1U);
        const bool useCompactDistances = maximumDistance <
            BackwardDistanceField::kCompactUnreachable;

        const std::size_t plane =
            static_cast<std::size_t>(field.width) * field.height;
        const double corridorRadius =
            std::isfinite(topologyHintTolerance)
                ? std::max(0.0F, topologyHintTolerance)
                : 0.0F;
        const double corridorRadiusSquared =
            corridorRadius * corridorRadius;

        if (useCompactDistances) {
            field.compactDistances.assign(
                map.voxelCount(),
                BackwardDistanceField::kCompactUnreachable);
            field.compactDistances[goalVoxel] = 0U;
        } else {
            field.wideDistances.assign(
                map.voxelCount(),
                BackwardDistanceField::kWideUnreachable);
            field.wideDistances[goalVoxel] = 0U;
        }

        std::vector<std::uint32_t> frontier;
        std::vector<std::uint32_t> nextFrontier;
        frontier.reserve(4096U);
        nextFrontier.reserve(4096U);
        frontier.push_back(static_cast<std::uint32_t>(goalVoxel));
        field.reachableVoxels = 1U;

        const auto withinTopologyCorridor =
            [&](std::size_t candidate) noexcept {
                if (topologyHint == nullptr || topologyHint->empty() ||
                    candidate == startVoxel || candidate == goalVoxel) {
                    return true;
                }
                const Point3D candidatePoint = pointFromVoxelIndex(
                    map,
                    candidate);
                for (const Point3D& hintPoint : *topologyHint) {
                    const double dx = static_cast<double>(
                        candidatePoint.x - hintPoint.x);
                    const double dy = static_cast<double>(
                        candidatePoint.y - hintPoint.y);
                    const double dz = static_cast<double>(
                        candidatePoint.z - hintPoint.z);
                    if (dx * dx + dy * dy + dz * dz <=
                        corridorRadiusSquared) {
                        return true;
                    }
                }
                return false;
            };
        std::uint32_t nextDistance = 1U;
        while (!frontier.empty()) {
            nextFrontier.clear();
            const auto tryVisit = [&](std::size_t candidate) {
                const bool unvisited = useCompactDistances
                    ? field.compactDistances[candidate] ==
                        BackwardDistanceField::kCompactUnreachable
                    : field.wideDistances[candidate] ==
                        BackwardDistanceField::kWideUnreachable;
                if (!unvisited) {
                    return;
                }
                if (!withinTopologyCorridor(candidate)) {
                    return;
                }

                if (candidate != startVoxel && candidate != goalVoxel &&
                    map.isObstacle(candidate)) {
                    const Point3D candidatePoint = pointFromVoxelIndex(
                        map,
                        candidate);
                    if (!isNearEndpoint(
                            candidatePoint,
                            rawStart,
                            rawGoal,
                            endpointImmunityRadius)) {
                        return;
                    }
                }

                if (useCompactDistances) {
                    field.compactDistances[candidate] =
                        static_cast<std::uint16_t>(nextDistance);
                } else {
                    field.wideDistances[candidate] = nextDistance;
                }
                nextFrontier.push_back(
                    static_cast<std::uint32_t>(candidate));
                ++field.reachableVoxels;
            };

            for (const std::uint32_t compactIndex : frontier) {
                const std::size_t voxelIndex = compactIndex;
                const std::size_t planeIndex = voxelIndex % plane;
                const std::uint32_t z = static_cast<std::uint32_t>(
                    voxelIndex / plane);
                const std::uint32_t y = static_cast<std::uint32_t>(
                    planeIndex / field.width);
                const std::uint32_t x = static_cast<std::uint32_t>(
                    planeIndex % field.width);
                if (x > 0U) {
                    tryVisit(voxelIndex - 1U);
                }
                if (x + 1U < field.width) {
                    tryVisit(voxelIndex + 1U);
                }
                if (y > 0U) {
                    tryVisit(voxelIndex - field.width);
                }
                if (y + 1U < field.height) {
                    tryVisit(voxelIndex + field.width);
                }
                if (z > 0U) {
                    tryVisit(voxelIndex - plane);
                }
                if (z + 1U < field.depth) {
                    tryVisit(voxelIndex + plane);
                }
            }
            frontier.swap(nextFrontier);
            ++nextDistance;
        }
        return field;
    }

    static ActionCatalog buildActionCatalog(
        const std::vector<voxel_planner::BusbarPose>& poses,
        const voxel_planner::PlannerConfig& config) {
        ActionCatalog catalog;
        catalog.straightSweeps.resize(poses.size());
        catalog.straightSweepBounds.resize(poses.size());
        catalog.straightLengths.resize(poses.size());
        catalog.bends.resize(poses.size());
        catalog.bendSweeps.resize(poses.size());
        catalog.bendSweepBounds.resize(poses.size());
        catalog.twists.resize(poses.size());
        catalog.twistSweeps.resize(poses.size());
        catalog.twistSweepBounds.resize(poses.size());
        catalog.poseEquivalenceClass.resize(poses.size());

        for (const voxel_planner::BusbarPose& pose : poses) {
            const float distance = static_cast<float>(std::sqrt(
                static_cast<double>(
                    pose.tx * pose.tx +
                    pose.ty * pose.ty +
                    pose.tz * pose.tz)));
            catalog.straightLengths[pose.poseId] = distance;
            catalog.straightSweeps[pose.poseId] =
                voxel_planner::generateStraightSweep(
                    pose,
                    config,
                    distance);
            catalog.straightSweepBounds[pose.poseId] = offsetBounds(
                catalog.straightSweeps[pose.poseId]);
        }

        for (const voxel_planner::BusbarPose& pose : poses) {
            std::uint32_t equivalenceClass = kInvalidPoseId;
            const Vector normal = unitNormal(pose);
            for (const voxel_planner::BusbarPose& previous : poses) {
                if (previous.poseId >= pose.poseId) {
                    continue;
                }
                if (previous.tx == pose.tx &&
                    previous.ty == pose.ty &&
                    previous.tz == pose.tz &&
                    std::abs(dot(normal, unitNormal(previous))) >=
                        1.0 - kGeometryTolerance) {
                    equivalenceClass =
                        catalog.poseEquivalenceClass[previous.poseId];
                    break;
                }
            }
            if (equivalenceClass == kInvalidPoseId) {
                equivalenceClass = static_cast<std::uint32_t>(
                    catalog.poseEquivalenceClassCount++);
            }
            catalog.poseEquivalenceClass[pose.poseId] = equivalenceClass;

            for (int direction : {-1, 1}) {
                int targetRoll =
                    pose.roll_angle_deg + direction * config.angle_step_deg;
                targetRoll %= 360;
                if (targetRoll < 0) {
                    targetRoll += 360;
                }
                for (const voxel_planner::BusbarPose& candidate : poses) {
                    if (candidate.tx == pose.tx &&
                        candidate.ty == pose.ty &&
                        candidate.tz == pose.tz &&
                        candidate.roll_angle_deg == targetRoll &&
                        candidate.poseId != pose.poseId) {
                        catalog.twists[pose.poseId].push_back({
                            candidate.poseId,
                            static_cast<float>(
                                direction * config.angle_step_deg *
                                kPi / 180.0)});
                        break;
                    }
                }
            }
        }

        const Vector globalAxes[3] = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}};
        for (const voxel_planner::BusbarPose& from : poses) {
            const Vector fromTangent = unitTangent(from);
            const Vector fromNormal = unitNormal(from);
            const Vector fromBinormal = normalized(cross(
                fromTangent,
                fromNormal));

            for (const Vector& globalAxis : globalAxes) {
                const double parallel = dot(globalAxis, fromTangent);
                const Vector projectedAxis{
                    globalAxis.x - parallel * fromTangent.x,
                    globalAxis.y - parallel * fromTangent.y,
                    globalAxis.z - parallel * fromTangent.z};
                const double projectedLength = std::sqrt(dot(
                    projectedAxis,
                    projectedAxis));
                if (projectedLength <= kGeometryTolerance) {
                    continue;
                }
                const Vector baseAxis{
                    projectedAxis.x / projectedLength,
                    projectedAxis.y / projectedLength,
                    projectedAxis.z / projectedLength};

                for (int direction : {-1, 1}) {
                    const Vector axis{
                        direction * baseAxis.x,
                        direction * baseAxis.y,
                        direction * baseAxis.z};
                    for (int bendDegrees = config.angle_step_deg;
                         bendDegrees <= 90;
                         bendDegrees += config.angle_step_deg) {
                        const double angle = bendDegrees * kPi / 180.0;
                        const Vector expectedTangent = rotate(
                            fromTangent,
                            axis,
                            angle);
                        const Vector expectedNormal = rotate(
                            fromNormal,
                            axis,
                            angle);

                        std::uint32_t toPoseId = kInvalidPoseId;
                        double bestTangentAlignment =
                            -std::numeric_limits<double>::infinity();
                        double bestNormalAlignment =
                            -std::numeric_limits<double>::infinity();
                        for (const voxel_planner::BusbarPose& candidate :
                             poses) {
                            const double tangentAlignment = dot(
                                expectedTangent,
                                unitTangent(candidate));
                            const double normalAlignment = dot(
                                expectedNormal,
                                unitNormal(candidate));
                            if (tangentAlignment >
                                    bestTangentAlignment +
                                        kGeometryTolerance ||
                                (std::abs(
                                     tangentAlignment -
                                     bestTangentAlignment) <=
                                     kGeometryTolerance &&
                                 normalAlignment > bestNormalAlignment)) {
                                bestTangentAlignment = tangentAlignment;
                                bestNormalAlignment = normalAlignment;
                                toPoseId = candidate.poseId;
                            }
                        }
                        if (toPoseId == kInvalidPoseId) {
                            continue;
                        }

                        const double flatAlignment = std::abs(dot(
                            axis,
                            fromNormal));
                        const double verticalAlignment = std::abs(dot(
                            axis,
                            fromBinormal));
                        const float radius =
                            flatAlignment >= verticalAlignment
                                ? config.flat_bend_factor *
                                      config.busbar_thickness
                                : config.vertical_bend_factor *
                                      config.busbar_width;
                        if (radius <= 0.0F) {
                            throw std::invalid_argument(
                                "Derived bend radius must be positive.");
                        }

                        const Vector radialStart = cross(
                            fromTangent,
                            axis);
                        const Vector radialEnd = rotate(
                            radialStart,
                            axis,
                            angle);
                        const double physicalLength = radius * angle;
                        const int dx = roundVoxel(
                            radius * (radialEnd.x - radialStart.x));
                        const int dy = roundVoxel(
                            radius * (radialEnd.y - radialStart.y));
                        const int dz = roundVoxel(
                            radius * (radialEnd.z - radialStart.z));
                        if (dx == 0 && dy == 0 && dz == 0) {
                            continue;
                        }

                        catalog.bends[from.poseId].push_back({
                            toPoseId,
                            dx,
                            dy,
                            dz,
                            radius,
                            static_cast<float>(angle),
                            static_cast<float>(physicalLength),
                            static_cast<float>(axis.x),
                            static_cast<float>(axis.y),
                            static_cast<float>(axis.z),
                            bendDegrees});
                        const double latticeDistance = std::sqrt(
                            static_cast<double>(
                                dx * dx + dy * dy + dz * dz));
                        catalog.heuristicScale = std::min(
                            catalog.heuristicScale,
                            static_cast<float>(
                                radius * angle / latticeDistance));
                    }
                }
            }
        }
        for (const voxel_planner::BusbarPose& pose : poses) {
            catalog.bendSweeps[pose.poseId].resize(
                catalog.bends[pose.poseId].size());
            catalog.bendSweepBounds[pose.poseId].resize(
                catalog.bends[pose.poseId].size());
            catalog.twistSweeps[pose.poseId].resize(
                catalog.twists[pose.poseId].size());
            catalog.twistSweepBounds[pose.poseId].resize(
                catalog.twists[pose.poseId].size());
        }
        return catalog;
    }

    static SearchOutcome runSingleSearch(
        const VoxelGrid& map,
        std::uint64_t startKey,
        std::uint64_t goalKey,
        const std::vector<voxel_planner::BusbarPose>& poses,
        const voxel_planner::PlannerConfig& config,
        ActionCatalog& actions,
        const BackwardDistanceField& backwardDistance,
        const CostPenaltyMap& penalties,
        const Point3D& rawStart,
        const Point3D& rawGoal,
        float endpointImmunityRadius) {
        SearchOutcome outcome;
        std::unordered_map<std::uint64_t, StateRecord> records;
        records.reserve(131072U);
        std::unordered_map<std::uint64_t, DominanceRecord> dominanceRecords;
        dominanceRecords.reserve(65536U);
        std::vector<Node> openStorage;
        openStorage.reserve(65536U);
        std::priority_queue<Node, std::vector<Node>, NodeCompare> openList(
            NodeCompare{},
            std::move(openStorage));

        const std::size_t startVoxel = startKey / poses.size();
        const std::uint32_t startPoseId = static_cast<std::uint32_t>(
            startKey % poses.size());
        const std::size_t goalVoxel = goalKey / poses.size();
        const std::uint32_t goalPoseId = static_cast<std::uint32_t>(
            goalKey % poses.size());
        const Vector goalTangent = unitTangent(poses[goalPoseId]);
        const Point3D goalPoint = pointFromVoxelIndex(map, goalVoxel);
        const Point3D startPoint = pointFromVoxelIndex(map, startVoxel);
        outcome.backwardReachableVoxels =
            backwardDistance.reachableVoxels;
        outcome.backwardStartDistance =
            backwardDistance.distance(startVoxel);
        outcome.closestBackwardDistance =
            outcome.backwardStartDistance;
        const float startH = heuristic(
            startVoxel,
            startPoint,
            goalPoint,
            unitTangent(poses[startPoseId]),
            goalTangent,
            actions.heuristicScale,
            backwardDistance);
        records[startKey] = {0.0F, kNoParent, kNoAction, kOpen};
        const std::uint64_t startDominanceKey = stateKey(
            startVoxel,
            actions.poseEquivalenceClass[startPoseId],
            actions.poseEquivalenceClassCount);
        dominanceRecords[startDominanceKey] = {0.0F, startKey};
        openList.push({
            startKey,
            static_cast<std::uint32_t>(startVoxel),
            startPoseId,
            0.0F,
            kHeuristicWeight * startH,
            startH});

        std::size_t maximumCandidateCount = 1U;
        for (std::size_t poseId = 0U;
             poseId < actions.bends.size();
             ++poseId) {
            maximumCandidateCount = std::max(
                maximumCandidateCount,
                actions.bends[poseId].size() +
                    actions.twists[poseId].size() + 1U);
        }
        std::vector<NeighborCandidate> candidates(maximumCandidateCount);

        std::size_t expansions = 0U;
        const std::size_t maxExpansions = penalties.values.empty()
            ? kMaxExpansionsPerPath
            : kMaxExpansionsPerPenalizedPath;
        std::uint64_t bestGoalKey = kNoParent;
        float bestGoalCost = std::numeric_limits<float>::infinity();
        std::uint64_t coarseGoalKey = kNoParent;
        float coarseGoalDistance = std::numeric_limits<float>::infinity();
        float coarseGoalCost = std::numeric_limits<float>::infinity();
        while (!openList.empty() && expansions < maxExpansions) {
            const Node current = openList.top();
            openList.pop();
            auto currentRecord = records.find(current.stateKey);
            if (currentRecord == records.end() ||
                currentRecord->second.status != kOpen ||
                current.g > currentRecord->second.g + 1e-5F) {
                continue;
            }
            const std::uint64_t currentDominanceKey = stateKey(
                current.voxelIndex,
                actions.poseEquivalenceClass[current.poseId],
                actions.poseEquivalenceClassCount);
            const auto currentDominance = dominanceRecords.find(
                currentDominanceKey);
            if (currentDominance == dominanceRecords.end() ||
                currentDominance->second.stateKey != current.stateKey ||
                current.g > currentDominance->second.g +
                    kDominanceCostTolerance) {
                continue;
            }
            currentRecord->second.status = kClosed;
            ++expansions;

            const Point3D currentPoint = pointFromVoxelIndex(
                map,
                current.voxelIndex);
            const double currentDx = static_cast<double>(
                goalPoint.x - currentPoint.x);
            const double currentDy = static_cast<double>(
                goalPoint.y - currentPoint.y);
            const double currentDz = static_cast<double>(
                goalPoint.z - currentPoint.z);
            const float currentGoalDistance = static_cast<float>(std::sqrt(
                currentDx * currentDx +
                currentDy * currentDy +
                currentDz * currentDz));
            outcome.closestGoalDistance = std::min(
                outcome.closestGoalDistance,
                currentGoalDistance);
            outcome.closestBackwardDistance = std::min(
                outcome.closestBackwardDistance,
                backwardDistance.distance(current.voxelIndex));

            const voxel_planner::BusbarPose& currentPose =
                poses[current.poseId];

            if (current.voxelIndex == goalVoxel) {
                if (current.g < bestGoalCost) {
                    bestGoalKey = current.stateKey;
                    bestGoalCost = current.g;
                }
                if (penalties.values.empty()) {
                    reconstructPath(
                        map,
                        records,
                        startKey,
                        current.stateKey,
                        poses,
                        actions,
                        outcome);
                    outcome.path.cost = current.g;
                    outcome.error_code = ErrorCode::NONE;
                    return outcome;
                }
                continue;
            }
            if (currentGoalDistance < 2.0F &&
                current.g + currentGoalDistance < coarseGoalCost) {
                coarseGoalKey = current.stateKey;
                coarseGoalDistance = currentGoalDistance;
                coarseGoalCost = current.g + currentGoalDistance;
            }

            if (penalties.values.empty() &&
                dot(unitTangent(currentPose), goalTangent) > 0.99) {
                const int terminalSteps = straightStepCount(
                    currentPoint,
                    goalPoint,
                    currentPose);
                bool connectorValid = terminalSteps > 0;
                Point3D connectorAnchor = currentPoint;
                for (int step = 0;
                     connectorValid && step < terminalSteps;
                     ++step) {
                    const Point3D connectorTarget{
                        connectorAnchor.x + currentPose.tx,
                        connectorAnchor.y + currentPose.ty,
                        connectorAnchor.z + currentPose.tz};
                    const std::size_t connectorVoxel = map.index(
                        static_cast<std::uint32_t>(connectorTarget.x),
                        static_cast<std::uint32_t>(connectorTarget.y),
                        static_cast<std::uint32_t>(connectorTarget.z));
                    const bool endpointImmune = isNearEndpoint(
                            connectorAnchor,
                            rawStart,
                            rawGoal,
                            endpointImmunityRadius) ||
                        isNearEndpoint(
                            connectorTarget,
                            rawStart,
                            rawGoal,
                            endpointImmunityRadius);
                    if (!endpointImmune &&
                        (map.isObstacle(connectorVoxel) ||
                         !map.isPoseAllowedForSearch(
                             connectorVoxel,
                             current.poseId) ||
                         !isSweepCollisionFreeForSearch(
                             map,
                             connectorAnchor,
                             actions.straightSweeps[current.poseId],
                             actions.straightSweepBounds[
                                 current.poseId]))) {
                        connectorValid = false;
                        break;
                    }
                    connectorAnchor = connectorTarget;
                }
                if (connectorValid) {
                    std::uint64_t parentKey = current.stateKey;
                    float connectorG = current.g;
                    Point3D connectorPoint = currentPoint;
                    for (int step = 0; step < terminalSteps; ++step) {
                        connectorPoint.x += currentPose.tx;
                        connectorPoint.y += currentPose.ty;
                        connectorPoint.z += currentPose.tz;
                        const std::size_t connectorVoxel = map.index(
                            static_cast<std::uint32_t>(connectorPoint.x),
                            static_cast<std::uint32_t>(connectorPoint.y),
                            static_cast<std::uint32_t>(connectorPoint.z));
                        const std::uint64_t connectorKey = stateKey(
                            connectorVoxel,
                            current.poseId,
                            poses.size());
                        connectorG +=
                            actions.straightLengths[current.poseId];
                        const bool connectorPenaltyExempt = isNearEndpoint(
                            connectorPoint,
                            rawStart,
                            rawGoal,
                            endpointImmunityRadius);
                        connectorG += connectorPenaltyExempt
                            ? 0.0F
                            : penalties.get(connectorVoxel);
                        records[connectorKey] = {
                            connectorG,
                            parentKey,
                            kStraightAction,
                            kOpen};
                        parentKey = connectorKey;
                    }
                    reconstructPath(
                        map,
                        records,
                        startKey,
                        parentKey,
                        poses,
                        actions,
                        outcome);
                    outcome.path.cost = connectorG;
                    outcome.error_code = ErrorCode::NONE;
                    return outcome;
                }
            }

            const std::vector<BendAction>& poseBends =
                actions.bends[current.poseId];
            const std::vector<TwistAction>& poseTwists =
                actions.twists[current.poseId];
            const std::size_t bendCount = poseBends.size();
            const std::size_t candidateCount =
                bendCount + poseTwists.size() + 1U;
            const bool currentNearEndpoint = isNearEndpoint(
                currentPoint,
                rawStart,
                rawGoal,
                endpointImmunityRadius);
            std::atomic<bool> allocationFailed(false);

            for (std::size_t index = 0U;
                 index < candidateCount;
                 ++index) {
                NeighborCandidate prepared;
                const bool straight = index == 0U;
                const bool twist = index > bendCount;
                const std::size_t twistIndex = twist
                    ? index - bendCount - 1U
                    : 0U;
                const BendAction* bend = !straight && !twist
                    ? &poseBends[index - 1U]
                    : nullptr;
                const TwistAction* twistAction = twist
                    ? &poseTwists[twistIndex]
                    : nullptr;
                prepared.targetPoseId = straight
                    ? current.poseId
                    : (twist
                        ? twistAction->toPoseId
                        : bend->toPoseId);
                prepared.transitionAction = straight
                    ? kStraightAction
                    : (twist
                        ? kTwistActionBase |
                            static_cast<std::uint32_t>(twistIndex)
                        : static_cast<std::uint32_t>(index - 1U));
                prepared.physicalLength = straight
                    ? actions.straightLengths[current.poseId]
                    : (twist ? 0.0F : bend->physicalLength);
                prepared.targetPoint = {
                    currentPoint.x +
                        (straight ? currentPose.tx : (twist ? 0 : bend->dx)),
                    currentPoint.y +
                        (straight ? currentPose.ty : (twist ? 0 : bend->dy)),
                    currentPoint.z +
                        (straight ? currentPose.tz : (twist ? 0 : bend->dz))};
                if (!map.isValid(
                        prepared.targetPoint.x,
                        prepared.targetPoint.y,
                        prepared.targetPoint.z)) {
                    candidates[index] = prepared;
                    continue;
                }
                prepared.targetVoxel = map.index(
                    static_cast<std::uint32_t>(prepared.targetPoint.x),
                    static_cast<std::uint32_t>(prepared.targetPoint.y),
                    static_cast<std::uint32_t>(prepared.targetPoint.z));
                prepared.endpointImmune = currentNearEndpoint ||
                    isNearEndpoint(
                        prepared.targetPoint,
                        rawStart,
                        rawGoal,
                        endpointImmunityRadius);

                const std::uint64_t targetKey = stateKey(
                    prepared.targetVoxel,
                    prepared.targetPoseId,
                    poses.size());
                const float tentativeG =
                    current.g + prepared.physicalLength +
                    (prepared.targetVoxel == current.voxelIndex
                        ? 0.0F
                        : penalties.get(prepared.targetVoxel));
                const auto existing = records.find(targetKey);
                if (existing != records.end()) {
                    if (existing->second.status == kClosed ||
                        tentativeG + 1e-5F >= existing->second.g) {
                        candidates[index] = prepared;
                        continue;
                    }
                }
                const std::uint64_t targetDominanceKey = stateKey(
                    prepared.targetVoxel,
                    actions.poseEquivalenceClass[prepared.targetPoseId],
                    actions.poseEquivalenceClassCount);
                const auto dominant = dominanceRecords.find(
                    targetDominanceKey);
                if (dominant != dominanceRecords.end() &&
                    tentativeG + kDominanceCostTolerance >=
                        dominant->second.g) {
                    candidates[index] = prepared;
                    continue;
                }
                bool duplicate = false;
                for (std::size_t previous = 0U;
                     previous < index;
                     ++previous) {
                    if (candidates[previous].evaluate &&
                        candidates[previous].targetVoxel ==
                            prepared.targetVoxel &&
                        actions.poseEquivalenceClass[
                            candidates[previous].targetPoseId] ==
                        actions.poseEquivalenceClass[
                            prepared.targetPoseId]) {
                        duplicate = true;
                        break;
                    }
                }
                prepared.evaluate = !duplicate;
                candidates[index] = prepared;
            }

#pragma omp parallel for schedule(static) num_threads(kNeighborThreadCount) \
    if(candidateCount >= kNeighborParallelThreshold)
            for (std::ptrdiff_t candidateIndex = 0;
                 candidateIndex <
                     static_cast<std::ptrdiff_t>(candidateCount);
                 ++candidateIndex) {
                const std::size_t index =
                    static_cast<std::size_t>(candidateIndex);
                NeighborCandidate evaluated = candidates[index];
                if (!evaluated.evaluate) {
                    continue;
                }
                const bool straight = index == 0U;
                const bool twist = index > bendCount;
                const std::size_t twistIndex = twist
                    ? index - bendCount - 1U
                    : 0U;
                const BendAction* bend = !straight && !twist
                    ? &poseBends[index - 1U]
                    : nullptr;
                const TwistAction* twistAction = twist
                    ? &poseTwists[twistIndex]
                    : nullptr;
                // O(1) center fail precedes footprint and sweep work.
                if (!evaluated.endpointImmune &&
                    (map.isObstacle(evaluated.targetVoxel) ||
                     !map.isPoseAllowedForSearch(
                         evaluated.targetVoxel,
                         evaluated.targetPoseId))) {
                    candidates[index] = evaluated;
                    continue;
                }

                const std::vector<voxel_planner::VoxelOffset>* sweep = nullptr;
                const voxel_planner::VoxelOffsetBounds* sweepBounds = nullptr;
                if (straight) {
                    sweep = &actions.straightSweeps[current.poseId];
                    sweepBounds =
                        &actions.straightSweepBounds[current.poseId];
                } else if (!twist) {
                    std::vector<voxel_planner::VoxelOffset>& cachedSweep =
                        actions.bendSweeps[current.poseId][index - 1U];
                    voxel_planner::VoxelOffsetBounds& cachedBounds =
                        actions.bendSweepBounds[current.poseId][index - 1U];
                    if (cachedSweep.empty()) {
                        try {
                            cachedSweep =
                                voxel_planner::generateExplicitBendSweep(
                                    currentPose,
                                    config,
                                    bend->radius,
                                    bend->angleRadians,
                                    bend->axisX,
                                    bend->axisY,
                                    bend->axisZ);
                            cachedBounds = offsetBounds(cachedSweep);
                        } catch (const std::bad_alloc&) {
                            allocationFailed.store(
                                true,
                                std::memory_order_relaxed);
                        } catch (const std::invalid_argument&) {
                        }
                    }
                    sweep = &cachedSweep;
                    sweepBounds = &cachedBounds;
                } else {
                    std::vector<voxel_planner::VoxelOffset>& cachedSweep =
                        actions.twistSweeps[current.poseId][twistIndex];
                    voxel_planner::VoxelOffsetBounds& cachedBounds =
                        actions.twistSweepBounds[current.poseId][twistIndex];
                    if (cachedSweep.empty()) {
                        try {
                            cachedSweep = voxel_planner::generateTwistSweep(
                                currentPose,
                                poses[twistAction->toPoseId],
                                config,
                                twistAction->angleRadians);
                            cachedBounds = offsetBounds(cachedSweep);
                        } catch (const std::bad_alloc&) {
                            allocationFailed.store(
                                true,
                                std::memory_order_relaxed);
                        } catch (const std::invalid_argument&) {
                        }
                    }
                    sweep = &cachedSweep;
                    sweepBounds = &cachedBounds;
                }

                if (allocationFailed.load(std::memory_order_relaxed) ||
                    sweep->empty() ||
                    sweepBounds == nullptr ||
                    (!evaluated.endpointImmune &&
                     !isSweepCollisionFreeForSearch(
                         map,
                         currentPoint,
                         *sweep,
                         *sweepBounds))) {
                    candidates[index] = evaluated;
                    continue;
                }
                evaluated.valid = true;
                candidates[index] = evaluated;
            }

            if (allocationFailed.load(std::memory_order_relaxed)) {
                throw std::bad_alloc();
            }
            for (std::size_t candidateIndex = 0U;
                 candidateIndex < candidateCount;
                 ++candidateIndex) {
                const NeighborCandidate& candidate =
                    candidates[candidateIndex];
                if (!candidate.valid) {
                    continue;
                }
                relaxState(
                    current,
                    candidate.targetVoxel,
                    candidate.targetPoseId,
                    candidate.physicalLength,
                    candidate.targetPoint,
                    goalPoint,
                    unitTangent(poses[candidate.targetPoseId]),
                    goalTangent,
                    actions.heuristicScale,
                    backwardDistance,
                    penalties,
                    poses.size(),
                    candidate.transitionAction,
                    records,
                    actions.poseEquivalenceClass,
                    actions.poseEquivalenceClassCount,
                    dominanceRecords,
                    openList);
                const float candidateGoalDistance = euclideanDistance(
                    candidate.targetPoint,
                    goalPoint);
                if (candidate.targetVoxel == goalVoxel) {
                    const std::uint64_t reachedGoalKey = stateKey(
                        candidate.targetVoxel,
                        candidate.targetPoseId,
                        poses.size());
                    const auto reachedGoal = records.find(reachedGoalKey);
                    if (reachedGoal != records.end()) {
                        if (reachedGoal->second.g < bestGoalCost) {
                            bestGoalKey = reachedGoalKey;
                            bestGoalCost = reachedGoal->second.g;
                        }
                        if (penalties.values.empty()) {
                            reconstructPath(
                                map,
                                records,
                                startKey,
                                reachedGoalKey,
                                poses,
                                actions,
                                outcome);
                            outcome.path.cost = reachedGoal->second.g;
                            outcome.error_code = ErrorCode::NONE;
                            return outcome;
                        }
                    }
                } else if (candidateGoalDistance < 2.0F) {
                    const std::uint64_t reachedCoarseGoalKey = stateKey(
                        candidate.targetVoxel,
                        candidate.targetPoseId,
                        poses.size());
                    const auto reachedCoarseGoal = records.find(
                        reachedCoarseGoalKey);
                    if (reachedCoarseGoal != records.end() &&
                        reachedCoarseGoal->second.g + candidateGoalDistance <
                            coarseGoalCost) {
                        coarseGoalKey = reachedCoarseGoalKey;
                        coarseGoalDistance = candidateGoalDistance;
                        coarseGoalCost = reachedCoarseGoal->second.g +
                            candidateGoalDistance;
                    }
                }
            }
        }
        if (bestGoalKey != kNoParent) {
            reconstructPath(
                map,
                records,
                startKey,
                bestGoalKey,
                poses,
                actions,
                outcome);
            outcome.path.cost = bestGoalCost;
            outcome.error_code = ErrorCode::NONE;
            outcome.expansions = expansions;
            outcome.generatedStates = records.size();
            outcome.remainingOpenStates = openList.size();
            return outcome;
        }
        if (coarseGoalKey != kNoParent) {
            reconstructPath(
                map,
                records,
                startKey,
                coarseGoalKey,
                poses,
                actions,
                outcome);
            outcome.path.path.push_back(goalPoint);
            appendConditionalPose(
                map,
                goalPoint,
                unitNormal(poses[goalPoseId]),
                goalTangent,
                outcome.path.path.size() - 1U,
                outcome.path.pose_description);
            outcome.path.cost = coarseGoalCost;
            outcome.error_code = ErrorCode::NONE;
            outcome.expansions = expansions;
            outcome.generatedStates = records.size();
            outcome.remainingOpenStates = openList.size();
            outcome.closestGoalDistance = std::min(
                outcome.closestGoalDistance,
                coarseGoalDistance);
            return outcome;
        }
        if (expansions >= maxExpansions) {
            outcome.error_code = ErrorCode::COMPUTATION_LIMIT_EXCEEDED;
        }
        outcome.expansions = expansions;
        outcome.generatedStates = records.size();
        outcome.remainingOpenStates = openList.size();
        for (const voxel_planner::BusbarPose& pose : poses) {
            const std::uint64_t candidateGoalKey = stateKey(
                goalVoxel,
                pose.poseId,
                poses.size());
            if (records.find(candidateGoalKey) != records.end()) {
                outcome.bestGeneratedGoalAlignment = std::max(
                    outcome.bestGeneratedGoalAlignment,
                    static_cast<float>(dot(
                        unitTangent(pose),
                        goalTangent)));
            }
        }
        return outcome;
    }

    template <typename Queue>
    static void relaxState(
        const Node& current,
        std::size_t targetVoxel,
        std::uint32_t targetPoseId,
        float physicalLength,
        const Point3D& targetPoint,
        const Point3D& goalPoint,
        const Vector& targetTangent,
        const Vector& goalTangent,
        float heuristicScale,
        const BackwardDistanceField& backwardDistance,
        const CostPenaltyMap& penalties,
        std::size_t poseCount,
        std::uint32_t transitionAction,
        std::unordered_map<std::uint64_t, StateRecord>& records,
        const std::vector<std::uint32_t>& poseEquivalenceClass,
        std::size_t poseEquivalenceClassCount,
        std::unordered_map<std::uint64_t, DominanceRecord>& dominanceRecords,
        Queue& openList) {
        const std::uint64_t targetKey = stateKey(
            targetVoxel,
            targetPoseId,
            poseCount);

        const float reusePenalty = targetVoxel == current.voxelIndex
            ? 0.0F
            : penalties.get(targetVoxel);
        const float tentativeG = current.g + physicalLength + reusePenalty;
        auto existing = records.find(targetKey);
        if (existing != records.end()) {
            if (existing->second.status == kClosed ||
                tentativeG + 1e-5F >= existing->second.g) {
                return;
            }
        }
        const std::uint64_t targetDominanceKey = stateKey(
            targetVoxel,
            poseEquivalenceClass[targetPoseId],
            poseEquivalenceClassCount);
        const auto dominant = dominanceRecords.find(targetDominanceKey);
        if (dominant != dominanceRecords.end() &&
            tentativeG + kDominanceCostTolerance >= dominant->second.g) {
            return;
        }

        records[targetKey] = {
            tentativeG,
            current.stateKey,
            transitionAction,
            kOpen};
        dominanceRecords[targetDominanceKey] = {tentativeG, targetKey};
        const float h = heuristic(
            targetVoxel,
            targetPoint,
            goalPoint,
            targetTangent,
            goalTangent,
            heuristicScale,
            backwardDistance);
        openList.push({
            targetKey,
            static_cast<std::uint32_t>(targetVoxel),
            targetPoseId,
            tentativeG,
            tentativeG + kHeuristicWeight * h,
            h});
    }

    static std::uint64_t stateKey(
        std::size_t voxelIndex,
        std::uint32_t poseId,
        std::size_t poseCount) {
        if (poseCount == 0U ||
            voxelIndex >
                (std::numeric_limits<std::uint64_t>::max() - poseId) /
                    poseCount) {
            throw std::overflow_error("SE(3) state key overflow.");
        }
        return static_cast<std::uint64_t>(voxelIndex) * poseCount + poseId;
    }

    static Point3D pointFromVoxelIndex(
        const VoxelGrid& map,
        std::size_t voxelIndex) {
        const std::size_t plane =
            static_cast<std::size_t>(map.width()) * map.height();
        const int z = static_cast<int>(voxelIndex / plane);
        const std::size_t remainder = voxelIndex % plane;
        const int y = static_cast<int>(remainder / map.width());
        const int x = static_cast<int>(remainder % map.width());
        return {x, y, z};
    }

    static float heuristic(
        std::size_t voxelIndex,
        const Point3D& from,
        const Point3D& goal,
        const Vector& currentTangent,
        const Vector& goalTangent,
        float scale,
        const BackwardDistanceField& backwardDistance) {
        const double dx = static_cast<double>(goal.x - from.x);
        const double dy = static_cast<double>(goal.y - from.y);
        const double dz = static_cast<double>(goal.z - from.z);
        const float spatialDistance = static_cast<float>(
            std::sqrt(dx * dx + dy * dy + dz * dz));
        const std::uint32_t obstacleDistance =
            backwardDistance.distance(voxelIndex);
        const float guidedDistance = obstacleDistance ==
                BackwardDistanceField::kWideUnreachable
            ? spatialDistance
            : static_cast<float>(obstacleDistance);
        const float goalAlignment = static_cast<float>(std::max(
            -1.0,
            std::min(1.0, dot(currentTangent, goalTangent))));
        const float orientationPenalty =
            (1.0F - goalAlignment) * scale *
            kGoalOrientationPenaltyDistance;

        float travelHeadingPenalty = 0.0F;
        if (spatialDistance > 1e-5F) {
            float travelAlignment = -1.0F;
            const std::size_t plane =
                static_cast<std::size_t>(backwardDistance.width) *
                backwardDistance.height;
            const std::size_t planeIndex = voxelIndex % plane;
            const int z = static_cast<int>(voxelIndex / plane);
            const int y = static_cast<int>(
                planeIndex / backwardDistance.width);
            const int x = static_cast<int>(
                planeIndex % backwardDistance.width);
            if (obstacleDistance !=
                BackwardDistanceField::kWideUnreachable) {
                for (int dzStep = -1; dzStep <= 1; ++dzStep) {
                    for (int dyStep = -1; dyStep <= 1; ++dyStep) {
                        for (int dxStep = -1; dxStep <= 1; ++dxStep) {
                            if (dxStep == 0 && dyStep == 0 && dzStep == 0) {
                                continue;
                            }
                            const int nx = x + dxStep;
                            const int ny = y + dyStep;
                            const int nz = z + dzStep;
                            if (nx < 0 || ny < 0 || nz < 0 ||
                                nx >= static_cast<int>(
                                    backwardDistance.width) ||
                                ny >= static_cast<int>(
                                    backwardDistance.height) ||
                                nz >= static_cast<int>(
                                    backwardDistance.depth)) {
                                continue;
                            }
                            const std::size_t neighbor =
                                static_cast<std::size_t>(nx) +
                                static_cast<std::size_t>(ny) *
                                    backwardDistance.width +
                                static_cast<std::size_t>(nz) * plane;
                            if (backwardDistance.distance(neighbor) >=
                                obstacleDistance) {
                                continue;
                            }
                            const double directionLength = std::sqrt(
                                static_cast<double>(
                                    dxStep * dxStep +
                                    dyStep * dyStep +
                                    dzStep * dzStep));
                            const Vector direction{
                                dxStep / directionLength,
                                dyStep / directionLength,
                                dzStep / directionLength};
                            travelAlignment = std::max(
                                travelAlignment,
                                static_cast<float>(dot(
                                    currentTangent,
                                    direction)));
                        }
                    }
                }
            }
            if (travelAlignment < -0.5F) {
                const Vector directHeading{
                    dx / spatialDistance,
                    dy / spatialDistance,
                    dz / spatialDistance};
                travelAlignment = static_cast<float>(dot(
                    currentTangent,
                    directHeading));
            }
            travelAlignment = std::max(
                -1.0F,
                std::min(1.0F, travelAlignment));
            travelHeadingPenalty =
                (1.0F - travelAlignment) * scale *
                kTravelHeadingPenaltyDistance;
        }
        return scale * guidedDistance +
            orientationPenalty +
            travelHeadingPenalty;
    }

    static void reconstructPath(
        const VoxelGrid& map,
        const std::unordered_map<std::uint64_t, StateRecord>& records,
        std::uint64_t startKey,
        std::uint64_t goalKey,
        const std::vector<voxel_planner::BusbarPose>& poses,
        const ActionCatalog& actions,
        SearchOutcome& outcome) {
        std::vector<std::uint64_t> reversed;
        std::uint64_t current = goalKey;
        while (true) {
            reversed.push_back(current);
            if (current == startKey) {
                break;
            }
            const auto record = records.find(current);
            if (record == records.end() || record->second.parent == kNoParent) {
                outcome.path.path.clear();
                outcome.path.pose_description.clear();
                outcome.statePath.clear();
                return;
            }
            current = record->second.parent;
        }
        outcome.statePath.assign(reversed.rbegin(), reversed.rend());

        for (std::size_t stateIndex = 0U;
             stateIndex < outcome.statePath.size();
             ++stateIndex) {
            const std::uint64_t key = outcome.statePath[stateIndex];
            const std::size_t voxelIndex = key / poses.size();
            const std::uint32_t poseId = static_cast<std::uint32_t>(
                key % poses.size());
            const Point3D point = pointFromVoxelIndex(map, voxelIndex);
            if (stateIndex == 0U) {
                outcome.path.path.push_back(point);
                appendConditionalPose(
                    map,
                    point,
                    unitNormal(poses[poseId]),
                    unitTangent(poses[poseId]),
                    0U,
                    outcome.path.pose_description);
            } else {
                const std::uint64_t previousKey =
                    outcome.statePath[stateIndex - 1U];
                const auto transitionRecord = records.find(key);
                if (transitionRecord == records.end()) {
                    outcome.path.path.clear();
                    outcome.path.pose_description.clear();
                    outcome.statePath.clear();
                    return;
                }
                appendTransitionCenterline(
                    map,
                    previousKey,
                    key,
                    transitionRecord->second.parentAction,
                    poses,
                    actions,
                    outcome.path.path,
                    outcome.path.pose_description);
            }
        }
    }

    static void appendConditionalPose(
        const VoxelGrid& map,
        const Point3D& point,
        const Vector& normal,
        const Vector& tangent,
        std::size_t waypointIndex,
        std::vector<Pose>& poseDescription) {
        if (map.isValid(point.x, point.y, point.z) &&
            map.getState(point.x, point.y, point.z) ==
                VoxelState::POSE_CONDITIONAL) {
            poseDescription.push_back({
                {normal.x, normal.y, normal.z},
                {tangent.x, tangent.y, tangent.z},
                waypointIndex});
        }
    }

    static void appendPoint26WithPose(
        const VoxelGrid& map,
        std::vector<Point3D>& path,
        const Point3D& target,
        const Vector& normal,
        const Vector& tangent,
        std::vector<Pose>& poseDescription) {
        const std::size_t firstNewPoint = path.size();
        appendPoint26(path, target);
        for (std::size_t index = firstNewPoint;
             index < path.size();
             ++index) {
            appendConditionalPose(
                map,
                path[index],
                normal,
                tangent,
                index,
                poseDescription);
        }
    }

    static void appendTransitionCenterline(
        const VoxelGrid& map,
        std::uint64_t fromKey,
        std::uint64_t toKey,
        std::uint32_t transitionAction,
        const std::vector<voxel_planner::BusbarPose>& poses,
        const ActionCatalog& actions,
        std::vector<Point3D>& path,
        std::vector<Pose>& poseDescription) {
        const std::size_t fromVoxel = fromKey / poses.size();
        const std::size_t toVoxel = toKey / poses.size();
        const std::uint32_t fromPoseId = static_cast<std::uint32_t>(
            fromKey % poses.size());
        const std::uint32_t toPoseId = static_cast<std::uint32_t>(
            toKey % poses.size());
        const Point3D from = pointFromVoxelIndex(map, fromVoxel);
        const Point3D to = pointFromVoxelIndex(map, toVoxel);
        if (transitionAction == kStraightAction) {
            appendPoint26WithPose(
                map,
                path,
                to,
                unitNormal(poses[toPoseId]),
                unitTangent(poses[toPoseId]),
                poseDescription);
            return;
        }

        if (isTwistAction(transitionAction)) {
            const std::size_t twistIndex = twistActionIndex(
                transitionAction);
            if (twistIndex < actions.twists[fromPoseId].size() &&
                !path.empty()) {
                appendConditionalPose(
                    map,
                    to,
                    unitNormal(poses[toPoseId]),
                    unitTangent(poses[toPoseId]),
                    path.size() - 1U,
                    poseDescription);
            }
            return;
        }

        if (transitionAction >= actions.bends[fromPoseId].size()) {
            appendPoint26WithPose(
                map,
                path,
                to,
                unitNormal(poses[toPoseId]),
                unitTangent(poses[toPoseId]),
                poseDescription);
            return;
        }
        const BendAction& selected =
            actions.bends[fromPoseId][transitionAction];

        const Vector startTangent = unitTangent(poses[fromPoseId]);
        const Vector startNormal = unitNormal(poses[fromPoseId]);
        const Vector endNormal = unitNormal(poses[toPoseId]);
        const Vector endTangent = unitTangent(poses[toPoseId]);
        const Vector axis{
            selected.axisX,
            selected.axisY,
            selected.axisZ};
        const Vector radialStart = cross(startTangent, axis);
        const int samples = std::max(
            1,
            static_cast<int>(std::ceil(
                selected.physicalLength / 0.49F)));
        for (int sample = 1; sample <= samples; ++sample) {
            const double angle =
                selected.angleRadians * sample / samples;
            const Vector radial = rotate(radialStart, axis, angle);
            const Vector tangent = rotate(startTangent, axis, angle);
            const Vector normal = normalized(rotate(
                startNormal,
                axis,
                angle));
            const Point3D samplePoint{
                from.x + roundVoxel(
                    selected.radius * (radial.x - radialStart.x)),
                from.y + roundVoxel(
                    selected.radius * (radial.y - radialStart.y)),
                from.z + roundVoxel(
                    selected.radius * (radial.z - radialStart.z))};
            appendPoint26WithPose(
                map,
                path,
                samplePoint,
                normal,
                tangent,
                poseDescription);
        }
        appendPoint26WithPose(
            map,
            path,
            to,
            endNormal,
            endTangent,
            poseDescription);
    }

    static void appendPoint26(
        std::vector<Point3D>& path,
        const Point3D& target) {
        if (path.empty()) {
            path.push_back(target);
            return;
        }
        const Point3D start = path.back();
        const int dx = target.x - start.x;
        const int dy = target.y - start.y;
        const int dz = target.z - start.z;
        const int steps = std::max({
            std::abs(dx),
            std::abs(dy),
            std::abs(dz)});
        if (steps == 0) {
            return;
        }
        for (int step = 1; step <= steps; ++step) {
            const Point3D point{
                start.x + roundVoxel(
                    static_cast<double>(dx) * step / steps),
                start.y + roundVoxel(
                    static_cast<double>(dy) * step / steps),
                start.z + roundVoxel(
                    static_cast<double>(dz) * step / steps)};
            if (!(path.back() == point)) {
                path.push_back(point);
            }
        }
    }

    static std::int64_t squaredDistance(
        const Point3D& lhs,
        const Point3D& rhs) {
        const std::int64_t dx =
            static_cast<std::int64_t>(lhs.x) - rhs.x;
        const std::int64_t dy =
            static_cast<std::int64_t>(lhs.y) - rhs.y;
        const std::int64_t dz =
            static_cast<std::int64_t>(lhs.z) - rhs.z;
        return dx * dx + dy * dy + dz * dz;
    }

    static double pathOverlapRatio(
        const VoxelGrid& map,
        const std::vector<Point3D>& candidate,
        const std::vector<PathResult>& acceptedPaths) {
        std::unordered_set<std::size_t> accepted;
        for (const PathResult& path : acceptedPaths) {
            for (const Point3D& point : path.path) {
                if (map.isValid(point.x, point.y, point.z)) {
                    accepted.insert(map.index(
                        static_cast<std::uint32_t>(point.x),
                        static_cast<std::uint32_t>(point.y),
                        static_cast<std::uint32_t>(point.z)));
                }
            }
        }
        std::unordered_set<std::size_t> candidateVoxels;
        std::size_t shared = 0U;
        for (const Point3D& point : candidate) {
            if (!map.isValid(point.x, point.y, point.z)) {
                continue;
            }
            const std::size_t index = map.index(
                static_cast<std::uint32_t>(point.x),
                static_cast<std::uint32_t>(point.y),
                static_cast<std::uint32_t>(point.z));
            if (candidateVoxels.insert(index).second &&
                accepted.find(index) != accepted.end()) {
                ++shared;
            }
        }
        return candidateVoxels.empty()
            ? 1.0
            : static_cast<double>(shared) /
                  static_cast<double>(candidateVoxels.size());
    }

    static void addSoftPathPenalty(
        const VoxelGrid& map,
        const std::vector<Point3D>& centerline,
        const Point3D& rawStart,
        const Point3D& rawGoal,
        float endpointProtectionRadius,
        CostPenaltyMap& penalties) {
        if (centerline.empty()) {
            return;
        }
        const double endpointRadiusSquared =
            static_cast<double>(endpointProtectionRadius) *
            endpointProtectionRadius;
        for (const Point3D& point : centerline) {
            for (int dz = -1; dz <= 1; ++dz) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        const Point3D candidate{
                            point.x + dx,
                            point.y + dy,
                            point.z + dz};
                        if (!map.isValid(
                                candidate.x,
                                candidate.y,
                                candidate.z) ||
                            static_cast<double>(squaredDistance(
                                candidate,
                                rawStart)) <= endpointRadiusSquared ||
                            static_cast<double>(squaredDistance(
                                candidate,
                                rawGoal)) <= endpointRadiusSquared) {
                            continue;
                        }
                        const std::size_t index = map.index(
                            static_cast<std::uint32_t>(candidate.x),
                            static_cast<std::uint32_t>(candidate.y),
                            static_cast<std::uint32_t>(candidate.z));
                        const float penalty = dx == 0 && dy == 0 && dz == 0
                            ? kCenterlineReusePenalty
                            : kNeighborReusePenalty;
                        penalties.values[index] += penalty;
                    }
                }
            }
        }
    }

};

} // namespace module3_astar
