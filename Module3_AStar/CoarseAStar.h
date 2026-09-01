#pragma once

#include "Module1_Main/Types.h"
#include "Module2_Morphology/PoseGenerator.h"
#include "Module2_Morphology/VoxelMorphologyEngine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace module3_astar {

/**
 * @brief Lightweight pose-aware A* over a 26-neighborhood voxel lattice.
 *
 * The search deliberately contains no Bend or Twist actions. Each edge is a
 * one-voxel translation, filtered first by the 0/1/2 pose semantics and then
 * by the cached physical swept-volume template.
 */
class CoarseAStar {
    struct PoseLut;

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
        const std::vector<Point3D>* topologyHint = nullptr,
        bool allowFullSearch = false,
        std::size_t maxExpansionsOverride = 0U) {
        PlanningResult result;
        const int boundedMaxPaths = std::min(
            maxPaths,
            kMaxPathsPerRequest);
        const bool pathLimitClamped = boundedMaxPaths != maxPaths;
        if (maxPaths <= 0 ||
            poses.empty() ||
            config.busbar_width <= 0.0F ||
            config.busbar_thickness <= 0.0F ||
            !std::isfinite(endpointImmunityRadius) ||
            endpointImmunityRadius < 0.0F ||
            !std::isfinite(config.max_overlap_ratio) ||
            config.max_overlap_ratio < 0.0F ||
            config.max_overlap_ratio > 1.0F ||
            (boundedMaxPaths > 1 &&
             (!std::isfinite(config.path_blocking_radius) ||
              config.path_blocking_radius <= 0.0F))) {
            result.error_code = ErrorCode::INVALID_ARGUMENT;
            result.message = "Invalid lightweight planner arguments.";
            return result;
        }
        if (!poseTableMatchesMap(map, poses)) {
            result.error_code = ErrorCode::POSE_MASK_UNAVAILABLE;
            result.message =
                "Pose table does not match the populated voxel pose masks.";
            return result;
        }

        PoseLut poseLut;
        try {
            poseLut = buildPoseLut(
                poses,
                config.angle_step_deg);
        } catch (const std::bad_alloc&) {
            result.error_code = ErrorCode::MEMORY_ALLOCATION_FAILED;
            result.message =
                "Memory allocation failed while building pose LUTs.";
            return result;
        } catch (const std::exception& error) {
            result.error_code = ErrorCode::INVALID_ARGUMENT;
            result.message = error.what();
            return result;
        }

        if (!map.isValid(
                startPose.position.x,
                startPose.position.y,
                startPose.position.z)) {
            result.error_code = ErrorCode::START_OUT_OF_BOUNDS;
            result.message = "Start pose is outside the voxel map.";
            return result;
        }
        if (!map.isValid(
                endPose.position.x,
                endPose.position.y,
                endPose.position.z)) {
            result.error_code = ErrorCode::END_OUT_OF_BOUNDS;
            result.message = "End pose is outside the voxel map.";
            return result;
        }
        if (startPose.poseId >= poses.size() ||
            startPose.poseId >= map.poseCount()) {
            result.error_code = ErrorCode::START_POSE_INVALID;
            result.message = "Start poseId is outside the generated pose table.";
            return result;
        }
        if (endPose.poseId >= poses.size() ||
            endPose.poseId >= map.poseCount()) {
            result.error_code = ErrorCode::END_POSE_INVALID;
            result.message = "End poseId is outside the generated pose table.";
            return result;
        }

        const std::size_t startIndex = map.index(
            static_cast<std::uint32_t>(startPose.position.x),
            static_cast<std::uint32_t>(startPose.position.y),
            static_cast<std::uint32_t>(startPose.position.z));
        const std::size_t goalIndex = map.index(
            static_cast<std::uint32_t>(endPose.position.x),
            static_cast<std::uint32_t>(endPose.position.y),
            static_cast<std::uint32_t>(endPose.position.z));
        startPose.poseId = resolveEndpointPose(
            map,
            startIndex,
            startPose.poseId,
            poses);
        endPose.poseId = resolveEndpointPose(
            map,
            goalIndex,
            endPose.poseId,
            poses);
        if (startPose.poseId == kInvalidPoseId) {
            result.error_code = ErrorCode::START_POINT_BLOCKED;
            result.message = "No valid pose is available at the start point.";
            return result;
        }
        if (endPose.poseId == kInvalidPoseId) {
            result.error_code = ErrorCode::END_POINT_BLOCKED;
            result.message = "No valid pose is available at the goal point.";
            return result;
        }

        if (startPose.position == endPose.position) {
            PathResult path;
            path.path.push_back(startPose.position);
            appendConditionalPose(
                map,
                startPose.position,
                poses[startPose.poseId],
                0U,
                path.pose_description);
            path.cost = 0.0;
            result.paths.push_back(std::move(path));
            result.status = PlannerStatus::OK;
            result.error_code = ErrorCode::NONE;
            result.message = "Returned the zero-distance path.";
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
            result.message = "Memory allocation failed while building "
                             "straight sweep templates.";
            return result;
        } catch (const std::exception& error) {
            result.error_code = ErrorCode::INVALID_ARGUMENT;
            result.message = error.what();
            return result;
        }

        const std::uint64_t startKey = stateKey(
            startIndex,
            startPose.poseId,
            poses.size());
        const std::uint64_t goalKey = stateKey(
            goalIndex,
            endPose.poseId,
            poses.size());

        const float topologyHintTolerance =
            std::max(config.busbar_width * 2.0F, 30.0F);
        HintDistanceLut hintLut;
        BackwardDistanceField backwardDistance;
        try {
            hintLut = buildHintDistanceLut(
                map,
                topologyHint,
                topologyHintTolerance);
            backwardDistance = buildBackwardDistanceField(
                map,
                startIndex,
                goalIndex,
                rawStart,
                rawGoal,
                endpointImmunityRadius,
                hintLut,
                allowFullSearch);
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
        bool goalToleranceAccepted = false;
        float bestGoalToleranceDistance =
            std::numeric_limits<float>::infinity();

        while (static_cast<int>(result.paths.size()) < boundedMaxPaths) {
            SearchOutcome outcome;
            try {
                outcome = runSingleSearch(
                    map,
                    startKey,
                    goalKey,
                    poses,
                    actions,
                    poseLut,
                    backwardDistance,
                    penalties,
                    rawStart,
                    rawGoal,
                    endpointImmunityRadius,
                    maxExpansionsOverride);
            } catch (const std::bad_alloc&) {
                result.error_code = ErrorCode::MEMORY_ALLOCATION_FAILED;
                result.message =
                    "Memory allocation failed during lightweight search.";
                return result;
            }

            if (outcome.path.path.empty()) {
                if (outcome.error_code ==
                    ErrorCode::COMPUTATION_LIMIT_EXCEEDED) {
                    result.error_code = result.paths.empty()
                        ? ErrorCode::COMPUTATION_LIMIT_EXCEEDED
                        : ErrorCode::NONE;
                    result.status = result.paths.empty()
                        ? PlannerStatus::NO_PATH
                        : PlannerStatus::OK;
                    result.message =
                        "26-neighborhood expansion limit exceeded after " +
                        std::to_string(outcome.expansions) +
                        " expansions; closest goal distance=" +
                        std::to_string(outcome.closestGoalDistance) +
                        ", generated states=" +
                        std::to_string(outcome.generatedStates) +
                        "; returning paths completed so far.";
                    return result;
                }
                break;
            }

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
            goalToleranceAccepted =
                goalToleranceAccepted || outcome.goalToleranceAccepted;
            bestGoalToleranceDistance = std::min(
                bestGoalToleranceDistance,
                outcome.closestGoalDistance);
            result.paths.push_back(std::move(outcome.path));
            if (static_cast<int>(result.paths.size()) < boundedMaxPaths) {
                addSoftPathPenalty(
                    map,
                    result.paths.back().path,
                    rawStart,
                    rawGoal,
                    endpointImmunityRadius,
                    penalties);
            }
        }

        if (result.paths.empty()) {
            result.status = PlannerStatus::NO_PATH;
            result.error_code = ErrorCode::PATH_NOT_FOUND;
            result.message = "No feasible 26-neighborhood path exists.";
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
            " lightweight path(s).";
        if (rejectedOverlapCandidates != 0) {
            result.message += " Rejected " +
                std::to_string(rejectedOverlapCandidates) +
                " overlapping candidate(s).";
        }
        if (diversityRetryLimitReached) {
            result.message += " Diversity retry limit reached.";
        }
        if (goalToleranceAccepted) {
            result.message += " Goal tolerance accepted at distance " +
                std::to_string(bestGoalToleranceDistance) + ".";
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
        const std::uint32_t startPoseId = resolveEndpointPose(
            map,
            startIndex,
            0U,
            poses);
        const std::uint32_t goalPoseId = resolveEndpointPose(
            map,
            goalIndex,
            startPoseId,
            poses);
        if (startPoseId == kInvalidPoseId) {
            result.error_code = ErrorCode::START_POINT_BLOCKED;
            result.message = "No valid pose is available at the start point.";
            return result;
        }
        if (goalPoseId == kInvalidPoseId) {
            result.error_code = ErrorCode::END_POINT_BLOCKED;
            result.message = "No valid pose is available at the goal point.";
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

    static PlanningResult findPaths(
        const VoxelGrid&,
        Point3D,
        Point3D,
        int) {
        PlanningResult result;
        result.error_code = ErrorCode::POSE_MASK_UNAVAILABLE;
        result.message =
            "The lightweight planner requires populated pose masks and a "
            "pose table.";
        return result;
    }

private:
    struct Vector {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
    };

    struct Delta {
        int dx = 0;
        int dy = 0;
        int dz = 0;
    };

    static constexpr std::size_t kDirectionCount = 26U;
    static constexpr std::size_t kMaxDirectionalCandidates = 8U;
    static constexpr std::size_t kMaxTransitionCandidates = 128U;

    struct PoseCandidateList {
        std::array<std::uint32_t, kMaxDirectionalCandidates> ids{};
        std::size_t size = 0U;
    };

    struct PoseTransitionList {
        std::array<std::uint32_t, kMaxTransitionCandidates> ids{};
        std::size_t size = 0U;
    };

    /**
     * @brief Immutable pose pruning tables built before the search starts.
     *
     * dir_to_poses contains the generated roll variants for each of the
     * 26 exact tangent directions. directional_transition is the direct
     * intersection of that table with pose_transition, so the hot path does
     * not allocate or scan the complete pose pool.
     */
    struct PoseLut {
        std::array<std::vector<std::uint32_t>, kDirectionCount>
            dir_to_poses;
        std::vector<PoseTransitionList> pose_transition;
        std::vector<std::array<PoseCandidateList, kDirectionCount>>
            directional_transition;
    };

    struct ActionCatalog {
        std::vector<std::vector<std::vector<voxel_planner::VoxelOffset>>>
            straightSweeps;
        std::vector<std::vector<voxel_planner::VoxelOffsetBounds>>
            straightSweepBounds;
    };

    static constexpr std::uint32_t kInvalidPoseId =
        std::numeric_limits<std::uint32_t>::max();
    static constexpr std::uint64_t kNoParent =
        std::numeric_limits<std::uint64_t>::max();
    static constexpr std::uint8_t kOpen = 1U;
    static constexpr std::uint8_t kClosed = 2U;
    static constexpr int kMaxPathsPerRequest = 10;
    static constexpr int kMaxPenaltyRetriesPerPath = 10;
    static constexpr std::size_t kMaxExpansionsPerPath = 50000U;
    static constexpr float kGoalToleranceDistance = 5.0F;
    static constexpr float kHeuristicWeight = 1.5F;
    static constexpr double kGeometryTolerance = 1e-4;
    static constexpr double kPi = 3.14159265358979323846;
    static constexpr float kCenterlineReusePenalty = 100000.0F;
    static constexpr float kNeighborReusePenalty = 1000.0F;
    static constexpr double kMinTangentContinuity = 0.5;
    static constexpr double kMinNormalContinuity = 0.95;
    static constexpr int kMaxRollStepDistance = 2;

    struct Node {
        std::uint64_t stateKey = 0U;
        std::size_t voxelIndex = 0U;
        std::uint32_t poseId = 0U;
        float g = std::numeric_limits<float>::infinity();
        float f = std::numeric_limits<float>::infinity();
        float h = std::numeric_limits<float>::infinity();
    };

    struct StateRecord {
        float g = std::numeric_limits<float>::infinity();
        std::uint64_t parent = kNoParent;
        std::uint8_t status = 0U;
    };

    struct NodeCompare {
        bool operator()(const Node& lhs, const Node& rhs) const noexcept {
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

    struct SearchOutcome {
        PathResult path;
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
        std::uint64_t closestGoalStateKey = kNoParent;
        float closestGoalCost = std::numeric_limits<float>::infinity();
        bool goalToleranceAccepted = false;
    };

    struct BackwardDistanceField {
        static constexpr std::uint32_t kUnreachable =
            std::numeric_limits<std::uint32_t>::max();
        std::uint32_t width = 0U;
        std::uint32_t height = 0U;
        std::uint32_t depth = 0U;
        std::vector<std::uint32_t> distances;
        std::size_t reachableVoxels = 0U;

        std::uint32_t distance(std::size_t index) const noexcept {
            return index < distances.size()
                ? distances[index]
                : kUnreachable;
        }
    };

    /**
     * @brief Per-request O(1) membership map for the topology hint band.
     *
     * A value of 0 denotes a voxel on the hint centerline. A value of 1
     * denotes a voxel inside the configured hint tolerance. -1 means that
     * the voxel is outside the hint band.
     */
    struct HintDistanceLut {
        std::vector<std::int32_t> hint_distance_map;

        bool enabled() const noexcept {
            return !hint_distance_map.empty();
        }

        bool contains(std::size_t voxelIndex) const noexcept {
            return voxelIndex < hint_distance_map.size() &&
                hint_distance_map[voxelIndex] >= 0;
        }
    };

    struct CostPenaltyMap {
        std::unordered_map<std::size_t, float> values;

        float get(std::size_t index) const noexcept {
            const auto entry = values.find(index);
            return entry == values.end() ? 0.0F : entry->second;
        }
    };

    static const std::array<Delta, 26>& neighborDeltas() {
        static const std::array<Delta, 26> deltas = [] {
            std::array<Delta, 26> values{};
            std::size_t index = 0U;
            for (int dz = -1; dz <= 1; ++dz) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        if (dx == 0 && dy == 0 && dz == 0) {
                            continue;
                        }
                        values[index++] = {dx, dy, dz};
                    }
                }
            }
            return values;
        }();
        return deltas;
    }

    static std::size_t neighborIndex(const Delta& delta) noexcept {
        std::size_t index = 0U;
        for (const Delta& candidate : neighborDeltas()) {
            if (candidate.dx == delta.dx &&
                candidate.dy == delta.dy &&
                candidate.dz == delta.dz) {
                return index;
            }
            ++index;
        }
        return 0U;
    }

    static Vector unitTangent(
        const voxel_planner::BusbarPose& pose) {
        const double length = std::sqrt(
            static_cast<double>(pose.tx) * pose.tx +
            static_cast<double>(pose.ty) * pose.ty +
            static_cast<double>(pose.tz) * pose.tz);
        if (length <= kGeometryTolerance) {
            throw std::invalid_argument("Pose tangent must not be zero.");
        }
        return {
            pose.tx / length,
            pose.ty / length,
            pose.tz / length};
    }

    static Vector unitNormal(
        const voxel_planner::BusbarPose& pose) {
        const double length = std::sqrt(
            static_cast<double>(pose.nx) * pose.nx +
            static_cast<double>(pose.ny) * pose.ny +
            static_cast<double>(pose.nz) * pose.nz);
        if (length <= kGeometryTolerance) {
            throw std::invalid_argument("Pose normal must not be zero.");
        }
        return {
            pose.nx / length,
            pose.ny / length,
            pose.nz / length};
    }

    static double dot(const Vector& lhs, const Vector& rhs) noexcept {
        return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
    }

    static bool poseTransitionIsAllowed(
        const voxel_planner::BusbarPose& currentPose,
        const voxel_planner::BusbarPose& candidatePose,
        const Vector& currentTangent,
        const Vector& currentNormal,
        const Vector& candidateTangent,
        const Vector& candidateNormal,
        int rollStepDegrees) noexcept {
        const int rollDifference = std::abs(
            currentPose.roll_angle_deg - candidatePose.roll_angle_deg);
        const int circularRollDifference = std::min(
            rollDifference,
            360 - rollDifference);
        return dot(currentTangent, candidateTangent) >=
                kMinTangentContinuity &&
            dot(currentNormal, candidateNormal) >=
                kMinNormalContinuity &&
            circularRollDifference <=
                kMaxRollStepDistance * rollStepDegrees;
    }

    static void appendPoseTransition(
        PoseTransitionList& list,
        std::uint32_t poseId) {
        if (list.size >= list.ids.size()) {
            throw std::overflow_error(
                "Pose transition LUT candidate capacity exceeded.");
        }
        list.ids[list.size++] = poseId;
    }

    static void appendDirectionalCandidate(
        PoseCandidateList& list,
        std::uint32_t poseId) {
        if (list.size >= list.ids.size()) {
            throw std::overflow_error(
                "Directional pose LUT candidate capacity exceeded.");
        }
        list.ids[list.size++] = poseId;
    }

    static PoseLut buildPoseLut(
        const std::vector<voxel_planner::BusbarPose>& poses,
        int rollStepDegrees) {
        const int normalizedRollStepDegrees = std::max(1, rollStepDegrees);
        PoseLut lut;
        lut.pose_transition.resize(poses.size());
        lut.directional_transition.resize(poses.size());

        std::vector<Vector> tangents(poses.size());
        std::vector<Vector> normals(poses.size());
        std::vector<std::size_t> poseDirections(poses.size());
        for (const voxel_planner::BusbarPose& pose : poses) {
            if (pose.poseId >= poses.size()) {
                throw std::invalid_argument(
                    "Pose identifiers must be contiguous for LUT creation.");
            }
            const std::size_t direction = neighborIndex(
                Delta{pose.tx, pose.ty, pose.tz});
            if (neighborDeltas()[direction].dx != pose.tx ||
                neighborDeltas()[direction].dy != pose.ty ||
                neighborDeltas()[direction].dz != pose.tz) {
                throw std::invalid_argument(
                    "Pose tangent is not one of the 26 voxel directions.");
            }
            lut.dir_to_poses[direction].push_back(pose.poseId);
            poseDirections[pose.poseId] = direction;
            tangents[pose.poseId] = unitTangent(pose);
            normals[pose.poseId] = unitNormal(pose);
        }

        for (std::size_t currentPoseId = 0U;
             currentPoseId < poses.size();
             ++currentPoseId) {
            PoseTransitionList& transitions =
                lut.pose_transition[currentPoseId];
            for (std::size_t candidatePoseId = 0U;
                 candidatePoseId < poses.size();
                 ++candidatePoseId) {
                if (poseTransitionIsAllowed(
                        poses[currentPoseId],
                        poses[candidatePoseId],
                        tangents[currentPoseId],
                        normals[currentPoseId],
                        tangents[candidatePoseId],
                        normals[candidatePoseId],
                        normalizedRollStepDegrees)) {
                    appendPoseTransition(
                        transitions,
                        static_cast<std::uint32_t>(candidatePoseId));
                }
            }

            for (std::size_t direction = 0U;
                 direction < kDirectionCount;
                 ++direction) {
                PoseCandidateList& candidates =
                    lut.directional_transition[
                        currentPoseId][direction];
                for (std::size_t transitionIndex = 0U;
                     transitionIndex < transitions.size;
                     ++transitionIndex) {
                    const std::uint32_t candidatePoseId =
                        transitions.ids[transitionIndex];
                    if (poseDirections[candidatePoseId] == direction) {
                        appendDirectionalCandidate(
                            candidates,
                            candidatePoseId);
                    }
                }
            }
        }
        return lut;
    }

    static std::uint64_t stateKey(
        std::size_t voxelIndex,
        std::uint32_t poseId,
        std::size_t poseCount) {
        if (poseCount == 0U ||
            voxelIndex >
                (std::numeric_limits<std::uint64_t>::max() - poseId) /
                    poseCount) {
            throw std::overflow_error("Lightweight state key overflow.");
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

    static HintDistanceLut buildHintDistanceLut(
        const VoxelGrid& map,
        const std::vector<Point3D>* topologyHint,
        float topologyHintTolerance) {
        HintDistanceLut lut;
        if (topologyHint == nullptr ||
            topologyHint->empty() ||
            map.voxelCount() == 0U ||
            !std::isfinite(topologyHintTolerance) ||
            topologyHintTolerance < 0.0F) {
            return lut;
        }

        lut.hint_distance_map.assign(
            map.voxelCount(),
            -1);
        const double radius =
            static_cast<double>(topologyHintTolerance);
        const double radiusSquared = radius * radius;
        const int radiusVoxels = static_cast<int>(std::ceil(radius));
        for (const Point3D& hint : *topologyHint) {
            if (!map.isValid(hint.x, hint.y, hint.z)) {
                continue;
            }
            for (int dz = -radiusVoxels;
                 dz <= radiusVoxels;
                 ++dz) {
                for (int dy = -radiusVoxels;
                     dy <= radiusVoxels;
                     ++dy) {
                    for (int dx = -radiusVoxels;
                         dx <= radiusVoxels;
                         ++dx) {
                        const double distanceSquared =
                            static_cast<double>(dx) * dx +
                            static_cast<double>(dy) * dy +
                            static_cast<double>(dz) * dz;
                        if (distanceSquared > radiusSquared) {
                            continue;
                        }
                        const int x = hint.x + dx;
                        const int y = hint.y + dy;
                        const int z = hint.z + dz;
                        if (!map.isValid(x, y, z)) {
                            continue;
                        }
                        const std::size_t voxelIndex = map.index(
                            static_cast<std::uint32_t>(x),
                            static_cast<std::uint32_t>(y),
                            static_cast<std::uint32_t>(z));
                        if (dx == 0 && dy == 0 && dz == 0) {
                            lut.hint_distance_map[voxelIndex] = 0;
                        } else if (lut.hint_distance_map[voxelIndex] < 0) {
                            lut.hint_distance_map[voxelIndex] = 1;
                        }
                    }
                }
            }
        }
        return lut;
    }

    static bool poseTableMatchesMap(
        const VoxelGrid& map,
        const std::vector<voxel_planner::BusbarPose>& poses) {
        if (poses.empty() || map.poseCount() != poses.size()) {
            return false;
        }
        std::vector<std::uint8_t> seen(poses.size(), 0U);
        for (const voxel_planner::BusbarPose& pose : poses) {
            if (pose.poseId >= poses.size() ||
                seen[pose.poseId] != 0U) {
                return false;
            }
            seen[pose.poseId] = 1U;
        }
        return true;
    }

    static std::uint32_t firstAllowedPose(
        const VoxelGrid& map,
        std::size_t voxelIndex) {
        if (voxelIndex >= map.voxelCount() ||
            map.isObstacle(voxelIndex)) {
            return kInvalidPoseId;
        }
        for (std::uint32_t poseId = 0U;
             poseId < map.poseCount();
             ++poseId) {
            if (map.isPoseAllowedForSearch(voxelIndex, poseId)) {
                return poseId;
            }
        }
        return kInvalidPoseId;
    }

    static std::uint32_t resolveEndpointPose(
        const VoxelGrid& map,
        std::size_t voxelIndex,
        std::uint32_t preferredPoseId,
        const std::vector<voxel_planner::BusbarPose>& poses) {
        if (voxelIndex >= map.voxelCount() ||
            map.isObstacle(voxelIndex)) {
            return kInvalidPoseId;
        }
        if (map.getState(voxelIndex) == VoxelState::UNCONDITIONAL) {
            return preferredPoseId < poses.size()
                ? preferredPoseId
                : 0U;
        }
        if (preferredPoseId < poses.size() &&
            map.isPoseAllowedForSearch(voxelIndex, preferredPoseId)) {
            return preferredPoseId;
        }
        const std::uint32_t fallback = firstAllowedPose(map, voxelIndex);
        if (fallback == kInvalidPoseId ||
            preferredPoseId >= poses.size()) {
            return fallback;
        }
        const Vector preferredTangent = unitTangent(poses[preferredPoseId]);
        const Vector preferredNormal = unitNormal(poses[preferredPoseId]);
        double bestScore = -std::numeric_limits<double>::infinity();
        std::uint32_t selected = fallback;
        for (const voxel_planner::BusbarPose& pose : poses) {
            if (!map.isPoseAllowedForSearch(voxelIndex, pose.poseId)) {
                continue;
            }
            const double score =
                dot(preferredTangent, unitTangent(pose)) +
                dot(preferredNormal, unitNormal(pose));
            if (score > bestScore) {
                bestScore = score;
                selected = pose.poseId;
            }
        }
        return selected;
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

    static std::vector<voxel_planner::VoxelOffset> makeUnitSweep(
        const std::vector<voxel_planner::VoxelOffset>& footprint,
        const Delta& delta) {
        std::unordered_set<std::uint64_t> seen;
        std::vector<voxel_planner::VoxelOffset> sweep;
        sweep.reserve(footprint.size() * 3U);
        for (int sample = 0; sample <= 2; ++sample) {
            const int tx = static_cast<int>(std::lround(
                static_cast<double>(delta.dx) * sample / 2.0));
            const int ty = static_cast<int>(std::lround(
                static_cast<double>(delta.dy) * sample / 2.0));
            const int tz = static_cast<int>(std::lround(
                static_cast<double>(delta.dz) * sample / 2.0));
            for (const voxel_planner::VoxelOffset& offset : footprint) {
                const int x = offset.dx + tx;
                const int y = offset.dy + ty;
                const int z = offset.dz + tz;
                const std::uint64_t key =
                    (static_cast<std::uint64_t>(
                        static_cast<std::uint32_t>(x)) << 42U) ^
                    (static_cast<std::uint64_t>(
                        static_cast<std::uint32_t>(y)) << 21U) ^
                    static_cast<std::uint64_t>(
                        static_cast<std::uint32_t>(z));
                if (seen.insert(key).second) {
                    sweep.push_back({x, y, z});
                }
            }
        }
        return sweep;
    }

    static ActionCatalog buildActionCatalog(
        const std::vector<voxel_planner::BusbarPose>& poses,
        const voxel_planner::PlannerConfig& config) {
        ActionCatalog actions;
        actions.straightSweeps.resize(poses.size());
        actions.straightSweepBounds.resize(poses.size());
        for (const voxel_planner::BusbarPose& pose : poses) {
            const std::vector<voxel_planner::VoxelOffset> footprint =
                voxel_planner::generatePoseFootprint(pose, config);
            actions.straightSweeps[pose.poseId].resize(26U);
            actions.straightSweepBounds[pose.poseId].resize(26U);
            std::size_t direction = 0U;
            for (const Delta& delta : neighborDeltas()) {
                actions.straightSweeps[pose.poseId][direction] =
                    makeUnitSweep(footprint, delta);
                actions.straightSweepBounds[pose.poseId][direction] =
                    offsetBounds(
                        actions.straightSweeps[pose.poseId][direction]);
                ++direction;
            }
        }
        return actions;
    }

    static bool isSweepCollisionFreeForSearch(
        const VoxelGrid& map,
        const Point3D& anchor,
        const std::vector<voxel_planner::VoxelOffset>& sweep,
        const voxel_planner::VoxelOffsetBounds& bounds) noexcept {
        if (sweep.empty() || !bounds.valid) {
            return false;
        }
        if (map.morphologyPrefix_ &&
            module2_morphology::VoxelMorphologyEngine::isBoxCollisionFree(
                *map.morphologyPrefix_,
                map.prefixWidth_,
                map.prefixHeight_,
                map.prefixDepth_,
                anchor.x,
                anchor.y,
                anchor.z,
                bounds)) {
            return true;
        }
        for (const voxel_planner::VoxelOffset& offset : sweep) {
            const std::int64_t x =
                static_cast<std::int64_t>(anchor.x) + offset.dx;
            const std::int64_t y =
                static_cast<std::int64_t>(anchor.y) + offset.dy;
            const std::int64_t z =
                static_cast<std::int64_t>(anchor.z) + offset.dz;
            if (x < 0 || y < 0 || z < 0 ||
                x >= static_cast<std::int64_t>(map.width()) ||
                y >= static_cast<std::int64_t>(map.height()) ||
                z >= static_cast<std::int64_t>(map.depth())) {
                return false;
            }
            if (map.isObstacle(map.index(
                    static_cast<std::uint32_t>(x),
                    static_cast<std::uint32_t>(y),
                    static_cast<std::uint32_t>(z)))) {
                return false;
            }
        }
        return true;
    }

    static BackwardDistanceField buildBackwardDistanceField(
        const VoxelGrid& map,
        std::size_t startVoxel,
        std::size_t goalVoxel,
        const Point3D& rawStart,
        const Point3D& rawGoal,
        float endpointImmunityRadius,
        const HintDistanceLut& hintLut,
        bool allowFullSearch) {
        BackwardDistanceField field;
        field.width = map.width();
        field.height = map.height();
        field.depth = map.depth();
        field.distances.assign(
            map.voxelCount(),
            BackwardDistanceField::kUnreachable);
        if (goalVoxel >= map.voxelCount()) {
            return field;
        }
        const std::size_t plane =
            static_cast<std::size_t>(field.width) * field.height;
        std::vector<std::size_t> frontier;
        frontier.push_back(goalVoxel);
        field.distances[goalVoxel] = 0U;
        field.reachableVoxels = 1U;
        while (!frontier.empty()) {
            std::vector<std::size_t> next;
            next.reserve(frontier.size() * 2U);
            for (const std::size_t current : frontier) {
                const std::size_t planeIndex = current % plane;
                const std::uint32_t z = static_cast<std::uint32_t>(
                    current / plane);
                const std::uint32_t y = static_cast<std::uint32_t>(
                    planeIndex / field.width);
                const std::uint32_t x = static_cast<std::uint32_t>(
                    planeIndex % field.width);
                const std::array<std::size_t, 6> neighbors = {
                    x > 0U ? current - 1U : current,
                    x + 1U < field.width ? current + 1U : current,
                    y > 0U ? current - field.width : current,
                    y + 1U < field.height ? current + field.width : current,
                    z > 0U ? current - plane : current,
                    z + 1U < field.depth ? current + plane : current};
                for (const std::size_t candidate : neighbors) {
                    if (candidate == current ||
                        field.distances[candidate] !=
                            BackwardDistanceField::kUnreachable) {
                        continue;
                    }
                    if (!allowFullSearch &&
                        hintLut.enabled() &&
                        candidate != startVoxel &&
                        candidate != goalVoxel) {
                        if (!hintLut.contains(candidate)) {
                            continue;
                        }
                    }
                    if (map.isObstacle(candidate) &&
                        candidate != startVoxel &&
                        candidate != goalVoxel &&
                        !isNearEndpoint(
                            pointFromVoxelIndex(map, candidate),
                            rawStart,
                            rawGoal,
                            endpointImmunityRadius)) {
                        continue;
                    }
                    field.distances[candidate] =
                        field.distances[current] + 1U;
                    next.push_back(candidate);
                    ++field.reachableVoxels;
                }
            }
            frontier.swap(next);
        }
        if (!allowFullSearch &&
            hintLut.enabled() &&
            field.distance(startVoxel) ==
                BackwardDistanceField::kUnreachable) {
            return buildBackwardDistanceField(
                map,
                startVoxel,
                goalVoxel,
                rawStart,
                rawGoal,
                endpointImmunityRadius,
                hintLut,
                true);
        }
        return field;
    }

    static bool isNearEndpoint(
        const Point3D& point,
        const Point3D& rawStart,
        const Point3D& rawGoal,
        float radius) {
        const double limit = static_cast<double>(radius) * radius;
        const auto distanceSquared = [&point](const Point3D& endpoint) {
            const double dx = static_cast<double>(point.x) - endpoint.x;
            const double dy = static_cast<double>(point.y) - endpoint.y;
            const double dz = static_cast<double>(point.z) - endpoint.z;
            return dx * dx + dy * dy + dz * dz;
        };
        return radius > 0.0F &&
            (distanceSquared(rawStart) < limit ||
             distanceSquared(rawGoal) < limit);
    }

    static float heuristic(
        std::size_t voxelIndex,
        const Point3D& point,
        const Point3D& goal,
        const Vector& currentTangent,
        const Vector& goalTangent,
        const BackwardDistanceField& backwardDistance) {
        const double dx = static_cast<double>(goal.x) - point.x;
        const double dy = static_cast<double>(goal.y) - point.y;
        const double dz = static_cast<double>(goal.z) - point.z;
        const float euclidean = static_cast<float>(std::sqrt(
            dx * dx + dy * dy + dz * dz));
        const std::uint32_t guided = backwardDistance.distance(voxelIndex);
        const float topological = guided ==
                BackwardDistanceField::kUnreachable
            ? euclidean
            : static_cast<float>(guided);
        const float alignment = static_cast<float>(std::max(
            -1.0,
            std::min(1.0, dot(currentTangent, goalTangent))));
        const float orientationPenalty =
            (1.0F - alignment) * 5.0F;
        const std::uint64_t mixed =
            static_cast<std::uint64_t>(voxelIndex) *
            std::uint64_t{11400714819323198485ULL};
        const float tieBreaker = static_cast<float>(
            mixed % 1000U) * 1.0e-6F;
        return std::max(euclidean, topological) +
            orientationPenalty +
            tieBreaker;
    }

    static std::uint32_t resolveTargetPose(
        const VoxelGrid& map,
        std::size_t targetVoxel,
        std::uint32_t currentPoseId,
        const std::vector<voxel_planner::BusbarPose>& poses,
        const Delta& movement,
        const PoseLut& poseLut) {
        if (targetVoxel >= map.voxelCount() ||
            map.isObstacle(targetVoxel) ||
            currentPoseId >= poses.size()) {
            return kInvalidPoseId;
        }
        if (map.getState(targetVoxel) == VoxelState::UNCONDITIONAL ||
            map.isPoseAllowedForSearch(targetVoxel, currentPoseId)) {
            return currentPoseId;
        }

        const std::size_t direction = neighborIndex(movement);
        const PoseCandidateList& candidates =
            poseLut.directional_transition[currentPoseId][direction];
        const Vector currentNormal = unitNormal(poses[currentPoseId]);
        const Vector currentTangent = unitTangent(poses[currentPoseId]);
        double bestScore = -std::numeric_limits<double>::infinity();
        std::uint32_t selected = kInvalidPoseId;
        for (std::size_t index = 0U;
             index < candidates.size;
             ++index) {
            const std::uint32_t poseId = candidates.ids[index];
            if (poseId >= poses.size() ||
                !map.isPoseAllowedForSearch(targetVoxel, poseId)) {
                continue;
            }
            const double score =
                dot(
                    currentNormal,
                    unitNormal(poses[poseId])) +
                dot(
                    currentTangent,
                    unitTangent(poses[poseId]));
            if (score > bestScore) {
                bestScore = score;
                selected = poseId;
            }
        }
        return selected;
    }

    static SearchOutcome runSingleSearch(
        const VoxelGrid& map,
        std::uint64_t startKey,
        std::uint64_t goalKey,
        const std::vector<voxel_planner::BusbarPose>& poses,
        const ActionCatalog& actions,
        const PoseLut& poseLut,
        const BackwardDistanceField& backwardDistance,
        const CostPenaltyMap& penalties,
        const Point3D& rawStart,
        const Point3D& rawGoal,
        float endpointImmunityRadius,
        std::size_t maxExpansionsOverride) {
        SearchOutcome outcome;
        std::unordered_map<std::uint64_t, StateRecord> records;
        records.reserve(131072U);
        std::priority_queue<Node, std::vector<Node>, NodeCompare> openList;

        const std::size_t startVoxel = startKey / poses.size();
        const std::uint32_t startPoseId = static_cast<std::uint32_t>(
            startKey % poses.size());
        const std::size_t goalVoxel = goalKey / poses.size();
        const std::uint32_t goalPoseId = static_cast<std::uint32_t>(
            goalKey % poses.size());
        const Point3D goalPoint = pointFromVoxelIndex(map, goalVoxel);
        const Vector goalTangent = unitTangent(poses[goalPoseId]);
        const std::size_t maxExpansions = maxExpansionsOverride == 0U
            ? kMaxExpansionsPerPath
            : std::max(kMaxExpansionsPerPath, maxExpansionsOverride);

        const Point3D startPoint = pointFromVoxelIndex(map, startVoxel);
        const float startH = heuristic(
            startVoxel,
            startPoint,
            goalPoint,
            unitTangent(poses[startPoseId]),
            goalTangent,
            backwardDistance);
        records[startKey] = {0.0F, kNoParent, kOpen};
        openList.push({
            startKey,
            startVoxel,
            startPoseId,
            0.0F,
            kHeuristicWeight * startH,
            startH});
        outcome.backwardReachableVoxels = backwardDistance.reachableVoxels;
        outcome.backwardStartDistance = backwardDistance.distance(startVoxel);
        outcome.closestBackwardDistance = outcome.backwardStartDistance;

        std::size_t expansions = 0U;
        while (!openList.empty() && expansions < maxExpansions) {
            const Node current = openList.top();
            openList.pop();
            auto record = records.find(current.stateKey);
            if (record == records.end() ||
                record->second.status != kOpen ||
                current.g > record->second.g + 1e-5F) {
                continue;
            }
            record->second.status = kClosed;
            ++expansions;

            const Point3D currentPoint = pointFromVoxelIndex(
                map,
                current.voxelIndex);
            const float currentGoalDistance = static_cast<float>(std::sqrt(
                std::pow(
                    static_cast<double>(goalPoint.x - currentPoint.x),
                    2.0) +
                std::pow(
                    static_cast<double>(goalPoint.y - currentPoint.y),
                    2.0) +
                std::pow(
                    static_cast<double>(goalPoint.z - currentPoint.z),
                    2.0)));
            if (currentGoalDistance < outcome.closestGoalDistance) {
                outcome.closestGoalDistance = currentGoalDistance;
                outcome.closestGoalStateKey = current.stateKey;
                outcome.closestGoalCost = current.g;
            }
            outcome.closestBackwardDistance = std::min(
                outcome.closestBackwardDistance,
                backwardDistance.distance(current.voxelIndex));

            if (current.voxelIndex == goalVoxel) {
                reconstructPath(
                    map,
                    records,
                    startKey,
                    current.stateKey,
                    poses,
                    outcome);
                outcome.path.cost = current.g;
                outcome.error_code = ErrorCode::NONE;
                outcome.expansions = expansions;
                outcome.generatedStates = records.size();
                outcome.remainingOpenStates = openList.size();
                return outcome;
            }

            for (const Delta& delta : neighborDeltas()) {
                const Point3D targetPoint{
                    currentPoint.x + delta.dx,
                    currentPoint.y + delta.dy,
                    currentPoint.z + delta.dz};
                if (!map.isValid(
                        targetPoint.x,
                        targetPoint.y,
                        targetPoint.z)) {
                    continue;
                }

                const std::size_t targetVoxel = map.index(
                    static_cast<std::uint32_t>(targetPoint.x),
                    static_cast<std::uint32_t>(targetPoint.y),
                    static_cast<std::uint32_t>(targetPoint.z));
                const std::size_t direction = neighborIndex(delta);
                const std::uint32_t targetPoseId = resolveTargetPose(
                    map,
                    targetVoxel,
                    current.poseId,
                    poses,
                    delta,
                    poseLut);
                if (targetPoseId == kInvalidPoseId) {
                    continue;
                }
                const bool endpointImmune = isNearEndpoint(
                    targetPoint,
                    rawStart,
                    rawGoal,
                    endpointImmunityRadius);
                if (!endpointImmune &&
                    !isSweepCollisionFreeForSearch(
                        map,
                        currentPoint,
                        actions.straightSweeps[targetPoseId][direction],
                        actions.straightSweepBounds[targetPoseId][direction])) {
                    continue;
                }

                const std::uint64_t targetKey = stateKey(
                    targetVoxel,
                    targetPoseId,
                    poses.size());
                const float stepCost = static_cast<float>(std::sqrt(
                    static_cast<double>(
                        delta.dx * delta.dx +
                        delta.dy * delta.dy +
                        delta.dz * delta.dz)));
                const float tentativeG = current.g +
                    stepCost +
                    (targetVoxel == current.voxelIndex
                        ? 0.0F
                        : penalties.get(targetVoxel));
                const auto existing = records.find(targetKey);
                if (existing != records.end() &&
                    (existing->second.status == kClosed ||
                     tentativeG + 1e-5F >= existing->second.g)) {
                    continue;
                }
                records[targetKey] = {
                    tentativeG,
                    current.stateKey,
                    kOpen};
                const float h = heuristic(
                    targetVoxel,
                    targetPoint,
                    goalPoint,
                    unitTangent(poses[goalPoseId]),
                    unitTangent(poses[targetPoseId]),
                    backwardDistance);
                openList.push({
                    targetKey,
                    targetVoxel,
                    targetPoseId,
                    tentativeG,
                    tentativeG + kHeuristicWeight * h,
                    h});
            }
        }

        outcome.expansions = expansions;
        outcome.generatedStates = records.size();
        outcome.remainingOpenStates = openList.size();
        if (expansions >= maxExpansions) {
            if (outcome.closestGoalStateKey != kNoParent &&
                outcome.closestGoalDistance <= kGoalToleranceDistance) {
                reconstructPath(
                    map,
                    records,
                    startKey,
                    outcome.closestGoalStateKey,
                    poses,
                    outcome);
                if (!outcome.path.path.empty()) {
                    outcome.path.cost = outcome.closestGoalCost;
                    outcome.path.goal_tolerance_accepted = true;
                    outcome.goalToleranceAccepted = true;
                    outcome.error_code = ErrorCode::NONE;
                    return outcome;
                }
            }
            outcome.error_code = ErrorCode::COMPUTATION_LIMIT_EXCEEDED;
        }
        return outcome;
    }

    static void reconstructPath(
        const VoxelGrid& map,
        const std::unordered_map<std::uint64_t, StateRecord>& records,
        std::uint64_t startKey,
        std::uint64_t goalKey,
        const std::vector<voxel_planner::BusbarPose>& poses,
        SearchOutcome& outcome) {
        std::vector<std::uint64_t> reversed;
        std::uint64_t current = goalKey;
        while (true) {
            reversed.push_back(current);
            if (current == startKey) {
                break;
            }
            const auto record = records.find(current);
            if (record == records.end() ||
                record->second.parent == kNoParent) {
                outcome.path.path.clear();
                outcome.path.pose_description.clear();
                return;
            }
            current = record->second.parent;
        }
        std::reverse(reversed.begin(), reversed.end());
        outcome.path.path.reserve(reversed.size());
        for (std::size_t index = 0U;
             index < reversed.size();
             ++index) {
            const std::uint64_t key = reversed[index];
            const std::size_t voxelIndex = key / poses.size();
            const std::uint32_t poseId = static_cast<std::uint32_t>(
                key % poses.size());
            const Point3D point = pointFromVoxelIndex(map, voxelIndex);
            outcome.path.path.push_back(point);
            appendConditionalPose(
                map,
                point,
                poses[poseId],
                index,
                outcome.path.pose_description);
        }
    }

    static void appendConditionalPose(
        const VoxelGrid& map,
        const Point3D& point,
        const voxel_planner::BusbarPose& pose,
        std::size_t waypointIndex,
        std::vector<Pose>& output) {
        if (!map.isValid(point.x, point.y, point.z) ||
            map.getState(point.x, point.y, point.z) !=
                VoxelState::POSE_CONDITIONAL) {
            return;
        }
        const Vector normal = unitNormal(pose);
        const Vector tangent = unitTangent(pose);
        output.push_back({
            {normal.x, normal.y, normal.z},
            {tangent.x, tangent.y, tangent.z},
            waypointIndex});
    }

    static double squaredDistance(
        const Point3D& lhs,
        const Point3D& rhs) noexcept {
        const double dx = static_cast<double>(lhs.x) - rhs.x;
        const double dy = static_cast<double>(lhs.y) - rhs.y;
        const double dz = static_cast<double>(lhs.z) - rhs.z;
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
        std::unordered_set<std::size_t> uniqueCandidate;
        std::size_t shared = 0U;
        for (const Point3D& point : candidate) {
            if (!map.isValid(point.x, point.y, point.z)) {
                continue;
            }
            const std::size_t index = map.index(
                static_cast<std::uint32_t>(point.x),
                static_cast<std::uint32_t>(point.y),
                static_cast<std::uint32_t>(point.z));
            if (uniqueCandidate.insert(index).second &&
                accepted.find(index) != accepted.end()) {
                ++shared;
            }
        }
        return uniqueCandidate.empty()
            ? 1.0
            : static_cast<double>(shared) /
                static_cast<double>(uniqueCandidate.size());
    }

    static void addSoftPathPenalty(
        const VoxelGrid& map,
        const std::vector<Point3D>& centerline,
        const Point3D& rawStart,
        const Point3D& rawGoal,
        float endpointProtectionRadius,
        CostPenaltyMap& penalties) {
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
                            squaredDistance(candidate, rawStart) <=
                                endpointRadiusSquared ||
                            squaredDistance(candidate, rawGoal) <=
                                endpointRadiusSquared) {
                            continue;
                        }
                        const std::size_t index = map.index(
                            static_cast<std::uint32_t>(candidate.x),
                            static_cast<std::uint32_t>(candidate.y),
                            static_cast<std::uint32_t>(candidate.z));
                        penalties.values[index] +=
                            dx == 0 && dy == 0 && dz == 0
                                ? kCenterlineReusePenalty
                                : kNeighborReusePenalty;
                    }
                }
            }
        }
    }
};

} // namespace module3_astar
