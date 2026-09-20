#pragma once

#include "VoxelPlannerAPI.h"
#include "environment/SdfVolume.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace optimizer {

struct TrajectoryNode {
    environment::Vec3 position{};
    environment::Vec3 tangent{1.0, 0.0, 0.0};
    environment::Vec3 normal{0.0, 0.0, 1.0};
};

class PoseTrajectory {
public:
    std::size_t size() const noexcept { return nodes_.size(); }
    bool empty() const noexcept { return nodes_.empty(); }

    const TrajectoryNode& operator[](std::size_t index) const {
        return nodes_[index];
    }

    TrajectoryNode& operator[](std::size_t index) {
        return nodes_[index];
    }

    const std::vector<TrajectoryNode>& nodes() const noexcept {
        return nodes_;
    }

    std::vector<TrajectoryNode>& nodes() noexcept {
        return nodes_;
    }

    void pushBack(const TrajectoryNode& node) {
        nodes_.push_back(node);
    }

    void recomputeFrames() {
        if (nodes_.empty()) {
            return;
        }
        if (nodes_.size() == 1U) {
            nodes_.front().tangent = {1.0, 0.0, 0.0};
            nodes_.front().normal = {0.0, 0.0, 1.0};
            return;
        }

        for (std::size_t index = 0U; index < nodes_.size(); ++index) {
            environment::Vec3 tangent{};
            if (index == 0U) {
                tangent = subtract(
                    nodes_[1U].position,
                    nodes_[0U].position);
            } else if (index + 1U == nodes_.size()) {
                tangent = subtract(
                    nodes_[index].position,
                    nodes_[index - 1U].position);
            } else {
                tangent = subtract(
                    nodes_[index + 1U].position,
                    nodes_[index - 1U].position);
            }
            nodes_[index].tangent = normalizedOrFallback(
                tangent,
                {1.0, 0.0, 0.0});
            nodes_[index].normal = stableNormal(nodes_[index].tangent);
        }
    }

    double polylineLength() const noexcept {
        double length = 0.0;
        for (std::size_t index = 1U; index < nodes_.size(); ++index) {
            length += norm(subtract(
                nodes_[index].position,
                nodes_[index - 1U].position));
        }
        return length;
    }

    static environment::Vec3 add(
        const environment::Vec3& lhs,
        const environment::Vec3& rhs) noexcept {
        return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
    }

    static environment::Vec3 subtract(
        const environment::Vec3& lhs,
        const environment::Vec3& rhs) noexcept {
        return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
    }

    static environment::Vec3 scale(
        const environment::Vec3& value,
        double factor) noexcept {
        return {value.x * factor, value.y * factor, value.z * factor};
    }

    static double dot(
        const environment::Vec3& lhs,
        const environment::Vec3& rhs) noexcept {
        return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
    }

    static environment::Vec3 cross(
        const environment::Vec3& lhs,
        const environment::Vec3& rhs) noexcept {
        return {
            lhs.y * rhs.z - lhs.z * rhs.y,
            lhs.z * rhs.x - lhs.x * rhs.z,
            lhs.x * rhs.y - lhs.y * rhs.x};
    }

    static double squaredNorm(const environment::Vec3& value) noexcept {
        return dot(value, value);
    }

    static double norm(const environment::Vec3& value) noexcept {
        return std::sqrt(squaredNorm(value));
    }

    static environment::Vec3 normalizedOrFallback(
        const environment::Vec3& value,
        const environment::Vec3& fallback) noexcept {
        const double length = norm(value);
        if (length <= 1.0e-9 || !std::isfinite(length)) {
            return fallback;
        }
        return scale(value, 1.0 / length);
    }

    static environment::Vec3 stableNormal(
        const environment::Vec3& tangent) noexcept {
        const environment::Vec3 absT{
            std::abs(tangent.x),
            std::abs(tangent.y),
            std::abs(tangent.z)};
        const environment::Vec3 reference =
            absT.x <= absT.y && absT.x <= absT.z
                ? environment::Vec3{1.0, 0.0, 0.0}
                : (absT.y <= absT.z
                    ? environment::Vec3{0.0, 1.0, 0.0}
                    : environment::Vec3{0.0, 0.0, 1.0});
        return normalizedOrFallback(
            cross(tangent, reference),
            {0.0, 0.0, 1.0});
    }

private:
    std::vector<TrajectoryNode> nodes_;
};

struct RefinementOptions {
    int smoothing_iterations = 8;
    double smoothing_factor = 0.35;
};

class ConfigurationSpaceRefiner {
public:
    explicit ConfigurationSpaceRefiner(RefinementOptions options = {})
        : options_(options) {}

    PoseTrajectory fromVoxelPath(
        const std::vector<voxel_planner::Point3D>& path,
        bool smooth = true) const {
        PoseTrajectory trajectory;
        trajectory.nodes().reserve(path.size());
        for (const voxel_planner::Point3D& point : path) {
            trajectory.pushBack({
                {static_cast<double>(point.x),
                 static_cast<double>(point.y),
                 static_cast<double>(point.z)},
                {},
                {}});
        }
        if (smooth) {
            smoothInterior(trajectory);
        }
        trajectory.recomputeFrames();
        return trajectory;
    }

    std::vector<voxel_planner::Point3D> toVoxelPath(
        const PoseTrajectory& trajectory) const {
        std::vector<voxel_planner::Point3D> path;
        for (const TrajectoryNode& node : trajectory.nodes()) {
            const voxel_planner::Point3D rounded{
                roundToInt(node.position.x),
                roundToInt(node.position.y),
                roundToInt(node.position.z)};
            appendContinuous26(path, rounded);
        }
        return path;
    }

private:
    static int roundToInt(double value) noexcept {
        if (value <= static_cast<double>(std::numeric_limits<int>::min())) {
            return std::numeric_limits<int>::min();
        }
        if (value >= static_cast<double>(std::numeric_limits<int>::max())) {
            return std::numeric_limits<int>::max();
        }
        return static_cast<int>(std::floor(value + 0.5));
    }

    static int stepToward(int current, int target) noexcept {
        return (target > current) - (target < current);
    }

    static void appendContinuous26(
        std::vector<voxel_planner::Point3D>& path,
        const voxel_planner::Point3D& target) {
        if (path.empty()) {
            path.push_back(target);
            return;
        }
        while (!(path.back() == target)) {
            voxel_planner::Point3D next = path.back();
            next.x += stepToward(next.x, target.x);
            next.y += stepToward(next.y, target.y);
            next.z += stepToward(next.z, target.z);
            if (next == path.back()) {
                break;
            }
            path.push_back(next);
        }
    }

    void smoothInterior(PoseTrajectory& trajectory) const {
        if (trajectory.size() <= 2U ||
            options_.smoothing_iterations <= 0 ||
            options_.smoothing_factor <= 0.0) {
            return;
        }
        const double alpha = std::clamp(
            options_.smoothing_factor,
            0.0,
            0.95);
        std::vector<environment::Vec3> updated(trajectory.size());
        for (int iteration = 0;
             iteration < options_.smoothing_iterations;
             ++iteration) {
            for (std::size_t index = 0U;
                 index < trajectory.size();
                 ++index) {
                updated[index] = trajectory[index].position;
            }
            for (std::size_t index = 1U;
                 index + 1U < trajectory.size();
                 ++index) {
                const environment::Vec3 average = PoseTrajectory::scale(
                    PoseTrajectory::add(
                        trajectory[index - 1U].position,
                        trajectory[index + 1U].position),
                    0.5);
                updated[index] = PoseTrajectory::add(
                    PoseTrajectory::scale(
                        trajectory[index].position,
                        1.0 - alpha),
                    PoseTrajectory::scale(average, alpha));
            }
            for (std::size_t index = 1U;
                 index + 1U < trajectory.size();
                 ++index) {
                trajectory[index].position = updated[index];
            }
        }
    }

    RefinementOptions options_{};
};

} // namespace optimizer
