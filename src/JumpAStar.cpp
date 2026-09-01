#include "Module3_AStar/JumpAStar.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>

namespace module3_astar {
namespace {

struct Direction {
    int dx = 0;
    int dy = 0;
    int dz = 0;
};

std::vector<Direction> makeDirections() {
    std::vector<Direction> directions;
    directions.reserve(26U);
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) {
                    continue;
                }
                directions.push_back({
                    dx,
                    dy,
                    dz});
            }
        }
    }
    return directions;
}

const std::vector<Direction>& directions() {
    static const std::vector<Direction> value = makeDirections();
    return value;
}

struct OpenNode {
    std::size_t index = 0U;
    double g = std::numeric_limits<double>::infinity();
    double f = std::numeric_limits<double>::infinity();
};

struct OpenNodeCompare {
    bool operator()(const OpenNode& lhs, const OpenNode& rhs) const noexcept {
        if (lhs.f != rhs.f) {
            return lhs.f > rhs.f;
        }
        return lhs.g < rhs.g;
    }
};

constexpr std::size_t kNoIndex =
    std::numeric_limits<std::size_t>::max();
constexpr double kHeuristicWeight = 2.0;
constexpr double kTieBreakScale = 1.001;

/**
 * @brief Dense position-only search records.
 *
 * JumpAStar has exactly one search state per voxel. Keeping the records in
 * direct-addressed arrays removes hashing from the inner relaxation loop
 * without changing the queue ordering or any search decision.
 */
struct FlatSearchRecords {
    explicit FlatSearchRecords(std::size_t voxelCount)
        : g(voxelCount),
          parent(voxelCount),
          closed(voxelCount, 0U) {
#pragma omp parallel sections
        {
#pragma omp section
            {
                std::fill(
                    g.begin(),
                    g.end(),
                    std::numeric_limits<double>::infinity());
            }
#pragma omp section
            {
                std::fill(parent.begin(), parent.end(), kNoIndex);
            }
        }
    }

    bool generated(std::size_t index) const noexcept {
        return index < g.size() &&
            g[index] != std::numeric_limits<double>::infinity();
    }

    std::vector<double> g;
    std::vector<std::size_t> parent;
    std::vector<std::uint8_t> closed;
    std::size_t generatedCount = 0U;
};

constexpr int kMaxPathsPerRequest = 10;
constexpr double kScoreEpsilon = 1e-12;
constexpr int kHardBlockRadius = 18;
constexpr int kMinimumHardBlockRadius = 4;
constexpr int kEndpointProtectionRadius = 3;
constexpr std::size_t kProtectedPathPrefix = 4U;
constexpr std::size_t kProtectedPathSuffix = 4U;

Point3D pointFromIndex(const VoxelGrid& map, std::size_t index) {
    const std::size_t plane =
        static_cast<std::size_t>(map.width()) * map.height();
    const int z = static_cast<int>(index / plane);
    const std::size_t remainder = index % plane;
    const int y = static_cast<int>(remainder / map.width());
    const int x = static_cast<int>(remainder % map.width());
    return {x, y, z};
}

std::size_t pointIndex(const VoxelGrid& map, const Point3D& point) {
    return map.index(
        static_cast<std::uint32_t>(point.x),
        static_cast<std::uint32_t>(point.y),
        static_cast<std::uint32_t>(point.z));
}

double distance(const Point3D& lhs, const Point3D& rhs) {
    const double dx = static_cast<double>(lhs.x - rhs.x);
    const double dy = static_cast<double>(lhs.y - rhs.y);
    const double dz = static_cast<double>(lhs.z - rhs.z);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double weightedHeuristic(
    const Point3D& point,
    const Point3D& goal) noexcept {
    const double originalH = distance(point, goal);
    const double tieBrokenH = originalH * kTieBreakScale;
    return kHeuristicWeight * tieBrokenH;
}

bool isHardBlocked(
    const std::vector<std::uint8_t>& hardBlocked,
    std::size_t index) noexcept {
    return index < hardBlocked.size() && hardBlocked[index] != 0U;
}

bool segmentAllowed(
    const VoxelGrid& map,
    const Point3D& from,
    const Point3D& to,
    const std::vector<std::uint8_t>& hardBlocked) {
    if (!map.isValid(from.x, from.y, from.z) ||
        !map.isValid(to.x, to.y, to.z)) {
        return false;
    }

    const int steps = std::max({
        std::abs(to.x - from.x),
        std::abs(to.y - from.y),
        std::abs(to.z - from.z)});
    if (steps <= 0) {
        return false;
    }

    Point3D current = from;
    for (int step = 0; step < steps; ++step) {
        const Point3D previous = current;
        Point3D next = current;
        next.x += (to.x > current.x) - (to.x < current.x);
        next.y += (to.y > current.y) - (to.y < current.y);
        next.z += (to.z > current.z) - (to.z < current.z);
        if (next == previous) {
            return false;
        }
        if (!map.isValid(next.x, next.y, next.z) ||
            map.isRawObstacle(pointIndex(map, next)) ||
            isHardBlocked(hardBlocked, pointIndex(map, next))) {
            return false;
        }
        current = next;
    }
    return current == to;
}

bool appendSegmentPath(
    const Point3D& to,
    std::vector<Point3D>& path) {
    if (path.empty()) {
        path.push_back(to);
        return true;
    }

    const Point3D from = path.back();
    const int steps = std::max({
        std::abs(to.x - from.x),
        std::abs(to.y - from.y),
        std::abs(to.z - from.z)});
    for (int step = 0; step < steps; ++step) {
        Point3D current = path.back();
        const Point3D previous = current;
        current.x += (to.x > current.x) - (to.x < current.x);
        current.y += (to.y > current.y) - (to.y < current.y);
        current.z += (to.z > current.z) - (to.z < current.z);
        if (current == previous) {
            return false;
        }
        path.push_back(current);
    }
    return path.back() == to;
}

std::vector<Point3D> makeSegmentPath(
    const Point3D& from,
    const Point3D& to) {
    std::vector<Point3D> path;
    path.push_back(from);
    if (!appendSegmentPath(to, path)) {
        path.clear();
    }
    return path;
}

std::vector<Point3D> runSearch(
    const VoxelGrid& map,
    const Point3D& start,
    const Point3D& goal,
    const std::vector<std::uint8_t>& hardBlocked) {
    if (!map.isValid(start.x, start.y, start.z) ||
        !map.isValid(goal.x, goal.y, goal.z)) {
        return {};
    }

    const std::size_t startIndex = pointIndex(map, start);
    const std::size_t goalIndex = pointIndex(map, goal);
    if (map.isRawObstacle(startIndex) ||
        map.isRawObstacle(goalIndex) ||
        isHardBlocked(hardBlocked, startIndex) ||
        isHardBlocked(hardBlocked, goalIndex)) {
        return {};
    }
    if (start == goal) {
        return {start};
    }
    if (segmentAllowed(map, start, goal, hardBlocked)) {
        return makeSegmentPath(start, goal);
    }

    FlatSearchRecords records(map.voxelCount());
    std::priority_queue<
        OpenNode,
        std::vector<OpenNode>,
        OpenNodeCompare> open;

    records.g[startIndex] = 0.0;
    records.parent[startIndex] = kNoIndex;
    records.closed[startIndex] = 0U;
    records.generatedCount = 1U;
    open.push({
        startIndex,
        0.0,
        weightedHeuristic(start, goal)});

    const auto relax = [&](std::size_t fromIndex,
                           const Point3D& to,
                           std::priority_queue<
                               OpenNode,
                               std::vector<OpenNode>,
                               OpenNodeCompare>& queue) {
        if (!map.isValid(to.x, to.y, to.z) ||
            map.isRawObstacle(pointIndex(map, to))) {
            return;
        }
        if (!records.generated(fromIndex)) {
            return;
        }
        const std::size_t targetIndex = pointIndex(map, to);
        if (isHardBlocked(hardBlocked, targetIndex)) {
            return;
        }
        const double tentative =
            records.g[fromIndex] +
            distance(pointFromIndex(map, fromIndex), to);
        if (records.generated(targetIndex) &&
            (records.closed[targetIndex] != 0U ||
             tentative + kScoreEpsilon >= records.g[targetIndex])) {
            return;
        }
        if (!records.generated(targetIndex)) {
            ++records.generatedCount;
        }
        records.g[targetIndex] = tentative;
        records.parent[targetIndex] = fromIndex;
        records.closed[targetIndex] = 0U;
        queue.push({
            targetIndex,
            tentative,
            tentative + weightedHeuristic(to, goal)});
    };

    while (!open.empty()) {
        const OpenNode currentNode = open.top();
        open.pop();
        if (!records.generated(currentNode.index) ||
            records.closed[currentNode.index] != 0U ||
            currentNode.g >
                records.g[currentNode.index] + kScoreEpsilon) {
            continue;
        }
        records.closed[currentNode.index] = 1U;
        if (currentNode.index == goalIndex) {
            break;
        }

        const Point3D current = pointFromIndex(map, currentNode.index);
        for (const Direction& direction : directions()) {
            const std::int64_t neighborX =
                static_cast<std::int64_t>(current.x) + direction.dx;
            const std::int64_t neighborY =
                static_cast<std::int64_t>(current.y) + direction.dy;
            const std::int64_t neighborZ =
                static_cast<std::int64_t>(current.z) + direction.dz;
            if (neighborX < 0 || neighborY < 0 || neighborZ < 0 ||
                neighborX >= static_cast<std::int64_t>(map.width()) ||
                neighborY >= static_cast<std::int64_t>(map.height()) ||
                neighborZ >= static_cast<std::int64_t>(map.depth())) {
                continue;
            }
            const Point3D neighbor{
                static_cast<int>(neighborX),
                static_cast<int>(neighborY),
                static_cast<int>(neighborZ)};
            relax(currentNode.index, neighbor, open);
        }
    }

    if (!records.generated(goalIndex) ||
        records.closed[goalIndex] == 0U) {
        return {};
    }

    std::vector<std::size_t> reversed;
    std::size_t currentIndex = goalIndex;
    for (std::size_t hops = 0U;
         currentIndex != kNoIndex && hops <= records.generatedCount;
         ++hops) {
        reversed.push_back(currentIndex);
        if (currentIndex == startIndex) {
            break;
        }
        if (!records.generated(currentIndex)) {
            return {};
        }
        currentIndex = records.parent[currentIndex];
    }
    if (reversed.empty() || reversed.back() != startIndex) {
        return {};
    }

    std::reverse(reversed.begin(), reversed.end());
    std::vector<Point3D> path;
    path.push_back(start);
    for (std::size_t i = 1U; i < reversed.size(); ++i) {
        const Point3D to = pointFromIndex(map, reversed[i]);
        if (!appendSegmentPath(to, path)) {
            return {};
        }
    }
    return path;
}

bool nearProtectedEndpoint(
    const Point3D& point,
    const Point3D& start,
    const Point3D& goal) {
    const auto insideRadius = [](const Point3D& lhs, const Point3D& rhs) {
        return std::max({
            std::abs(lhs.x - rhs.x),
            std::abs(lhs.y - rhs.y),
            std::abs(lhs.z - rhs.z)}) <= kEndpointProtectionRadius;
    };
    return insideRadius(point, start) || insideRadius(point, goal);
}

void hardBlockCube(
    const VoxelGrid& map,
    const Point3D& center,
    const Point3D& start,
    const Point3D& goal,
    int radius,
    std::vector<std::uint8_t>& hardBlocked) {
    for (int dz = -radius; dz <= radius; ++dz) {
        for (int dy = -radius; dy <= radius; ++dy) {
            for (int dx = -radius; dx <= radius; ++dx) {
                const Point3D blocked{
                    center.x + dx,
                    center.y + dy,
                    center.z + dz};
                if (!map.isValid(blocked.x, blocked.y, blocked.z) ||
                    nearProtectedEndpoint(blocked, start, goal)) {
                    continue;
                }
                hardBlocked[pointIndex(map, blocked)] = 1U;
            }
        }
    }
}

void hardBlockPathCorridor(
    const VoxelGrid& map,
    const std::vector<Point3D>& path,
    const Point3D& start,
    const Point3D& goal,
    int radius,
    std::vector<std::uint8_t>& hardBlocked) {
    if (hardBlocked.size() != map.voxelCount() ||
        path.size() <= 2U ||
        radius <= 0) {
        return;
    }
    std::size_t begin = std::min(kProtectedPathPrefix, path.size() - 1U);
    std::size_t end = path.size() > kProtectedPathSuffix
        ? path.size() - kProtectedPathSuffix
        : begin;
    if (begin >= end) {
        begin = path.size() / 2U;
        end = begin + 1U;
    }
    const std::size_t centerStride =
        static_cast<std::size_t>(std::max(1, radius));
    for (std::size_t pathIndex = begin;
         pathIndex < end;
         pathIndex += centerStride) {
        hardBlockCube(
            map,
            path[pathIndex],
            start,
            goal,
            radius,
            hardBlocked);
    }
    if (end > begin && (end - 1U - begin) % centerStride != 0U) {
        hardBlockCube(
            map,
            path[end - 1U],
            start,
            goal,
            radius,
            hardBlocked);
    }
}

bool samePath(
    const std::vector<Point3D>& lhs,
    const std::vector<Point3D>& rhs) {
    return lhs.size() == rhs.size() &&
        std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

bool alreadyReturned(
    const std::vector<std::vector<Point3D>>& paths,
    const std::vector<Point3D>& candidate) {
    for (const std::vector<Point3D>& path : paths) {
        if (samePath(path, candidate)) {
            return true;
        }
    }
    return false;
}

void rebuildHardBlockMask(
    const VoxelGrid& map,
    const std::vector<std::vector<Point3D>>& paths,
    const Point3D& start,
    const Point3D& goal,
    int radius,
    std::vector<std::uint8_t>& hardBlocked) {
    std::fill(hardBlocked.begin(), hardBlocked.end(), 0U);
    for (const std::vector<Point3D>& path : paths) {
        hardBlockPathCorridor(
            map,
            path,
            start,
            goal,
            radius,
            hardBlocked);
    }
}

} // namespace

std::vector<std::vector<Point3D>> JumpAStar::findPaths(
    const VoxelGrid& map,
    Point3D start,
    Point3D goal,
    int maxPaths) {
    if (maxPaths <= 0) {
        throw std::invalid_argument("maxPaths must be positive.");
    }

    std::vector<std::vector<Point3D>> paths;
    std::vector<std::uint8_t> hardBlocked(map.voxelCount(), 0U);
    const int boundedMaxPaths = std::min(maxPaths, kMaxPathsPerRequest);
    int blockadeRadius = kHardBlockRadius;
    while (static_cast<int>(paths.size()) < boundedMaxPaths) {
        rebuildHardBlockMask(
            map,
            paths,
            start,
            goal,
            blockadeRadius,
            hardBlocked);
        std::vector<Point3D> path =
            runSearch(map, start, goal, hardBlocked);
        if (path.empty()) {
            if (!paths.empty() && blockadeRadius > kMinimumHardBlockRadius) {
                blockadeRadius = std::max(
                    kMinimumHardBlockRadius,
                    blockadeRadius / 2);
                continue;
            }
            break;
        }
        if (alreadyReturned(paths, path)) {
            break;
        }
        paths.push_back(std::move(path));
    }
    return paths;
}

std::vector<std::vector<Point3D>> JumpAStar::findPaths(
    const VoxelGrid& map,
    Point3D start,
    Point3D goal,
    int maxPaths,
    const voxel_planner::EndpointPose&,
    const voxel_planner::EndpointPose&) {
    return findPaths(map, start, goal, maxPaths);
}

std::vector<Point3D> JumpAStar::findPath(
    const VoxelGrid& map,
    Point3D start,
    Point3D goal) {
    std::vector<std::vector<Point3D>> paths =
        findPaths(map, start, goal, 1);
    return paths.empty() ? std::vector<Point3D>{} : std::move(paths.front());
}

} // namespace module3_astar
