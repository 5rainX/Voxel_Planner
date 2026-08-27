#include "Module3_AStar/JumpAStar.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace module3_astar {
namespace {

struct Direction {
    int dx = 0;
    int dy = 0;
    int dz = 0;
    double cost = 0.0;
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
                    dz,
                    std::sqrt(static_cast<double>(
                        dx * dx + dy * dy + dz * dz))});
            }
        }
    }
    return directions;
}

const std::vector<Direction>& directions() {
    static const std::vector<Direction> value = makeDirections();
    return value;
}

struct EdgeKey {
    std::size_t from = 0U;
    std::size_t to = 0U;

    bool operator==(const EdgeKey& other) const noexcept {
        return from == other.from && to == other.to;
    }
};

struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& edge) const noexcept {
        const std::size_t mixed =
            edge.from ^ (edge.to + static_cast<std::size_t>(0x9e3779b9U) +
                         (edge.from << 6U) + (edge.from >> 2U));
        return mixed;
    }
};

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

bool segmentAllowed(
    const VoxelGrid& map,
    const Point3D& from,
    const Point3D& to,
    const std::unordered_set<EdgeKey, EdgeKeyHash>& bannedEdges) {
    const int steps = std::max({
        std::abs(to.x - from.x),
        std::abs(to.y - from.y),
        std::abs(to.z - from.z)});
    if (steps <= 0) {
        return false;
    }

    Point3D current = from;
    for (int step = 0; step < steps; ++step) {
        Point3D next = current;
        next.x += (to.x > current.x) - (to.x < current.x);
        next.y += (to.y > current.y) - (to.y < current.y);
        next.z += (to.z > current.z) - (to.z < current.z);
        if (!map.isValid(next.x, next.y, next.z) ||
            map.isRawObstacle(pointIndex(map, next)) ||
            bannedEdges.find({
                pointIndex(map, current),
                pointIndex(map, next)}) != bannedEdges.end()) {
            return false;
        }
        current = next;
    }
    return current == to;
}

bool hasAdjacentRawObstacle(
    const VoxelGrid& map,
    const Point3D& point) {
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) {
                    continue;
                }
                const Point3D neighbor{
                    point.x + dx,
                    point.y + dy,
                    point.z + dz};
                if (!map.isValid(neighbor.x, neighbor.y, neighbor.z) ||
                    map.isRawObstacle(pointIndex(map, neighbor))) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool jumpCandidate(
    const VoxelGrid& map,
    const Point3D& origin,
    const Direction& direction,
    const Point3D& goal,
    const std::unordered_set<EdgeKey, EdgeKeyHash>& bannedEdges,
    Point3D& result) {
    Point3D current = origin;
    while (true) {
        const Point3D next{
            current.x + direction.dx,
            current.y + direction.dy,
            current.z + direction.dz};
        if (!map.isValid(next.x, next.y, next.z) ||
            map.isRawObstacle(pointIndex(map, next))) {
            return current == origin
                ? false
                : segmentAllowed(map, origin, current, bannedEdges) &&
                    (result = current, true);
        }
        current = next;
        if (current == goal || hasAdjacentRawObstacle(map, current)) {
            if (!segmentAllowed(map, origin, current, bannedEdges)) {
                return false;
            }
            result = current;
            return true;
        }
    }
}

std::vector<Point3D> runSearch(
    const VoxelGrid& map,
    const Point3D& start,
    const Point3D& goal,
    const std::unordered_set<EdgeKey, EdgeKeyHash>& bannedEdges) {
    if (!map.isValid(start.x, start.y, start.z) ||
        !map.isValid(goal.x, goal.y, goal.z) ||
        map.isRawObstacle(pointIndex(map, start)) ||
        map.isRawObstacle(pointIndex(map, goal))) {
        return {};
    }
    if (start == goal) {
        return {start};
    }

    const std::size_t count = map.voxelCount();
    const std::size_t startIndex = pointIndex(map, start);
    const std::size_t goalIndex = pointIndex(map, goal);
    std::vector<double> gScore(count, std::numeric_limits<double>::infinity());
    std::vector<std::size_t> parent(count, kNoIndex);
    std::vector<bool> closed(count, false);
    std::priority_queue<
        OpenNode,
        std::vector<OpenNode>,
        OpenNodeCompare> open;

    gScore[startIndex] = 0.0;
    open.push({startIndex, 0.0, distance(start, goal)});

    const auto relax = [&](std::size_t fromIndex,
                           const Point3D& from,
                           const Point3D& to,
                           std::priority_queue<
                               OpenNode,
                               std::vector<OpenNode>,
                               OpenNodeCompare>& queue) {
        if (!map.isValid(to.x, to.y, to.z) ||
            map.isRawObstacle(pointIndex(map, to))) {
            return;
        }
        const std::size_t targetIndex = pointIndex(map, to);
        if (closed[targetIndex] ||
            !segmentAllowed(map, from, to, bannedEdges)) {
            return;
        }
        const double tentative = gScore[fromIndex] + distance(from, to);
        if (tentative + 1e-12 >= gScore[targetIndex]) {
            return;
        }
        gScore[targetIndex] = tentative;
        parent[targetIndex] = fromIndex;
        queue.push({
            targetIndex,
            tentative,
            tentative + distance(to, goal)});
    };

    while (!open.empty()) {
        const OpenNode currentNode = open.top();
        open.pop();
        if (closed[currentNode.index] ||
            currentNode.g > gScore[currentNode.index] + 1e-12) {
            continue;
        }
        closed[currentNode.index] = true;
        if (currentNode.index == goalIndex) {
            break;
        }

        const Point3D current = pointFromIndex(map, currentNode.index);
        for (const Direction& direction : directions()) {
            const Point3D neighbor{
                current.x + direction.dx,
                current.y + direction.dy,
                current.z + direction.dz};
            relax(currentNode.index, current, neighbor, open);
        }

        for (const Direction& direction : directions()) {
            Point3D jumpPoint{};
            if (jumpCandidate(
                    map,
                    current,
                    direction,
                    goal,
                    bannedEdges,
                    jumpPoint) &&
                !(jumpPoint == current)) {
                relax(currentNode.index, current, jumpPoint, open);
            }
        }
    }

    if (!closed[goalIndex]) {
        return {};
    }

    std::vector<std::size_t> reversed;
    for (std::size_t current = goalIndex;
         current != kNoIndex;
         current = parent[current]) {
        reversed.push_back(current);
        if (current == startIndex) {
            break;
        }
    }
    if (reversed.empty() || reversed.back() != startIndex) {
        return {};
    }

    std::reverse(reversed.begin(), reversed.end());
    std::vector<Point3D> path;
    path.push_back(start);
    for (std::size_t i = 1U; i < reversed.size(); ++i) {
        const Point3D from = pointFromIndex(map, reversed[i - 1U]);
        const Point3D to = pointFromIndex(map, reversed[i]);
        Point3D current = from;
        while (!(current == to)) {
            current.x += (to.x > current.x) - (to.x < current.x);
            current.y += (to.y > current.y) - (to.y < current.y);
            current.z += (to.z > current.z) - (to.z < current.z);
            path.push_back(current);
        }
    }
    return path;
}

void banPathEdges(
    const VoxelGrid& map,
    const std::vector<Point3D>& path,
    std::unordered_set<EdgeKey, EdgeKeyHash>& bannedEdges) {
    for (std::size_t i = 1U; i < path.size(); ++i) {
        const std::size_t from = pointIndex(map, path[i - 1U]);
        const std::size_t to = pointIndex(map, path[i]);
        bannedEdges.insert({from, to});
        bannedEdges.insert({to, from});
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
    std::unordered_set<EdgeKey, EdgeKeyHash> bannedEdges;
    for (int pathIndex = 0; pathIndex < maxPaths; ++pathIndex) {
        std::vector<Point3D> path =
            runSearch(map, start, goal, bannedEdges);
        if (path.empty()) {
            break;
        }
        paths.push_back(path);
        banPathEdges(map, path, bannedEdges);
    }
    return paths;
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
