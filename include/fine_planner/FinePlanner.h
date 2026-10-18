#pragma once

#include "environment/SdfVolume.h"

#include <cstddef>
#include <string>
#include <vector>

namespace fine_planner {

using Vec3 = environment::Vec3;

enum class SegmentType { STRAIGHT, HORIZONTAL_BEND, VERTICAL_BEND };

// Bends require all exact parameters. Names identify the requested process;
// the explicit axis, rather than the name, defines the actual bend plane.
struct SegmentSpec {
    SegmentType type = SegmentType::STRAIGHT;
    Vec3 start{};
    Vec3 end{};
    Vec3 center{};
    Vec3 axis{};
    double radius = 0.0;
    double sweep_rad = 0.0;
};

struct Pose {
    Vec3 position{};
    Vec3 tangent{};
    Vec3 normal{};
    Vec3 binormal{};
};

struct Segment {
    SegmentSpec spec{};
    std::size_t start_index = 0;
    std::size_t end_index = 0;
};

struct Request {
    std::vector<Vec3> coarse_path;
    std::vector<SegmentSpec> segments; // Empty permits only a proven straight path.
    Vec3 start_tangent{1, 0, 0};
    Vec3 start_normal{0, 0, 1};
    Vec3 end_tangent{1, 0, 0};
    Vec3 end_normal{0, 0, 1};
    double width = 0.0;
    double thickness = 0.0;
    double safe_margin = 0.0;
    double minimum_bend_radius = 0.0;
    double maximum_sample_spacing = 0.25; // Native voxel units.
    std::string sample_id;
    std::string scene_id;
    std::string map_sha256;
};

struct ValidationReport {
    bool ok = false;
    std::string reason;
    double minimum_conservative_clearance = 0.0;
    double maximum_spacing = 0.0;
};

struct Result {
    bool success = false;
    std::string reason;
    std::vector<Pose> poses;
    std::vector<Segment> segments;
    ValidationReport validation;
};

class FinePlanner {
public:
    static Result plan(const Request& request, const environment::GlobalSdf& sdf);
};

} // namespace fine_planner
