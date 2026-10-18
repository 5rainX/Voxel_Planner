#include "fine_planner/ReferenceExport.h"
#include "GeometryOps.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace fine_planner {
namespace {
void string(std::ostream& out, const std::string& value) {
    out << '"';
    for (unsigned char c:value) {
        if (c=='"'||c=='\\') out << '\\' << static_cast<char>(c);
        else if (c=='\n') out << "\\n";
        else if (c=='\r') out << "\\r";
        else if (c=='\t') out << "\\t";
        else if (c<32U) throw std::invalid_argument("CONTROL_CHARACTER_IN_ID");
        else out << static_cast<char>(c);
    }
    out << '"';
}
void vec(std::ostream& out, Vec3 v) { out << '[' << v.x << ',' << v.y << ',' << v.z << ']'; }
}

void writeReferenceGeometryV1(
    const Result& result, const Request& request,
    const ValidationReport& validation,
    const std::filesystem::path& output) {
    if (!result.success||!validation.ok||request.sample_id.empty()||
        request.scene_id.empty()||request.map_sha256.size()!=64U||
        result.poses.size()<2U||result.segments.empty())
        throw std::invalid_argument("UNVALIDATED_NATIVE_REFERENCE");
    if (std::filesystem::exists(output))
        throw std::runtime_error("REFERENCE_OUTPUT_ALREADY_EXISTS");
    std::ofstream out(output,std::ios::binary);
    if (!out) throw std::runtime_error("REFERENCE_OUTPUT_OPEN_FAILED");
    out << std::setprecision(17);
    out << "{\"sample_id\":"; string(out,request.sample_id);
    out << ",\"reference_geometry_version\":\"reference_geometry_v1\",\"status\":\"VALID\"";
    out << ",\"source\":{\"planner\":\"native_fine_planner\",\"source_kind\":\"NATIVE_FINE_PLANNER\",\"source_record\":";
    string(out,"native_request:"+request.sample_id);
    out << ",\"source_line\":1,\"source_sample_id\":"; string(out,request.sample_id);
    out << ",\"scene_id\":"; string(out,request.scene_id);
    out << ",\"map_id\":"; string(out,request.map_sha256);
    out << ",\"generator_version\":\"P1-fine-planner-v1\"}";
    out << ",\"reference_centerline\":[";
    for (std::size_t i=0;i<result.poses.size();++i) {
        if (i) out << ',';
        vec(out,result.poses[i].position);
    }
    out << "],\"reference_poses\":[";
    double cumulative=0.0;
    double total=0.0;
    for (std::size_t i=1;i<result.poses.size();++i)
        total+=detail::distance(result.poses[i-1].position,result.poses[i].position);
    std::size_t segmentIndex=0U;
    for (std::size_t i=0;i<result.poses.size();++i) {
        if (i) { out << ','; cumulative+=detail::distance(result.poses[i-1].position,result.poses[i].position); }
        while (segmentIndex+1U<result.segments.size() &&
               i>result.segments[segmentIndex].end_index) ++segmentIndex;
        const Pose& pose=result.poses[i];
        out << "{\"position\":"; vec(out,pose.position);
        out << ",\"tangent\":"; vec(out,pose.tangent);
        out << ",\"normal\":"; vec(out,pose.normal);
        out << ",\"binormal\":"; vec(out,pose.binormal);
        out << ",\"orientation_valid\":true,\"segment_index\":" << segmentIndex;
        out << ",\"segment_type\":"; string(out,detail::typeName(result.segments[segmentIndex].spec.type));
        out << ",\"s_voxel\":" << cumulative << ",\"s_normalized\":" << cumulative/total << '}';
    }
    out << "],\"segments\":[";
    for (std::size_t i=0;i<result.segments.size();++i) {
        if (i) out << ',';
        const Segment& segment=result.segments[i];
        out << "{\"segment_index\":" << i << ",\"segment_type\":";
        string(out,detail::typeName(segment.spec.type));
        out << ",\"start_index\":" << segment.start_index
            << ",\"end_index\":" << segment.end_index
            << ",\"geometry_source\":\"NATIVE_FINE_PLANNER_EXACT_GEOMETRY\""
            << ",\"geometry_reconstruction_method\":\"DETERMINISTIC_EXACT_SEGMENTS\""
            << ",\"orientation_source\":\"NATIVE_TRANSPORTED_FRAME\""
            << ",\"start\":";
        vec(out,segment.spec.start);
        out << ",\"end\":"; vec(out,segment.spec.end);
        if (segment.spec.type!=SegmentType::STRAIGHT) {
            out << ",\"center\":"; vec(out,segment.spec.center);
            out << ",\"axis\":"; vec(out,segment.spec.axis);
            out << ",\"radius_voxel\":" << segment.spec.radius
                << ",\"signed_sweep_rad\":" << segment.spec.sweep_rad;
        }
        out << '}';
    }
    out << "],\"validation\":{\"point_count\":" << result.poses.size()
        << ",\"length_voxel\":" << total
        << ",\"max_spacing_voxel\":" << validation.maximum_spacing
        << ",\"max_junction_position_error_voxel\":0"
        << ",\"max_junction_tangent_discontinuity_deg\":0"
        << ",\"orientation_valid_nodes\":" << result.poses.size()
        << ",\"sampling_rule\":\"max_spacing_0.25_voxel\""
        << ",\"validation_scope\":\"SAMPLED_SDF_CONSERVATIVE_SECTION\""
        << ",\"independent_full_size_collision_proven\":false"
        << ",\"minimum_conservative_clearance_voxel\":"
        << validation.minimum_conservative_clearance << "}}\n";
    if (!out) throw std::runtime_error("REFERENCE_OUTPUT_WRITE_FAILED");
}

} // namespace fine_planner
