#include "fine_planner/FinePlanner.h"
#include "fine_planner/GeometryValidator.h"
#include "GeometryOps.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fine_planner {
namespace {
constexpr double kPi = 3.14159265358979323846;

std::vector<SegmentSpec> specifications(const Request& request) {
    if (!request.segments.empty()) return request.segments;
    const Vec3 start=request.coarse_path.front(), end=request.coarse_path.back();
    const Vec3 direction=detail::unit(detail::sub(end,start));
    const double length=detail::distance(start,end);
    double previous=-1.0;
    for (Vec3 point: request.coarse_path) {
        const Vec3 offset=detail::sub(point,start);
        const double station=detail::dot(offset,direction);
        if (detail::norm(detail::sub(offset,detail::mul(direction,station)))>1e-4 ||
            station<previous-1e-4 || station<-1e-4 || station>length+1e-4)
            throw std::invalid_argument("INSUFFICIENT_GEOMETRY_FOR_BEND");
        previous=station;
    }
    return {{SegmentType::STRAIGHT,start,end,{},{},0.0,0.0}};
}

void appendSegment(Result& result, const SegmentSpec& spec, Vec3& normal,
                   double spacing) {
    const bool straight=spec.type==SegmentType::STRAIGHT;
    Vec3 tangent{}, axis{}, radial{};
    double length=0.0;
    if (straight) {
        tangent=detail::unit(detail::sub(spec.end,spec.start));
        length=detail::distance(spec.start,spec.end);
    } else {
        if (!std::isfinite(spec.radius)||spec.radius<=0.0 ||
            !std::isfinite(spec.sweep_rad)||std::abs(spec.sweep_rad)<1e-8||
            std::abs(spec.sweep_rad)>2*kPi)
            throw std::invalid_argument("INVALID_ARC_PARAMETERS");
        if (!detail::finite(spec.axis) || std::abs(detail::norm(spec.axis)-1.0)>1e-5)
            throw std::invalid_argument("ARC_AXIS_NOT_UNIT");
        axis=spec.axis;
        radial=detail::sub(spec.start,spec.center);
        if (std::abs(detail::norm(radial)-spec.radius)>1e-4 ||
            std::abs(detail::dot(radial,axis))>1e-4 ||
            !detail::near(detail::add(spec.center,detail::rotate(radial,axis,spec.sweep_rad)),spec.end))
            throw std::invalid_argument("ARC_ENDPOINT_MISMATCH");
        tangent=detail::unit(detail::mul(detail::cross(axis,radial),
                                        spec.sweep_rad>0.0?1.0:-1.0));
        length=std::abs(spec.sweep_rad)*spec.radius;
    }
    if (!result.poses.empty() &&
        (!detail::near(result.poses.back().position,spec.start) ||
         !detail::aligned(result.poses.back().tangent,tangent)))
        throw std::invalid_argument("SEGMENT_NOT_G1_CONTINUOUS");
    const std::size_t begin=result.poses.empty()?0U:result.poses.size()-1U;
    const double step=straight?spacing:std::min(spacing/2.0,spec.radius*0.1);
    const double countDouble=std::ceil(length/step);
    if (!std::isfinite(countDouble)||countDouble>1000000.0)
        throw std::invalid_argument("TOO_MANY_GEOMETRY_SAMPLES");
    const std::size_t count=std::max<std::size_t>(1U,static_cast<std::size_t>(countDouble));
    const Vec3 startNormal=normal;
    for (std::size_t i=result.poses.empty()?0U:1U; i<=count; ++i) {
        const double fraction=static_cast<double>(i)/static_cast<double>(count);
        const double angle=spec.sweep_rad*fraction;
        const Vec3 position=straight
            ?detail::add(spec.start,detail::mul(detail::sub(spec.end,spec.start),fraction))
            :detail::add(spec.center,detail::rotate(radial,axis,angle));
        const Vec3 currentTangent=straight?tangent:detail::rotate(tangent,axis,angle);
        const Vec3 currentNormal=straight?startNormal:detail::rotate(startNormal,axis,angle);
        result.poses.push_back(detail::makePose(position,currentTangent,currentNormal));
    }
    normal=result.poses.back().normal;
    result.segments.push_back({spec,begin,result.poses.size()-1U});
}
}

Result FinePlanner::plan(const Request& request, const environment::GlobalSdf& sdf) {
    Result result;
    try {
        if (request.coarse_path.size()<2U || request.width<=0.0 || request.thickness<=0.0 ||
            !std::isfinite(request.width)||!std::isfinite(request.thickness)||
            !std::isfinite(request.safe_margin)||request.safe_margin<0.0||
            !std::isfinite(request.minimum_bend_radius)||request.minimum_bend_radius<0.0||
            !std::isfinite(request.maximum_sample_spacing)||
            request.maximum_sample_spacing<=0.0||request.maximum_sample_spacing>0.25)
            throw std::invalid_argument("INVALID_REQUEST");
        for (Vec3 point:request.coarse_path)
            if (!detail::finite(point)) throw std::invalid_argument("NONFINITE_COARSE_PATH");
        const auto specs=specifications(request);
        if (!detail::near(specs.front().start,request.coarse_path.front())||
            !detail::near(specs.back().end,request.coarse_path.back()))
            throw std::invalid_argument("COARSE_ENDPOINT_MISMATCH");
        Vec3 normal=detail::unit(request.start_normal);
        for (const SegmentSpec& spec:specs) {
            if (spec.type!=SegmentType::STRAIGHT &&
                spec.radius+1e-9<request.minimum_bend_radius)
                throw std::invalid_argument("BEND_RADIUS_BELOW_MINIMUM");
            appendSegment(result,spec,normal,request.maximum_sample_spacing);
        }
        if (!detail::aligned(result.poses.front().tangent,request.start_tangent)||
            !detail::aligned(result.poses.back().tangent,request.end_tangent)||
            !detail::aligned(result.poses.front().normal,request.start_normal)||
            !detail::aligned(result.poses.back().normal,request.end_normal))
            throw std::invalid_argument("ENDPOINT_POSE_MISMATCH");
        result.validation=GeometryValidator::validate(request,result,sdf);
        result.success=result.validation.ok;
        result.reason=result.validation.reason;
    } catch (const std::exception& error) {
        result.success=false;
        result.reason=error.what();
    }
    return result;
}

} // namespace fine_planner
