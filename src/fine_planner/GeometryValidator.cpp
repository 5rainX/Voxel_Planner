#include "fine_planner/GeometryValidator.h"
#include "GeometryOps.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fine_planner {
namespace {
ValidationReport reject(const char* reason) { return {false,reason,0.0,0.0}; }
}

ValidationReport GeometryValidator::validate(
    const Request& request, const Result& result,
    const environment::GlobalSdf& sdf) {
    try {
        if (result.poses.size()<2U || result.segments.empty() ||
            request.coarse_path.size()<2U || request.width<=0.0 ||
            request.thickness<=0.0 || request.safe_margin<0.0 ||
            request.minimum_bend_radius<0.0 || request.maximum_sample_spacing<=0.0 ||
            request.maximum_sample_spacing>0.25 || !std::isfinite(request.width) ||
            !std::isfinite(request.thickness) || !std::isfinite(request.safe_margin) ||
            !std::isfinite(request.minimum_bend_radius) ||
            !std::isfinite(request.maximum_sample_spacing) ||
            !std::isfinite(request.maximum_twist_step_rad) ||
            request.maximum_twist_step_rad<=0.0 ||
            request.maximum_twist_step_rad>1.5707963267948966 ||
            !std::isfinite(request.minimum_twist_length) ||
            request.minimum_twist_length<0.0 ||
            !std::isfinite(request.maximum_twist_rate_rad_per_voxel) ||
            request.maximum_twist_rate_rad_per_voxel<0.0)
            return reject("INVALID_GEOMETRY");
        for (Vec3 point:request.coarse_path)
            if (!detail::finite(point)) return reject("NONFINITE_COARSE_PATH");
        if (request.segments.empty()) {
            if (result.segments.size()!=1U ||
                result.segments.front().spec.type!=SegmentType::STRAIGHT)
                return reject("INSUFFICIENT_GEOMETRY_FOR_BEND");
            const Vec3 start=request.coarse_path.front(), end=request.coarse_path.back();
            const Vec3 direction=detail::unit(detail::sub(end,start));
            const double length=detail::distance(start,end);
            double previous=-1.0;
            for (Vec3 point:request.coarse_path) {
                const Vec3 offset=detail::sub(point,start);
                const double station=detail::dot(offset,direction);
                if (detail::norm(detail::sub(offset,detail::mul(direction,station)))>1e-4 ||
                    station<previous-1e-4 || station<-1e-4 || station>length+1e-4)
                    return reject("INSUFFICIENT_GEOMETRY_FOR_BEND");
                previous=station;
            }
        }
        if (!detail::near(result.poses.front().position,request.coarse_path.front())||
            !detail::near(result.poses.back().position,request.coarse_path.back())||
            !detail::aligned(result.poses.front().tangent,request.start_tangent)||
            !detail::aligned(result.poses.back().tangent,request.end_tangent)||
            !detail::aligned(result.poses.front().normal,request.start_normal)||
            !detail::aligned(result.poses.back().normal,request.end_normal))
            return reject("ENDPOINT_POSE_MISMATCH");
        if (!request.segments.empty()) {
            if (request.segments.size()!=result.segments.size())
                return reject("REQUEST_SEGMENT_MISMATCH");
            for (std::size_t i=0;i<request.segments.size();++i) {
                const SegmentSpec& requested=request.segments[i];
                const SegmentSpec& actual=result.segments[i].spec;
                if (requested.type!=actual.type||
                    !detail::near(requested.start,actual.start)||
                    !detail::near(requested.end,actual.end)||
                    ((requested.type==SegmentType::HORIZONTAL_BEND ||
                      requested.type==SegmentType::VERTICAL_BEND) &&
                     (!detail::near(requested.center,actual.center)||
                      !detail::near(requested.axis,actual.axis)||
                      std::abs(requested.radius-actual.radius)>1e-8||
                      std::abs(requested.sweep_rad-actual.sweep_rad)>1e-8)) ||
                    (requested.type==SegmentType::STRAIGHT_TWIST &&
                     (!std::isfinite(requested.twist_angle_rad)||
                      !std::isfinite(requested.twist_length)||
                      !detail::near(requested.twist_start_normal,actual.twist_start_normal)||
                      !detail::near(requested.twist_end_normal,actual.twist_end_normal)||
                      std::abs(requested.twist_angle_rad-actual.twist_angle_rad)>1e-8||
                      std::abs(requested.twist_length-actual.twist_length)>1e-8||
                      requested.twist_profile!=actual.twist_profile)))
                    return reject("REQUEST_SEGMENT_MISMATCH");
            }
        }
        ValidationReport report;
        report.minimum_conservative_clearance=std::numeric_limits<double>::infinity();
        const double extent=0.5*std::hypot(request.width,request.thickness);
        if (!std::isfinite(extent)) return reject("INVALID_DIMENSIONS");
        for (std::size_t i=0;i<result.poses.size();++i) {
            const Pose& pose=result.poses[i];
            if (!detail::finite(pose.position)||!detail::finite(pose.tangent)||
                !detail::finite(pose.normal)||!detail::finite(pose.binormal)||
                std::abs(detail::norm(pose.tangent)-1.0)>1e-5||
                std::abs(detail::norm(pose.normal)-1.0)>1e-5||
                std::abs(detail::norm(pose.binormal)-1.0)>1e-5||
                std::abs(detail::dot(pose.tangent,pose.normal))>1e-5||
                std::abs(detail::dot(pose.tangent,pose.binormal))>1e-5||
                std::abs(detail::dot(pose.normal,pose.binormal))>1e-5||
                !detail::near(detail::cross(pose.tangent,pose.normal),pose.binormal,1e-5))
                return reject("INVALID_FRAME");
            if (i>0U) {
                const double step=detail::distance(result.poses[i-1U].position,pose.position);
                report.maximum_spacing=std::max(report.maximum_spacing,step);
                if (step<=1e-9||step>request.maximum_sample_spacing+1e-6)
                    return reject("INVALID_SAMPLE_ORDER_OR_SPACING");
            }
            if (!sdf.contains(pose.position)) return reject("OUTSIDE_SDF_DOMAIN");
            const double distance=sdf.getDistance(pose.position);
            if (!std::isfinite(distance)) return reject("NONFINITE_SDF");
            report.minimum_conservative_clearance=std::min(
                report.minimum_conservative_clearance,distance-extent-request.safe_margin);
            if (distance+1e-9<extent+request.safe_margin)
                return reject("INSUFFICIENT_SAMPLED_CLEARANCE");
        }
        std::size_t expectedStart=0U;
        for (const Segment& segment:result.segments) {
            const SegmentSpec& spec=segment.spec;
            if (segment.start_index!=expectedStart ||
                segment.end_index<=segment.start_index ||
                segment.end_index>=result.poses.size() ||
                !detail::near(result.poses[segment.start_index].position,spec.start)||
                !detail::near(result.poses[segment.end_index].position,spec.end))
                return reject("SEGMENT_COVERAGE_OR_ENDPOINT");
            Vec3 startTangent{};
            Vec3 axis{}, radial{};
            const bool twist=spec.type==SegmentType::STRAIGHT_TWIST;
            const bool straight=spec.type==SegmentType::STRAIGHT || twist;
            if (straight) {
                startTangent=detail::unit(detail::sub(spec.end,spec.start));
                if (twist) {
                    const double length=detail::distance(spec.start,spec.end);
                    const double peakFactor=spec.twist_profile==TwistProfile::SMOOTHSTEP?1.5:1.0;
                    if ((spec.twist_profile!=TwistProfile::LINEAR &&
                         spec.twist_profile!=TwistProfile::SMOOTHSTEP) ||
                        !std::isfinite(spec.twist_angle_rad)||
                        !std::isfinite(spec.twist_length)||spec.twist_length<=0.0||
                        std::abs(spec.twist_length-length)>1e-4||
                        spec.twist_length+1e-9<request.minimum_twist_length||
                        (request.maximum_twist_rate_rad_per_voxel>0.0 &&
                         peakFactor*std::abs(spec.twist_angle_rad)/spec.twist_length>
                             request.maximum_twist_rate_rad_per_voxel+1e-9)||
                        !detail::finite(spec.twist_start_normal)||
                        !detail::finite(spec.twist_end_normal)||
                        std::abs(detail::norm(spec.twist_start_normal)-1.0)>1e-5||
                        std::abs(detail::norm(spec.twist_end_normal)-1.0)>1e-5||
                        std::abs(detail::dot(spec.twist_start_normal,startTangent))>1e-5||
                        std::abs(detail::dot(spec.twist_end_normal,startTangent))>1e-5||
                        !detail::aligned(detail::rotate(spec.twist_start_normal,
                                                        startTangent,spec.twist_angle_rad),
                                         spec.twist_end_normal,1e-5))
                        return reject("INVALID_TWIST_PARAMETERS_OR_CONSTRAINT");
                }
            } else if (spec.type==SegmentType::HORIZONTAL_BEND ||
                       spec.type==SegmentType::VERTICAL_BEND) {
                if (!std::isfinite(spec.radius)||spec.radius+1e-9<request.minimum_bend_radius||
                    spec.radius<=0.0||!std::isfinite(spec.sweep_rad)||
                    std::abs(spec.sweep_rad)<1e-8)
                    return reject("INVALID_BEND_RADIUS_OR_SWEEP");
                if (!detail::finite(spec.axis) || std::abs(detail::norm(spec.axis)-1.0)>1e-5)
                    return reject("ARC_AXIS_NOT_UNIT");
                axis=spec.axis;
                radial=detail::sub(spec.start,spec.center);
                if (std::abs(detail::norm(radial)-spec.radius)>1e-4||
                    std::abs(detail::dot(radial,axis))>1e-4)
                    return reject("INVALID_ARC_PLANE_OR_RADIUS");
                startTangent=detail::unit(detail::mul(detail::cross(axis,radial),
                                                     spec.sweep_rad>0.0?1.0:-1.0));
            } else {
                return reject("UNKNOWN_SEGMENT_TYPE");
            }
            if (!detail::aligned(result.poses[segment.start_index].tangent,startTangent))
                return reject("SEGMENT_TANGENT_MISMATCH");
            const Vec3 startNormal=result.poses[segment.start_index].normal;
            if (twist && !detail::aligned(startNormal,spec.twist_start_normal,1e-5))
                return reject("TWIST_START_FRAME_MISMATCH");
            const std::size_t steps=segment.end_index-segment.start_index;
            for (std::size_t j=0U;j<=steps;++j) {
                const double fraction=static_cast<double>(j)/static_cast<double>(steps);
                const double angle=spec.sweep_rad*fraction;
                const double twistFraction=spec.twist_profile==TwistProfile::SMOOTHSTEP
                    ?fraction*fraction*(3.0-2.0*fraction):fraction;
                const Vec3 expectedPosition=straight
                    ?detail::add(spec.start,detail::mul(detail::sub(spec.end,spec.start),fraction))
                    :detail::add(spec.center,detail::rotate(radial,axis,angle));
                const Vec3 expectedTangent=straight
                    ?startTangent:detail::rotate(startTangent,axis,angle);
                const Vec3 expectedNormal=twist
                    ?detail::rotate(spec.twist_start_normal,startTangent,
                                   spec.twist_angle_rad*twistFraction)
                    :(straight?startNormal:detail::rotate(startNormal,axis,angle));
                const Pose& pose=result.poses[segment.start_index+j];
                if (!detail::near(pose.position,expectedPosition)||
                    !detail::aligned(pose.tangent,expectedTangent,1e-5)||
                    !detail::aligned(pose.normal,expectedNormal,1e-5))
                    return reject("SEGMENT_GEOMETRY_OR_FRAME_MISMATCH");
                if (twist && j>0U) {
                    const double previousFraction=static_cast<double>(j-1U)/steps;
                    const double previousTwistFraction=spec.twist_profile==TwistProfile::SMOOTHSTEP
                        ?previousFraction*previousFraction*(3.0-2.0*previousFraction)
                        :previousFraction;
                    if (std::abs(spec.twist_angle_rad)*
                            (twistFraction-previousTwistFraction)>
                        request.maximum_twist_step_rad+1e-6)
                        return reject("TWIST_SAMPLE_ANGLE_TOO_LARGE");
                }
            }
            if (twist && !detail::aligned(result.poses[segment.end_index].normal,
                                          spec.twist_end_normal,1e-5))
                return reject("TWIST_END_FRAME_MISMATCH");
            expectedStart=segment.end_index;
        }
        if (expectedStart!=result.poses.size()-1U)
            return reject("SEGMENT_COVERAGE_OR_ENDPOINT");
        report.ok=true;
        return report;
    } catch (const std::exception& error) {
        return {false,error.what(),0.0,0.0};
    }
}

} // namespace fine_planner
