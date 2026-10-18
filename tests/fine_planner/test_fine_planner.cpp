#include "fine_planner/FinePlanner.h"
#include "fine_planner/GeometryValidator.h"
#include "fine_planner/ReferenceExport.h"
#include "Module1_Main/Types.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {
void require(bool ok, const char* message) {
    if (!ok) { std::cerr << message << '\n'; std::exit(1); }
}

fine_planner::Request straightRequest() {
    fine_planner::Request request;
    request.coarse_path = {{2, 2, 2}, {6, 2, 2}};
    request.start_tangent = {1, 0, 0};
    request.end_tangent = {1, 0, 0};
    request.start_normal = {0, 0, 1};
    request.end_normal = {0, 0, 1};
    request.width = 0.5;
    request.thickness = 0.5;
    request.sample_id = "mvp:straight";
    request.scene_id = "mvp";
    request.map_sha256 = std::string(64, 'a');
    return request;
}

fine_planner::Request bendRequest(bool vertical) {
    auto request=straightRequest();
    request.sample_id=vertical?"mvp:vertical":"mvp:horizontal";
    request.coarse_path=vertical
        ?std::vector<fine_planner::Vec3>{{2,4,4},{4,4,4},{6,4,6}}
        :std::vector<fine_planner::Vec3>{{2,4,4},{4,4,4},{6,6,4}};
    request.segments.push_back({fine_planner::SegmentType::STRAIGHT,
                                {2,4,4},{4,4,4},{},{},0,0});
    request.segments.push_back(vertical
        ?fine_planner::SegmentSpec{fine_planner::SegmentType::VERTICAL_BEND,
                                  {4,4,4},{6,4,6},{4,4,6},{0,-1,0},2,1.5707963267948966}
        :fine_planner::SegmentSpec{fine_planner::SegmentType::HORIZONTAL_BEND,
                                  {4,4,4},{6,6,4},{4,6,4},{0,0,1},2,1.5707963267948966});
    request.start_normal=vertical?fine_planner::Vec3{0,1,0}:fine_planner::Vec3{0,0,1};
    request.end_normal=request.start_normal;
    request.end_tangent=vertical?fine_planner::Vec3{0,0,1}:fine_planner::Vec3{0,1,0};
    request.minimum_bend_radius=1.5;
    return request;
}

fine_planner::Request twistRequest(double angle,
                                   fine_planner::TwistProfile profile = fine_planner::TwistProfile::LINEAR) {
    auto request=straightRequest();
    request.sample_id="mvp:twist";
    request.minimum_twist_length=2.0;
    request.maximum_twist_rate_rad_per_voxel=0.7;
    const fine_planner::Vec3 endNormal{0,-std::sin(angle),std::cos(angle)};
    request.end_normal=endNormal;
    fine_planner::SegmentSpec twist;
    twist.type=fine_planner::SegmentType::STRAIGHT_TWIST;
    twist.start={2,2,2};
    twist.end={6,2,2};
    twist.twist_start_normal={0,0,1};
    twist.twist_end_normal=endNormal;
    twist.twist_angle_rad=angle;
    twist.twist_length=4.0;
    twist.twist_profile=profile;
    request.segments={twist};
    return request;
}

fine_planner::Request crossValidationRequest(const std::string& name, double spacing) {
    constexpr double pi=3.14159265358979323846;
    double angle=0.0;
    if (name=="positive90") angle=pi/2;
    else if (name=="negative90") angle=-pi/2;
    else if (name=="positive180") angle=pi;
    else if (name=="positive360") angle=2*pi;
    else if (name=="nonaxis") angle=1.2;
    else if (name=="arbitraryframe") angle=-0.7;
    else if (name=="short") angle=0.1;
    else if (name=="smoothstep") angle=1.1;
    else require(name=="zero","unknown twist cross-validation case");
    auto request=twistRequest(angle);
    request.maximum_twist_rate_rad_per_voxel=0.0;
    request.minimum_twist_length=0.2;
    request.maximum_twist_step_rad=0.5;
    request.maximum_sample_spacing=spacing;
    fine_planner::Vec3 start={2,2,2}, end={6,2,2};
    fine_planner::Vec3 tangent={1,0,0}, normal={0,0,1};
    if (name=="nonaxis") {
        constexpr double invRoot3=0.57735026918962576451;
        constexpr double invRoot2=0.70710678118654752440;
        start={3,3,3}; end={7,7,7};
        tangent={invRoot3,invRoot3,invRoot3};
        normal={invRoot2,-invRoot2,0};
        request.width=0.4; request.thickness=0.9;
    } else if (name=="arbitraryframe") {
        constexpr double invRoot2=0.70710678118654752440;
        normal={0,invRoot2,invRoot2};
        request.width=0.8; request.thickness=0.3;
    } else if (name=="short") {
        start={3,3,3}; end={3.25,3,3};
        request.width=0.2; request.thickness=0.4;
    } else if (name=="smoothstep") {
        request.width=0.6; request.thickness=0.2;
        request.segments[0].twist_profile=fine_planner::TwistProfile::SMOOTHSTEP;
    }
    const fine_planner::Vec3 binormal={
        tangent.y*normal.z-tangent.z*normal.y,
        tangent.z*normal.x-tangent.x*normal.z,
        tangent.x*normal.y-tangent.y*normal.x};
    const fine_planner::Vec3 endNormal={
        std::cos(angle)*normal.x+std::sin(angle)*binormal.x,
        std::cos(angle)*normal.y+std::sin(angle)*binormal.y,
        std::cos(angle)*normal.z+std::sin(angle)*binormal.z};
    request.sample_id="cross-validation:"+name;
    request.coarse_path={start,end};
    request.start_tangent=tangent;
    request.end_tangent=tangent;
    request.start_normal=normal;
    request.end_normal=endNormal;
    auto& segment=request.segments[0];
    segment.start=start;
    segment.end=end;
    segment.twist_start_normal=normal;
    segment.twist_end_normal=endNormal;
    segment.twist_angle_rad=angle;
    segment.twist_length=std::sqrt((end.x-start.x)*(end.x-start.x)+
                                   (end.y-start.y)*(end.y-start.y)+
                                   (end.z-start.z)*(end.z-start.z));
    return request;
}
}

int main(int argc, char** argv) {
    VoxelGrid grid(16U, 16U, 16U);
    const environment::GlobalSdf sdf(grid);
    const auto request = straightRequest();
    const auto result = fine_planner::FinePlanner::plan(request, sdf);
    require(result.success, "straight plan failed");
    require(result.poses.size() == 17U, "straight spacing/point count");
    require(result.segments.size() == 1U, "straight segment count");
    require(result.segments[0].spec.type == fine_planner::SegmentType::STRAIGHT, "straight type");
    require(std::abs(result.poses.front().position.x - 2.0) < 1e-9, "start position");
    require(std::abs(result.poses.back().position.x - 6.0) < 1e-9, "end position");
    const auto validation = fine_planner::GeometryValidator::validate(request, result, sdf);
    require(validation.ok, "independent straight validation failed");

    for (bool vertical:{false,true}) {
        const auto bend=bendRequest(vertical);
        const auto curved=fine_planner::FinePlanner::plan(bend,sdf);
        require(curved.success,vertical?"vertical bend failed":"horizontal bend failed");
        require(curved.segments.size()==2U,"bend segment count");
        const auto& arc=curved.segments[1];
        require(arc.start_index==curved.segments[0].end_index,"duplicate/missing junction");
        for (std::size_t i=arc.start_index;i<=arc.end_index;++i) {
            const auto& pose=curved.poses[i];
            const auto delta=fine_planner::Vec3{pose.position.x-arc.spec.center.x,
                                                pose.position.y-arc.spec.center.y,
                                                pose.position.z-arc.spec.center.z};
            require(std::abs(std::sqrt(delta.x*delta.x+delta.y*delta.y+delta.z*delta.z)-2.0)<1e-6,
                    "arc radius");
            require(std::abs(pose.tangent.x*pose.normal.x+pose.tangent.y*pose.normal.y+
                             pose.tangent.z*pose.normal.z)<1e-6,"frame orthogonality");
        }
        require(fine_planner::GeometryValidator::validate(bend,curved,sdf).ok,
                "independent bend validation failed");
        auto changedRequest=bend;
        changedRequest.segments[1].sweep_rad*=-1;
        require(!fine_planner::GeometryValidator::validate(changedRequest,curved,sdf).ok,
                "result was not bound to requested arc parameters");
        auto broken=curved;
        broken.poses[arc.start_index+1U].normal={1,0,0};
        require(!fine_planner::GeometryValidator::validate(bend,broken,sdf).ok,
                "tampered frame accepted");
        auto tooTight=bend;
        tooTight.minimum_bend_radius=3.0;
        require(!fine_planner::FinePlanner::plan(tooTight,sdf).success,
                "minimum bend radius ignored");
        auto nonUnitAxis=bend;
        nonUnitAxis.segments[1].axis.x*=2.0;
        nonUnitAxis.segments[1].axis.y*=2.0;
        nonUnitAxis.segments[1].axis.z*=2.0;
        require(!fine_planner::FinePlanner::plan(nonUnitAxis,sdf).success,
                "nonunit arc axis accepted");
        auto badEnd=bend;
        badEnd.end_tangent={-1,0,0};
        require(!fine_planner::FinePlanner::plan(badEnd,sdf).success,
                "endpoint tangent ignored");
        auto badJunction=bend;
        badJunction.segments[1].start={4,4.5,4};
        require(!fine_planner::FinePlanner::plan(badJunction,sdf).success,
                "discontinuous segment junction accepted");
        if (argc>1 && !vertical &&
            (argc<4 || (std::string(argv[3])!="twist" &&
                         std::string(argv[3])!="twist-case"))) {
            auto exported=curved;
            auto exportRequest=bend;
            if (argc>2) exportRequest.map_sha256=argv[2];
            fine_planner::writeReferenceGeometryV1(exported,exportRequest,
                fine_planner::GeometryValidator::validate(bend,curved,sdf),argv[1]);
        }
    }
    auto unknown=straightRequest();
    unknown.coarse_path={{2,2,2},{4,3,2},{6,2,2}};
    require(!fine_planner::FinePlanner::plan(unknown,sdf).success,
            "coarse zigzag was guessed into a curve");
    VoxelGrid blocked(16U,16U,16U);
    blocked.setRawObstacle(4,2,2);
    const environment::GlobalSdf blockedSdf(blocked);
    require(!fine_planner::FinePlanner::plan(request,blockedSdf).success,
            "occupied straight accepted");
    VoxelGrid adjacent(16U,16U,16U);
    adjacent.setRawObstacle(4,3,2);
    const environment::GlobalSdf adjacentSdf(adjacent);
    auto wide=straightRequest();
    wide.width=2.0;
    require(!fine_planner::FinePlanner::plan(wide,adjacentSdf).success,
            "section extent did not reject adjacent obstacle");
    constexpr double halfPi=1.5707963267948966;
    for (double angle:{0.0,halfPi,-halfPi,0.7}) {
        for (auto profile:{fine_planner::TwistProfile::LINEAR,
                           fine_planner::TwistProfile::SMOOTHSTEP}) {
            const auto twist=twistRequest(angle,profile);
            const auto twisted=fine_planner::FinePlanner::plan(twist,sdf);
            require(twisted.success,"straight twist failed");
            require(twisted.poses.size()>=17U,"twist sampling too sparse");
            const auto& middle=twisted.poses[twisted.poses.size()/2U];
            require(std::abs(middle.position.y-2.0)<1e-9 &&
                    std::abs(middle.position.z-2.0)<1e-9,"twist centerline bent");
            require(std::abs(middle.normal.y+std::sin(angle/2.0))<1e-6 &&
                    std::abs(middle.normal.z-std::cos(angle/2.0))<1e-6,
                    "twist midpoint orientation incorrect");
            const auto& quarter=twisted.poses[(twisted.poses.size()-1U)/4U];
            const double quarterFraction=static_cast<double>((twisted.poses.size()-1U)/4U)/
                                         static_cast<double>(twisted.poses.size()-1U);
            const double profileFraction=profile==fine_planner::TwistProfile::SMOOTHSTEP
                ?quarterFraction*quarterFraction*(3.0-2.0*quarterFraction)
                :quarterFraction;
            require(std::abs(quarter.normal.y+std::sin(angle*profileFraction))<1e-6 &&
                    std::abs(quarter.normal.z-std::cos(angle*profileFraction))<1e-6,
                    "twist profile not applied to section frame");
            require(fine_planner::GeometryValidator::validate(twist,twisted,sdf).ok,
                    "independent twist validation failed");
            auto tampered=twisted;
            tampered.poses[tampered.poses.size()/2U].normal={0,0,1};
            if (std::abs(angle)>0.01)
                require(!fine_planner::GeometryValidator::validate(twist,tampered,sdf).ok,
                        "validator accepted lost twist orientation");
        }
    }
    auto tooShort=twistRequest(halfPi);
    tooShort.minimum_twist_length=5.0;
    require(!fine_planner::FinePlanner::plan(tooShort,sdf).success,"minimum twist length ignored");
    auto tooFast=twistRequest(halfPi);
    tooFast.maximum_twist_rate_rad_per_voxel=0.2;
    require(!fine_planner::FinePlanner::plan(tooFast,sdf).success,"maximum twist rate ignored");
    auto profileRate=twistRequest(halfPi);
    profileRate.maximum_twist_rate_rad_per_voxel=0.5;
    require(fine_planner::FinePlanner::plan(profileRate,sdf).success,
            "linear twist incorrectly rejected by rate bound");
    profileRate.segments[0].twist_profile=fine_planner::TwistProfile::SMOOTHSTEP;
    require(!fine_planner::FinePlanner::plan(profileRate,sdf).success,
            "smoothstep peak twist rate ignored");
    auto badEndFrame=twistRequest(halfPi);
    badEndFrame.segments[0].twist_end_normal={0,0,1};
    require(!fine_planner::FinePlanner::plan(badEndFrame,sdf).success,"twist end frame ignored");
    auto badLength=twistRequest(halfPi);
    badLength.segments[0].twist_length=3.0;
    require(!fine_planner::FinePlanner::plan(badLength,sdf).success,"twist length/chord mismatch accepted");
    const auto validTwistRequest=twistRequest(halfPi);
    const auto validTwist=fine_planner::FinePlanner::plan(validTwistRequest,sdf);
    require(validTwist.success,"validator fixture failed");
    auto nonfiniteRequest=validTwistRequest;
    nonfiniteRequest.segments[0].twist_angle_rad=std::numeric_limits<double>::quiet_NaN();
    require(!fine_planner::GeometryValidator::validate(nonfiniteRequest,validTwist,sdf).ok,
            "validator accepted nonfinite requested twist angle");
    nonfiniteRequest=validTwistRequest;
    nonfiniteRequest.segments[0].twist_length=std::numeric_limits<double>::quiet_NaN();
    require(!fine_planner::GeometryValidator::validate(nonfiniteRequest,validTwist,sdf).ok,
            "validator accepted nonfinite requested twist length");
    auto badDimension=twistRequest(halfPi);
    badDimension.width=0.0;
    require(!fine_planner::FinePlanner::plan(badDimension,sdf).success,"invalid width accepted");
    auto blockedTwist=twistRequest(halfPi);
    require(!fine_planner::FinePlanner::plan(blockedTwist,blockedSdf).success,
            "twist passed occupied SDF sample");
    auto wideTwist=twistRequest(halfPi);
    wideTwist.width=2.0;
    require(!fine_planner::FinePlanner::plan(wideTwist,adjacentSdf).success,
            "twist conservative section extent ignored adjacent obstacle");
    auto zeroLength=twistRequest(0.0);
    zeroLength.segments[0].end=zeroLength.segments[0].start;
    zeroLength.coarse_path.back()=zeroLength.coarse_path.front();
    zeroLength.segments[0].twist_length=0.0;
    require(!fine_planner::FinePlanner::plan(zeroLength,sdf).success,
            "zero-length twist accepted");
    for (double angle:{3.14159265358979323846,6.28318530717958647692}) {
        auto fullTurn=twistRequest(angle);
        fullTurn.maximum_twist_rate_rad_per_voxel=0.0;
        const auto turned=fine_planner::FinePlanner::plan(fullTurn,sdf);
        require(turned.success,"large signed twist rejected");
        const auto quarterIndex=(turned.poses.size()-1U)/4U;
        const auto& quarter=turned.poses[quarterIndex];
        const double quarterAngle=angle*static_cast<double>(quarterIndex)/
                                  static_cast<double>(turned.poses.size()-1U);
        require(std::abs(quarter.normal.y+std::sin(quarterAngle))<1e-6 &&
                std::abs(quarter.normal.z-std::cos(quarterAngle))<1e-6,
                "large twist collapsed to endpoint frame");
    }
    if (argc>3 && std::string(argv[3])=="twist-case") {
        require(argc==6,"twist cross-validation arguments");
        auto exportRequest=crossValidationRequest(argv[4],std::stod(argv[5]));
        exportRequest.map_sha256=argv[2];
        const auto exportResult=fine_planner::FinePlanner::plan(exportRequest,sdf);
        require(exportResult.success,"cross-validation twist fixture failed");
        fine_planner::writeReferenceGeometryV1(exportResult,exportRequest,
                                               exportResult.validation,argv[1]);
    }
    if (argc>3 && std::string(argv[3])=="twist") {
        auto exportRequest=twistRequest(halfPi);
        exportRequest.map_sha256=argv[2];
        const auto exportResult=fine_planner::FinePlanner::plan(exportRequest,sdf);
        require(exportResult.success,"twist fixture failed");
        fine_planner::writeReferenceGeometryV1(exportResult,exportRequest,
                                               exportResult.validation,argv[1]);
    }
    if (argc > 1) {
        require(argc>2,"map hash required for export fixture");
    }
    return 0;
}
