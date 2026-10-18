#include "fine_planner/FinePlanner.h"
#include "fine_planner/GeometryValidator.h"
#include "fine_planner/ReferenceExport.h"
#include "Module1_Main/Types.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

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
        if (argc>1 && !vertical) {
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
    if (argc > 1) {
        require(argc>2,"map hash required for export fixture");
    }
    return 0;
}
