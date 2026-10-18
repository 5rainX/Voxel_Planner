#pragma once

#include "fine_planner/FinePlanner.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fine_planner::detail {

inline Vec3 add(Vec3 a, Vec3 b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline Vec3 sub(Vec3 a, Vec3 b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline Vec3 mul(Vec3 a, double s) { return {a.x*s, a.y*s, a.z*s}; }
inline double dot(Vec3 a, Vec3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
inline double norm(Vec3 a) { return std::sqrt(dot(a,a)); }
inline bool finite(Vec3 a) { return std::isfinite(a.x)&&std::isfinite(a.y)&&std::isfinite(a.z); }
inline Vec3 unit(Vec3 a) {
    const double n=norm(a);
    if (!finite(a)||!std::isfinite(n)||n<1e-9) throw std::invalid_argument("DEGENERATE_VECTOR");
    return mul(a,1.0/n);
}
inline Vec3 rotate(Vec3 v, Vec3 axis, double angle) {
    return add(add(mul(v,std::cos(angle)),mul(cross(axis,v),std::sin(angle))),
               mul(axis,dot(axis,v)*(1-std::cos(angle))));
}
inline double distance(Vec3 a, Vec3 b) { return norm(sub(a,b)); }
inline bool near(Vec3 a, Vec3 b, double tolerance=1e-4) { return distance(a,b)<=tolerance; }
inline bool aligned(Vec3 a, Vec3 b, double tolerance=1e-4) {
    return dot(unit(a),unit(b)) >= 1.0-tolerance;
}
inline Pose makePose(Vec3 position, Vec3 tangent, Vec3 normal) {
    tangent=unit(tangent);
    normal=unit(sub(normal,mul(tangent,dot(normal,tangent))));
    return {position,tangent,normal,unit(cross(tangent,normal))};
}
inline const char* typeName(SegmentType type) {
    switch(type) {
    case SegmentType::STRAIGHT: return "STRAIGHT";
    case SegmentType::HORIZONTAL_BEND: return "HORIZONTAL_BEND";
    case SegmentType::VERTICAL_BEND: return "VERTICAL_BEND";
    }
    return "UNKNOWN";
}

} // namespace fine_planner::detail
