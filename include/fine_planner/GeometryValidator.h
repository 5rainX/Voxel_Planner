#pragma once

#include "fine_planner/FinePlanner.h"

namespace fine_planner {

// Independent of FinePlanner::plan: validates an already constructed result.
class GeometryValidator {
public:
    static ValidationReport validate(
        const Request& request,
        const Result& result,
        const environment::GlobalSdf& sdf);
};

} // namespace fine_planner
