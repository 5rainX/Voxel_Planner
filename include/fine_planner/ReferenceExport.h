#pragma once

#include "fine_planner/FinePlanner.h"

#include <filesystem>

namespace fine_planner {

// Writes the existing Reference Geometry V1 object, not a legacy TXT sidecar.
// Caller supplies IDs and an independently validated result.
void writeReferenceGeometryV1(
    const Result& result,
    const Request& request,
    const ValidationReport& validation,
    const std::filesystem::path& output);

} // namespace fine_planner
