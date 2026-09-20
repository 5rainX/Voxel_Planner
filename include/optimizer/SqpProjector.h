#pragma once

#include "environment/SdfVolume.h"
#include "optimizer/Trajectory.h"

#include <cstddef>
#include <limits>

namespace optimizer {

struct SqpProjectorOptions {
    int max_iterations = 80;
    double safe_margin = 1.0;
    double busbar_width = 0.0;
    double busbar_thickness = 0.0;
    double max_curvature = 1.0;
    double validation_sample_voxel_fraction = 0.25;
    double tension_weight = 0.15;
    double collision_weight = 18.0;
    double smoothness_weight = 0.65;
    double gradient_tolerance = 1.0e-4;
    double initial_step = 0.12;
    double minimum_step = 1.0e-5;
    double armijo_factor = 1.0e-4;
    double shrink_factor = 0.5;
    bool lock_endpoints = true;
};

struct SqpProjectorResult {
    bool success = false;
    PoseTrajectory trajectory;
    double final_energy = 0.0;
    double gradient_norm = 0.0;
    bool constraints_satisfied = false;
    double minimum_clearance = -std::numeric_limits<double>::infinity();
    double maximum_curvature = 0.0;
    double effective_safe_margin = 0.0;
    int iterations = 0;
};

class SqpProjector {
public:
    explicit SqpProjector(SqpProjectorOptions options = {});

    SqpProjectorResult project(
        const PoseTrajectory& initial,
        const environment::GlobalSdf& sdf) const;

    bool satisfiesHardConstraints(
        const PoseTrajectory& trajectory,
        const environment::GlobalSdf& sdf) const;

private:
    struct EnergyGradient;
    struct ConstraintValidation;

    EnergyGradient evaluate(
        const PoseTrajectory& trajectory,
        const environment::GlobalSdf& sdf,
        double targetSpacing,
        bool withGradient) const;

    ConstraintValidation validateConstraints(
        const PoseTrajectory& trajectory,
        const environment::GlobalSdf& sdf) const;

    bool insideSdfDomain(
        const PoseTrajectory& trajectory,
        const environment::GlobalSdf& sdf) const noexcept;

    bool pointInsideSdfDomain(
        const environment::Vec3& position,
        const environment::GlobalSdf& sdf) const noexcept;

    PoseTrajectory stepTrajectory(
        const PoseTrajectory& current,
        const std::vector<environment::Vec3>& gradient,
        double step) const;

    double effectiveSafeMargin() const noexcept;

    SqpProjectorOptions options_{};
};

} // namespace optimizer
