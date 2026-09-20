#include "optimizer/SqpProjector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace optimizer {
namespace {

void addScaled(
    environment::Vec3& target,
    const environment::Vec3& value,
    double scale) noexcept {
    target.x += value.x * scale;
    target.y += value.y * scale;
    target.z += value.z * scale;
}

bool finiteVec(const environment::Vec3& value) noexcept {
    return std::isfinite(value.x) &&
        std::isfinite(value.y) &&
        std::isfinite(value.z);
}

double vectorNorm(const environment::Vec3& value) noexcept {
    return std::sqrt(
        value.x * value.x +
        value.y * value.y +
        value.z * value.z);
}

} // namespace

struct SqpProjector::EnergyGradient {
    double energy = 0.0;
    std::vector<environment::Vec3> gradient;
    double gradientNormSquared = 0.0;
};

struct SqpProjector::ConstraintValidation {
    bool ok = false;
    double minimumClearance = -std::numeric_limits<double>::infinity();
    double maximumCurvature = 0.0;
};

SqpProjector::SqpProjector(SqpProjectorOptions options)
    : options_(options) {
    options_.max_iterations = std::max(1, options_.max_iterations);
    options_.safe_margin = std::max(0.0, options_.safe_margin);
    options_.busbar_width = std::max(0.0, options_.busbar_width);
    options_.busbar_thickness = std::max(0.0, options_.busbar_thickness);
    if (std::isnan(options_.max_curvature) ||
        options_.max_curvature <= 0.0) {
        options_.max_curvature =
            std::numeric_limits<double>::infinity();
    }
    options_.validation_sample_voxel_fraction = std::clamp(
        options_.validation_sample_voxel_fraction,
        1.0e-4,
        0.25);
    options_.tension_weight = std::max(0.0, options_.tension_weight);
    options_.collision_weight = std::max(0.0, options_.collision_weight);
    options_.smoothness_weight = std::max(0.0, options_.smoothness_weight);
    options_.gradient_tolerance =
        std::max(0.0, options_.gradient_tolerance);
    options_.initial_step = std::max(1.0e-9, options_.initial_step);
    options_.minimum_step = std::max(1.0e-12, options_.minimum_step);
    options_.armijo_factor = std::clamp(
        options_.armijo_factor,
        1.0e-8,
        0.5);
    options_.shrink_factor = std::clamp(
        options_.shrink_factor,
        0.05,
        0.95);
}

SqpProjectorResult SqpProjector::project(
    const PoseTrajectory& initial,
    const environment::GlobalSdf& sdf) const {
    SqpProjectorResult result;
    result.trajectory = initial;
    result.effective_safe_margin = effectiveSafeMargin();
    if (initial.size() < 2U) {
        return result;
    }

    if (!insideSdfDomain(initial, sdf)) {
        return result;
    }

    if (initial.size() == 2U) {
        result.trajectory.recomputeFrames();
        const ConstraintValidation validation =
            validateConstraints(result.trajectory, sdf);
        result.constraints_satisfied = validation.ok;
        result.minimum_clearance = validation.minimumClearance;
        result.maximum_curvature = validation.maximumCurvature;
        result.success = validation.ok;
        return result;
    }

    const double targetSpacing =
        initial.polylineLength() /
        static_cast<double>(std::max<std::size_t>(1U, initial.size() - 1U));
    if (!std::isfinite(targetSpacing) || targetSpacing <= 1.0e-9) {
        return result;
    }

    PoseTrajectory current = initial;
    EnergyGradient state = evaluate(
        current,
        sdf,
        targetSpacing,
        true);
    double previousAcceptedStep = options_.initial_step;

    for (int iteration = 0;
         iteration < options_.max_iterations;
         ++iteration) {
        result.iterations = iteration + 1;
        result.final_energy = state.energy;
        result.gradient_norm = std::sqrt(state.gradientNormSquared);
        if (result.gradient_norm <= options_.gradient_tolerance) {
            break;
        }

        double step = previousAcceptedStep;
        bool accepted = false;
        while (step >= options_.minimum_step) {
            PoseTrajectory candidate = stepTrajectory(
                current,
                state.gradient,
                step);
            if (!insideSdfDomain(candidate, sdf)) {
                step *= options_.shrink_factor;
                continue;
            }
            const EnergyGradient candidateState = evaluate(
                candidate,
                sdf,
                targetSpacing,
                false);
            if (std::isfinite(candidateState.energy) &&
                candidateState.energy <=
                    state.energy -
                        options_.armijo_factor *
                        step *
                        state.gradientNormSquared) {
                current = std::move(candidate);
                state = evaluate(
                    current,
                    sdf,
                    targetSpacing,
                    true);
                previousAcceptedStep = std::min(
                    options_.initial_step,
                    step / options_.shrink_factor);
                accepted = true;
                break;
            }
            step *= options_.shrink_factor;
        }

        if (!accepted) {
            break;
        }
    }

    current.recomputeFrames();
    result.final_energy = state.energy;
    result.gradient_norm = std::sqrt(state.gradientNormSquared);
    const ConstraintValidation validation =
        validateConstraints(current, sdf);
    result.constraints_satisfied = validation.ok;
    result.minimum_clearance = validation.minimumClearance;
    result.maximum_curvature = validation.maximumCurvature;

    // A trajectory that only lowered the penalty energy is not a valid
    // projection. Keep the original A* trajectory so callers can safely
    // fall back without accidentally publishing a failed continuous result.
    if (!validation.ok || !std::isfinite(result.final_energy)) {
        result.trajectory = initial;
        result.success = false;
        return result;
    }

    result.trajectory = std::move(current);
    result.success = true;
    return result;
}

bool SqpProjector::satisfiesHardConstraints(
    const PoseTrajectory& trajectory,
    const environment::GlobalSdf& sdf) const {
    return validateConstraints(trajectory, sdf).ok;
}

SqpProjector::EnergyGradient SqpProjector::evaluate(
    const PoseTrajectory& trajectory,
    const environment::GlobalSdf& sdf,
    double targetSpacing,
    bool withGradient) const {
    EnergyGradient result;
    if (withGradient) {
        result.gradient.assign(trajectory.size(), environment::Vec3{});
    }

    const auto mayWriteGradient = [this, &trajectory](
                                      std::size_t index) {
        return !options_.lock_endpoints ||
            (index != 0U && index + 1U < trajectory.size());
    };

    for (std::size_t index = 1U;
         index < trajectory.size();
         ++index) {
        const environment::Vec3 edge = PoseTrajectory::subtract(
            trajectory[index].position,
            trajectory[index - 1U].position);
        const double length = PoseTrajectory::norm(edge);
        if (length <= 1.0e-9 || !std::isfinite(length)) {
            continue;
        }
        const double residual = length - targetSpacing;
        result.energy += options_.tension_weight * residual * residual;
        if (withGradient && options_.tension_weight > 0.0) {
            const environment::Vec3 direction =
                PoseTrajectory::scale(edge, 1.0 / length);
            const environment::Vec3 contribution =
                PoseTrajectory::scale(
                    direction,
                    2.0 * options_.tension_weight * residual);
            if (mayWriteGradient(index)) {
                addScaled(result.gradient[index], contribution, 1.0);
            }
            if (mayWriteGradient(index - 1U)) {
                addScaled(result.gradient[index - 1U], contribution, -1.0);
            }
        }
    }

    for (std::size_t index = 1U;
         index + 1U < trajectory.size();
         ++index) {
        const environment::Vec3 secondDifference =
            PoseTrajectory::add(
                PoseTrajectory::subtract(
                    trajectory[index - 1U].position,
                    PoseTrajectory::scale(
                        trajectory[index].position,
                        2.0)),
                trajectory[index + 1U].position);
        result.energy += options_.smoothness_weight *
            PoseTrajectory::squaredNorm(secondDifference);
        if (withGradient && options_.smoothness_weight > 0.0) {
            const environment::Vec3 contribution =
                PoseTrajectory::scale(
                    secondDifference,
                    2.0 * options_.smoothness_weight);
            if (mayWriteGradient(index - 1U)) {
                addScaled(result.gradient[index - 1U], contribution, 1.0);
            }
            if (mayWriteGradient(index)) {
                addScaled(result.gradient[index], contribution, -2.0);
            }
            if (mayWriteGradient(index + 1U)) {
                addScaled(result.gradient[index + 1U], contribution, 1.0);
            }
        }
    }

    for (std::size_t index = 0U;
         index < trajectory.size();
         ++index) {
        const double distance =
            sdf.getDistance(trajectory[index].position);
        const double requiredDistance = effectiveSafeMargin();
        if (!std::isfinite(distance) ||
            distance >= requiredDistance) {
            continue;
        }
        const double residual = requiredDistance - distance;
        result.energy += options_.collision_weight * residual * residual;
        if (withGradient && options_.collision_weight > 0.0 &&
            mayWriteGradient(index)) {
            const environment::Vec3 sdfGradient =
                sdf.getGradient(trajectory[index].position);
            if (finiteVec(sdfGradient)) {
                addScaled(
                    result.gradient[index],
                    sdfGradient,
                    -2.0 * options_.collision_weight * residual);
            }
        }
    }

    if (withGradient) {
        for (const environment::Vec3& value : result.gradient) {
            result.gradientNormSquared +=
                PoseTrajectory::squaredNorm(value);
        }
    }
    return result;
}

bool SqpProjector::insideSdfDomain(
    const PoseTrajectory& trajectory,
    const environment::GlobalSdf& sdf) const noexcept {
    for (const TrajectoryNode& node : trajectory.nodes()) {
        if (!pointInsideSdfDomain(node.position, sdf)) {
            return false;
        }
    }
    return true;
}

bool SqpProjector::pointInsideSdfDomain(
    const environment::Vec3& position,
    const environment::GlobalSdf& sdf) const noexcept {
    const environment::Vec3& origin = sdf.origin();
    const environment::Vec3& voxelSize = sdf.voxelSize();
    const std::array<std::size_t, 3>& dimensions = sdf.dimensions();
    const environment::Vec3 maximum{
        origin.x + voxelSize.x * static_cast<double>(dimensions[0] - 1U),
        origin.y + voxelSize.y * static_cast<double>(dimensions[1] - 1U),
        origin.z + voxelSize.z * static_cast<double>(dimensions[2] - 1U)};
    return finiteVec(position) &&
        position.x >= origin.x &&
        position.y >= origin.y &&
        position.z >= origin.z &&
        position.x <= maximum.x &&
        position.y <= maximum.y &&
        position.z <= maximum.z;
}

double SqpProjector::effectiveSafeMargin() const noexcept {
    if (!std::isfinite(options_.busbar_width) ||
        !std::isfinite(options_.busbar_thickness) ||
        options_.busbar_width <= 0.0 ||
        options_.busbar_thickness <= 0.0) {
        return options_.safe_margin;
    }

    const double halfDiagonal = 0.5 * std::hypot(
        options_.busbar_width,
        options_.busbar_thickness);
    if (!std::isfinite(halfDiagonal)) {
        return std::numeric_limits<double>::infinity();
    }
    return options_.safe_margin + halfDiagonal;
}

SqpProjector::ConstraintValidation SqpProjector::validateConstraints(
    const PoseTrajectory& trajectory,
    const environment::GlobalSdf& sdf) const {
    ConstraintValidation result;
    if (trajectory.size() < 2U ||
        !insideSdfDomain(trajectory, sdf)) {
        return result;
    }

    const double requiredDistance = effectiveSafeMargin();
    if (!std::isfinite(requiredDistance)) {
        return result;
    }

    const environment::Vec3& voxelSize = sdf.voxelSize();
    const double minimumVoxelSize = std::min({
        voxelSize.x,
        voxelSize.y,
        voxelSize.z});
    const double sampleStep = std::max(
        1.0e-9,
        options_.validation_sample_voxel_fraction * minimumVoxelSize);

    result.minimumClearance = std::numeric_limits<double>::infinity();
    for (std::size_t index = 1U;
         index < trajectory.size();
         ++index) {
        const environment::Vec3 start =
            trajectory[index - 1U].position;
        const environment::Vec3 end = trajectory[index].position;
        const environment::Vec3 edge =
            PoseTrajectory::subtract(end, start);
        const double length = vectorNorm(edge);
        if (!std::isfinite(length) || length <= 1.0e-9) {
            return result;
        }

        const std::size_t sampleCount = std::max<std::size_t>(
            1U,
            static_cast<std::size_t>(
                std::ceil(length / sampleStep)));
        for (std::size_t sample = 0U;
             sample <= sampleCount;
             ++sample) {
            const double fraction =
                static_cast<double>(sample) /
                static_cast<double>(sampleCount);
            const environment::Vec3 point =
                PoseTrajectory::add(
                    start,
                    PoseTrajectory::scale(edge, fraction));
            if (!pointInsideSdfDomain(point, sdf)) {
                return result;
            }
            const double distance = sdf.getDistance(point);
            if (!std::isfinite(distance)) {
                return result;
            }
            result.minimumClearance = std::min(
                result.minimumClearance,
                distance - requiredDistance);
            if (distance < requiredDistance) {
                return result;
            }
        }
    }

    for (std::size_t index = 1U;
         index + 1U < trajectory.size();
         ++index) {
        const environment::Vec3 firstEdge =
            PoseTrajectory::subtract(
                trajectory[index].position,
                trajectory[index - 1U].position);
        const environment::Vec3 secondEdge =
            PoseTrajectory::subtract(
                trajectory[index + 1U].position,
                trajectory[index].position);
        const environment::Vec3 chord =
            PoseTrajectory::add(firstEdge, secondEdge);
        const double firstLength = vectorNorm(firstEdge);
        const double secondLength = vectorNorm(secondEdge);
        const double chordLength = vectorNorm(chord);
        if (!std::isfinite(firstLength) ||
            !std::isfinite(secondLength) ||
            !std::isfinite(chordLength) ||
            firstLength <= 1.0e-9 ||
            secondLength <= 1.0e-9 ||
            chordLength <= 1.0e-9) {
            return result;
        }

        const environment::Vec3 crossProduct =
            PoseTrajectory::cross(firstEdge, secondEdge);
        const double curvature =
            2.0 * vectorNorm(crossProduct) /
            (firstLength * secondLength * chordLength);
        if (!std::isfinite(curvature)) {
            return result;
        }
        result.maximumCurvature = std::max(
            result.maximumCurvature,
            curvature);
        if (std::isfinite(options_.max_curvature) &&
            curvature > options_.max_curvature + 1.0e-9) {
            return result;
        }
    }

    result.ok = std::isfinite(result.minimumClearance) &&
        result.minimumClearance >= -1.0e-9;
    return result;
}

PoseTrajectory SqpProjector::stepTrajectory(
    const PoseTrajectory& current,
    const std::vector<environment::Vec3>& gradient,
    double step) const {
    PoseTrajectory candidate = current;
    for (std::size_t index = 0U;
         index < candidate.size();
         ++index) {
        if (options_.lock_endpoints &&
            (index == 0U || index + 1U == candidate.size())) {
            continue;
        }
        candidate[index].position = PoseTrajectory::subtract(
            candidate[index].position,
            PoseTrajectory::scale(gradient[index], step));
    }
    return candidate;
}

} // namespace optimizer
