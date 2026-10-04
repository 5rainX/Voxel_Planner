"""Dense geometry sampling from explicit, provenance-backed segment geometry."""

import math

import numpy as np

from path_transformer.dataset.geometry import SampleError, unit, vec


SPACING_VOXEL = 0.25  # matches SqpProjectorOptions.validation_sample_voxel_fraction


def sample_linear(start, end, spacing=SPACING_VOXEL):
    start, end = vec(start), vec(end)
    length = float(np.linalg.norm(end - start))
    if length < 1e-8 or spacing <= 0:
        raise SampleError("REFERENCE_GEOMETRY_INVALID")
    steps = max(1, math.ceil(length / spacing))
    return np.array([start + (end - start) * i / steps for i in range(steps + 1)])


def sample_explicit_arc(start, center, axis, radius, sweep, spacing=SPACING_VOXEL):
    """Only for a source that explicitly supplies all arc parameters.

    This helper is not invoked for historical VERTICAL_BEND summaries.
    """
    start, center, axis = vec(start), vec(center), unit(vec(axis))
    radial = start - center
    if not math.isfinite(radius) or radius <= 0 or not math.isfinite(sweep) or abs(sweep) < 1e-8:
        raise SampleError("REFERENCE_GEOMETRY_INVALID")
    if abs(np.linalg.norm(radial) - radius) > 1e-4 or abs(np.dot(radial, axis)) > 1e-4:
        raise SampleError("REFERENCE_GEOMETRY_INVALID", "arc start does not match radius/plane")
    steps = max(1, math.ceil(abs(sweep) * radius / spacing))
    points = []
    for i in range(steps + 1):
        angle = sweep * i / steps
        rotated = (radial * math.cos(angle) + np.cross(axis, radial) * math.sin(angle) +
                   axis * np.dot(axis, radial) * (1 - math.cos(angle)))
        points.append(center + rotated)
    return np.stack(points)


def orthogonal_frame(tangent, normal):
    t = unit(vec(tangent))
    n = vec(normal)
    u = unit(n - np.dot(n, t) * t)
    return t, u, unit(np.cross(t, u))


def angle_degrees(a, b):
    return math.degrees(math.acos(float(np.clip(np.dot(unit(a), unit(b)), -1, 1))))
