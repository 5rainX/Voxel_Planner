"""Arc-length alignment and local frames compatible with infer.py."""

import math

import numpy as np


class SampleError(ValueError):
    def __init__(self, code, detail=""):
        super().__init__(detail or code)
        self.code = code


def vec(value):
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (3,) or not np.isfinite(result).all():
        raise SampleError("NONFINITE_VALUE", "expected finite xyz vector")
    return result


def unit(value, code="FRAME_DEGENERATE"):
    value = np.asarray(value, dtype=np.float64)
    length = float(np.linalg.norm(value))
    if not math.isfinite(length) or length < 1e-8:
        raise SampleError(code)
    return value / length


def arclength(points, *, allow_duplicates=False):
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 2 or not np.isfinite(points).all():
        raise SampleError("NONFINITE_VALUE", "path needs at least two finite xyz points")
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    if not allow_duplicates and np.any(lengths < 1e-8):
        raise SampleError("ZERO_LENGTH_COARSE_PATH", "repeated coarse nodes")
    total = float(lengths.sum())
    if total < 1e-8:
        raise SampleError("ZERO_LENGTH_COARSE_PATH")
    return np.r_[0.0, np.cumsum(lengths)] / total, total


def align_paths(coarse, reference):
    """Interpolate reference at each coarse normalized cumulative arc length."""
    coarse_s, coarse_length = arclength(coarse)
    reference = np.asarray(reference, dtype=np.float64)
    if reference.ndim != 2 or reference.shape[1] != 3 or not np.isfinite(reference).all():
        raise SampleError("REFERENCE_ALIGNMENT_FAILED")
    # Adjacent segment endpoints may repeat; remove only those duplicates.
    keep = [0] + [i for i in range(1, len(reference))
                  if np.linalg.norm(reference[i] - reference[i - 1]) >= 1e-8]
    reference = reference[keep]
    reference_s, reference_length = arclength(reference)
    aligned = np.column_stack([
        np.interp(coarse_s, reference_s, reference[:, axis]) for axis in range(3)
    ])
    return aligned, coarse_s, coarse_length, reference_length


def reference_polyline(segments):
    """Accept only segments whose centerline is specified by their endpoints.

    A bend's end points, length and tangent do not uniquely describe its arc.
    A twist is accepted only if its documented centerline length equals its chord.
    """
    if not isinstance(segments, list) or not segments:
        raise SampleError("REFERENCE_PATH_EMPTY")
    points = []
    for segment in segments:
        kind = segment.get("type")
        if kind not in ("STRAIGHT", "TWIST"):
            raise SampleError("REFERENCE_GEOMETRY_UNDERDETERMINED", str(kind))
        start, end = vec(segment.get("start")), vec(segment.get("end"))
        chord = float(np.linalg.norm(end - start))
        length = float(segment.get("length", float("nan")))
        if chord < 1e-8 or not math.isfinite(length) or abs(chord - length) > 1e-4:
            raise SampleError("REFERENCE_GEOMETRY_UNDERDETERMINED", "segment length differs from chord")
        if points and np.linalg.norm(points[-1] - start) > 1e-4:
            raise SampleError("REFERENCE_DISCONTINUITY")
        if not points:
            points.append(start)
        points.append(end)
    return np.stack(points)


def frames(poses):
    """Same t, projected-normal, t-cross-u convention as infer.py."""
    result = []
    for pose in poses:
        pose = np.asarray(pose, dtype=np.float64)
        if pose.shape != (6,) or not np.isfinite(pose).all():
            raise SampleError("FRAME_DEGENERATE")
        t = unit(pose[3:])
        u = unit(pose[:3] - np.dot(pose[:3], t) * t)
        v = unit(np.cross(t, u))
        if result and np.dot(u, result[-1][1]) < 0:
            u, v = -u, -v
        result.append((t, u, v))
    return result


def target(coarse, aligned, frame):
    delta = aligned - coarse
    t, u, v = frame
    du, dv, longitudinal = (float(np.dot(delta, axis)) for axis in (u, v, t))
    rebuilt = coarse + du * u + dv * v
    return du, dv, longitudinal, float(np.linalg.norm(rebuilt - aligned))


def signed_theta(coarse_frame, reference_normal):
    t, u, v = coarse_frame
    normal = vec(reference_normal)
    normal = unit(normal - np.dot(normal, t) * t)
    return math.atan2(float(np.dot(normal, v)), float(np.dot(normal, u)))
