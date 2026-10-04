"""Read C++ SDF binaries or derive an equivalent field from a hashed text map."""

import hashlib
import math
from pathlib import Path

import numpy as np

from .geometry import SampleError


INFINITY = 1.0e30  # GlobalSdf's finite empty/full-map sentinel


def read_text_map(path, expected_sha256):
    path = Path(path)
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise SampleError("MAP_HASH_MISMATCH", str(path))
    lines = raw.decode("utf-8-sig").splitlines()
    header = lines[0].split()
    if len(header) < 4 or header[0] not in ("map:", "header"):
        raise SampleError("MAP_FORMAT_INVALID")
    try:
        width, height, depth = map(int, header[1:4])
        if min(width, height, depth) <= 0:
            raise ValueError
        occupied = np.zeros((depth, height, width), dtype=bool)
        for line in lines[1:]:
            cells = line.split()
            if len(cells) < 3:
                continue
            x, y, z = map(int, cells[:3])
            if 0 <= x < width and 0 <= y < height and 0 <= z < depth:
                occupied[z, y, x] = True
    except ValueError as exc:
        raise SampleError("MAP_FORMAT_INVALID") from exc
    return occupied


def _edt_line(values):
    """1D lower envelope of squared distances, as in SdfVolume.cpp."""
    sites = np.flatnonzero(values < INFINITY * 0.5)
    if not len(sites):
        return np.full(len(values), INFINITY, dtype=np.float64)
    envelope = []
    starts = []
    for site in sites:
        site = int(site)
        while envelope:
            previous = envelope[-1]
            crossing = (values[site] + site * site - values[previous] - previous * previous) / (2 * (site - previous))
            if crossing > starts[-1]:
                break
            envelope.pop()
            starts.pop()
        envelope.append(site)
        starts.append(-math.inf if len(envelope) == 1 else crossing)
    out = np.empty(len(values), dtype=np.float64)
    active = 0
    for i in range(len(values)):
        while active + 1 < len(envelope) and starts[active + 1] < i:
            active += 1
        out[i] = (i - envelope[active]) ** 2 + values[envelope[active]]
    return out


def _squared_distances(features):
    field = np.where(features, 0.0, INFINITY)
    # X, Y, Z; width/height/depth are voxel counts and spacing is 1.
    for axis in (2, 1, 0):
        field = np.apply_along_axis(_edt_line, axis, field)
    return field


def field_from_occupancy(occupied):
    occupied = np.asarray(occupied, dtype=bool)
    if occupied.ndim != 3 or not all(occupied.shape):
        raise SampleError("MAP_FORMAT_INVALID")
    if not occupied.any():
        return np.full(occupied.shape, INFINITY, dtype="<f4")
    if occupied.all():
        return np.full(occupied.shape, -INFINITY, dtype="<f4")
    to_obstacles = _squared_distances(occupied)
    to_free = _squared_distances(~occupied)
    magnitude = np.sqrt(np.where(occupied, to_free, to_obstacles))
    return np.where(occupied, -1.0, 1.0) * np.maximum(0.0, magnitude - 0.5).astype("<f4")


def read_sdf_binary(path, shape):
    data = np.fromfile(path, dtype="<f4")
    if data.size != math.prod(shape):
        raise SampleError("SDF_BINARY_SIZE_MISMATCH")
    return data.reshape(shape)  # z, y, x, like DatasetWriter.h


def query(field, point):
    """C++ trilinear interior query; reject outside rather than use soft wall."""
    x, y, z = map(float, point)
    depth, height, width = field.shape
    if not all(math.isfinite(a) for a in (x, y, z)) or not (0 <= x <= width - 1 and 0 <= y <= height - 1 and 0 <= z <= depth - 1):
        raise SampleError("SDF_OUT_OF_BOUNDS")
    x0, y0, z0 = math.floor(x), math.floor(y), math.floor(z)
    x1, y1, z1 = min(x0 + 1, width - 1), min(y0 + 1, height - 1), min(z0 + 1, depth - 1)
    fx, fy, fz = x - x0, y - y0, z - z0
    c00 = float(field[z0, y0, x0]) * (1 - fx) + float(field[z0, y0, x1]) * fx
    c10 = float(field[z0, y1, x0]) * (1 - fx) + float(field[z0, y1, x1]) * fx
    c01 = float(field[z1, y0, x0]) * (1 - fx) + float(field[z1, y0, x1]) * fx
    c11 = float(field[z1, y1, x0]) * (1 - fx) + float(field[z1, y1, x1]) * fx
    return ((c00 * (1 - fy) + c10 * fy) * (1 - fz) +
            (c01 * (1 - fy) + c11 * fy) * fz)


def features(field, point, normal, binormal, width, thickness):
    point, normal, binormal = map(lambda a: np.asarray(a, dtype=float), (point, normal, binormal))
    positions = (point, point + width / 2 * normal, point - width / 2 * normal,
                 point + thickness / 2 * binormal, point - thickness / 2 * binormal)
    values = [query(field, p) for p in positions]
    return dict(zip(("center", "u_pos", "u_neg", "v_pos", "v_neg"), values))
