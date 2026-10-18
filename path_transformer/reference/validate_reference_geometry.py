"""Validate versioned dense reference geometry and source correspondence."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from path_transformer.dataset.geometry import SampleError
from .geometry import angle_degrees, orthogonal_frame


VERSION = "reference_geometry_v1"


def require(condition, message):
    if not condition:
        raise SampleError("REFERENCE_SCHEMA_INVALID", message)


def validate_reference(record, source=None):
    require(record.get("reference_geometry_version") == VERSION, "version")
    require(record.get("status") == "VALID", "status")
    require(bool(record.get("sample_id")), "sample ID")
    provenance = record.get("source", {})
    for key in ("planner", "source_record", "source_line", "scene_id", "map_id"):
        require(provenance.get(key) is not None and provenance[key] != "", f"source.{key}")
    native = provenance.get("source_kind") == "NATIVE_FINE_PLANNER"
    if native:
        require(provenance.get("planner") == "native_fine_planner" and
                provenance.get("raw_fine_export") is None, "native source identity")
        require(record.get("validation", {}).get("validation_scope") ==
                "SAMPLED_SDF_CONSERVATIVE_SECTION" and
                record["validation"].get("independent_full_size_collision_proven") is False,
                "native validation scope")
        width = record.get("section_width_voxel")
        thickness = record.get("section_thickness_voxel")
        require((width is None) == (thickness is None), "native section dimensions")
        if width is not None:
            require(all(isinstance(value, (int, float)) and math.isfinite(value) and value > 0
                        for value in (width, thickness)), "native section dimensions")
            if source is not None:
                args = source.get("request", {}).get("args", [])
                try:
                    source_width = float(args[args.index("--width") + 1])
                    source_thickness = float(args[args.index("--thickness") + 1])
                except (AttributeError, ValueError, IndexError, TypeError):
                    require(False, "native source section dimensions")
                require(math.isfinite(source_width) and math.isfinite(source_thickness) and
                        abs(source_width - width) < 1e-8 and
                        abs(source_thickness - thickness) < 1e-8,
                        "native source section dimensions")
    else:
        require(bool(provenance.get("raw_fine_export")), "source.raw_fine_export")
    points = np.asarray(record.get("reference_centerline"), dtype=float)
    poses = record.get("reference_poses", [])
    require(points.ndim == 2 and points.shape[1] == 3 and len(points) >= 2, "centerline shape")
    require(np.isfinite(points).all(), "nonfinite centerline")
    require(len(poses) == len(points), "pose count")
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    require(np.all(lengths > 1e-9), "repeated points")
    cumulative = np.r_[0.0, np.cumsum(lengths)]
    total = float(cumulative[-1])
    require(total > 0, "zero length")
    segments = record.get("segments", [])
    require(bool(segments), "segment metadata")
    owner = {i for segment in segments for i in range(segment["start_index"], segment["end_index"] + 1)}
    require(owner == set(range(len(points))), "segment coverage")
    if native:
        expected_begin = 0
        for segment_index, segment in enumerate(segments):
            begin, end = segment["start_index"], segment["end_index"]
            require(segment["segment_index"] == segment_index and
                    begin == expected_begin and begin < end < len(points),
                    "native segment indices")
            expected_begin = end
            # The shared junction pose belongs to the preceding segment.
            for index in range(begin if begin == 0 else begin + 1, end + 1):
                require(poses[index].get("segment_index") == segment["segment_index"] and
                        poses[index].get("segment_type") == segment["segment_type"],
                        "native pose segment identity")
            require(np.linalg.norm(points[begin] - segment["start"]) < 1e-4 and
                    np.linalg.norm(points[end] - segment["end"]) < 1e-4,
                    "native segment endpoints")
            fractions = np.linspace(0.0, 1.0, end - begin + 1)
            kind = segment["segment_type"]
            if kind in ("STRAIGHT", "STRAIGHT_TWIST"):
                expected = (np.asarray(segment["start"])[None, :] * (1 - fractions[:, None]) +
                            np.asarray(segment["end"])[None, :] * fractions[:, None])
                if kind == "STRAIGHT_TWIST":
                    length = float(segment["twist_length_voxel"])
                    angle = float(segment["twist_angle_rad"])
                    profile = segment["twist_profile"]
                    start_normal = np.asarray(segment["twist_start_normal"], dtype=float)
                    end_normal = np.asarray(segment["twist_end_normal"], dtype=float)
                    tangent = (points[end] - points[begin]) / np.linalg.norm(points[end] - points[begin])
                    require(math.isfinite(length) and length > 0 and
                            abs(length - np.linalg.norm(points[end] - points[begin])) < 1e-4 and
                            math.isfinite(angle) and profile in ("LINEAR", "SMOOTHSTEP") and
                            start_normal.shape == (3,) and end_normal.shape == (3,) and
                            np.isfinite(start_normal).all() and np.isfinite(end_normal).all() and
                            abs(np.linalg.norm(start_normal) - 1) < 1e-5 and
                            abs(np.linalg.norm(end_normal) - 1) < 1e-5 and
                            abs(float(start_normal @ tangent)) < 1e-5 and
                            abs(float(end_normal @ tangent)) < 1e-5,
                            "native twist parameters")
                    shaped = fractions if profile == "LINEAR" else fractions**2 * (3 - 2 * fractions)
                    angles = angle * shaped
                    expected_normals = (np.cos(angles)[:, None] * start_normal +
                                        np.sin(angles)[:, None] * np.cross(tangent, start_normal))
                    require(np.linalg.norm(expected_normals[-1] - end_normal) < 1e-4,
                            "native twist end frame")
                    for index, expected_normal in zip(range(begin, end + 1), expected_normals):
                        pose = poses[index]
                        require(pose.get("orientation_valid") is True and
                                np.linalg.norm(np.asarray(pose["normal"], dtype=float) - expected_normal) < 1e-4 and
                                np.linalg.norm(np.asarray(pose["tangent"], dtype=float) - tangent) < 1e-4,
                                "native twist frame evolution")
                    peak_rate = (1.5 if profile == "SMOOTHSTEP" else 1.0) * abs(angle) / length
                    require(abs(float(segment["peak_twist_rate_rad_per_voxel"]) - peak_rate) < 1e-6,
                            "native twist rate metadata")
            elif kind in ("HORIZONTAL_BEND", "VERTICAL_BEND"):
                center = np.asarray(segment["center"], dtype=float)
                axis = np.asarray(segment["axis"], dtype=float)
                radius = float(segment["radius_voxel"])
                sweep = float(segment["signed_sweep_rad"])
                require(center.shape == (3,) and np.isfinite(center).all() and
                        axis.shape == (3,) and np.isfinite(axis).all() and
                        abs(np.linalg.norm(axis) - 1) < 1e-5 and
                        math.isfinite(radius) and radius > 0 and math.isfinite(sweep) and
                        abs(sweep) > 1e-8, "native arc parameters")
                radial = points[begin:end + 1] - center
                require(np.max(np.abs(np.linalg.norm(radial, axis=1) - radius)) < 1e-4 and
                        np.max(np.abs(radial @ axis)) < 1e-4,
                        "native arc radius/plane")
                radial_start = radial[0]
                angles = sweep * fractions
                expected = center + (np.cos(angles)[:, None] * radial_start +
                                     np.sin(angles)[:, None] * np.cross(axis, radial_start))
            else:
                require(False, "native segment type")
            require(np.max(np.linalg.norm(points[begin:end + 1] - expected, axis=1)) < 1e-4,
                    "native segment ordered geometry")
        require(expected_begin == len(points) - 1, "native segment final coverage")
    for i, pose in enumerate(poses):
        position = np.asarray(pose.get("position"), dtype=float)
        tangent = np.asarray(pose.get("tangent"), dtype=float)
        require(position.shape == (3,) and np.isfinite(position).all() and
                np.linalg.norm(position - points[i]) < 1e-6, "pose position")
        require(tangent.shape == (3,) and np.isfinite(tangent).all() and
                abs(np.linalg.norm(tangent) - 1) < 1e-5, "pose tangent")
        numerical = points[1] - points[0] if i == 0 else (points[-1] - points[-2] if i == len(points)-1 else points[i+1] - points[i-1])
        require(angle_degrees(tangent, numerical) < 5.0, "tangent/centerline mismatch")
        require(abs(pose["s_normalized"] - cumulative[i] / total) < 1e-6, "arc length mismatch")
        if pose.get("orientation_valid"):
            t, u, v = orthogonal_frame(tangent, pose.get("normal"))
            require(np.linalg.norm(u - pose["normal"]) < 1e-5 and
                    np.linalg.norm(v - pose["binormal"]) < 1e-5, "invalid orientation frame")
        else:
            require(pose.get("normal") is None and pose.get("binormal") is None,
                    "invalid orientation must be null")
    require(abs(poses[0]["s_normalized"]) < 1e-8 and
            abs(poses[-1]["s_normalized"] - 1) < 1e-8, "arc endpoints")
    if source is not None:
        require(record["sample_id"] == source["sample_id"], "source sample ID")
        require(provenance["map_id"] == source["map_sha256"], "source map ID")
        if source.get("fine_segments"):
            start = source["fine_segments"][0]["start"]
            goal = source["fine_segments"][-1]["end"]
        else:
            require(native, "missing fine segments for legacy source")
            start = source["request"]["start"]
            goal = source["request"]["goal"]
        require(np.linalg.norm(points[0] - start) < 1e-4 and
                np.linalg.norm(points[-1] - goal) < 1e-4, "reference endpoints")
    return True


def validate_file(path):
    total = valid = rejected = 0
    ids = set()
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        total += 1
        row = json.loads(line)
        reference = row.get(VERSION)
        require(reference is not None and row.get("sample_id") not in ids, "missing/duplicate reference")
        ids.add(row["sample_id"])
        if reference.get("status") == "VALID":
            validate_reference(reference, row)
            valid += 1
        else:
            require(reference.get("status") == "UNRECOVERABLE" and bool(reference.get("reason")),
                    "unrecoverable record must have reason")
            rejected += 1
    return {"total": total, "valid": valid, "unrecoverable": rejected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    print(json.dumps(validate_file(args.path)))


if __name__ == "__main__":
    main()
