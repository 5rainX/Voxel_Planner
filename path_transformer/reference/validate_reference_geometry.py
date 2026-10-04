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
    for key in ("planner", "source_record", "source_line", "scene_id", "map_id", "raw_fine_export"):
        require(provenance.get(key) is not None and provenance[key] != "", f"source.{key}")
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
        require(np.linalg.norm(points[0] - source["fine_segments"][0]["start"]) < 1e-4 and
                np.linalg.norm(points[-1] - source["fine_segments"][-1]["end"]) < 1e-4,
                "reference endpoints")
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
