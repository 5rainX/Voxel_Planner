"""Versioned, non-destructive enrichment of fine-planner accepted records."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from path_transformer.dataset.geometry import SampleError, unit, vec
from .geometry import SPACING_VOXEL, angle_degrees, orthogonal_frame, sample_explicit_arc, sample_linear
from .validate_reference_geometry import VERSION, validate_reference


GENERATOR_VERSION = "A06.1-v1"


def _segment_points(segment, direct_segment=None):
    kind = segment.get("type")
    start, end = vec(segment.get("start")), vec(segment.get("end"))
    if direct_segment is not None:
        poses = direct_segment.get("poses", [])
        if not isinstance(poses, list) or len(poses) < 2 or not all(isinstance(pose, dict) for pose in poses):
            raise SampleError("PLANNER_EXPORT_SCHEMA_INVALID", "poses")
        points = np.asarray([pose["position"] for pose in poses], dtype=float)
        if (points.ndim != 2 or points.shape[1] != 3 or len(points) < 2 or
                not np.isfinite(points).all()):
            raise SampleError("PLANNER_EXPORT_GEOMETRY_INVALID", "dense positions")
        steps = np.linalg.norm(np.diff(points, axis=0), axis=1)
        if (np.any(steps <= 1e-9) or np.any(steps > SPACING_VOXEL + 1e-6) or
                np.linalg.norm(points[0] - start) > 1e-4 or
                np.linalg.norm(points[-1] - end) > 1e-4):
            raise SampleError("PLANNER_EXPORT_GEOMETRY_INVALID", "spacing or endpoints")
        for side, pose in (("start", poses[0]), ("end", poses[-1])):
            if angle_degrees(vec(pose["tangent"]), vec(segment[side + "_tangent"])) > 5.0:
                raise SampleError("PLANNER_EXPORT_GEOMETRY_INVALID", "endpoint tangent")
        return points, "DIRECT_PLANNER_SAMPLES", poses
    if kind in ("STRAIGHT", "TWIST"):
        chord = float(np.linalg.norm(end - start))
        length = float(segment.get("length", float("nan")))
        if chord < 1e-8 or not math.isfinite(length) or abs(chord - length) > 1e-4:
            raise SampleError("REFERENCE_GEOMETRY_UNDERDETERMINED", "straight/twist length differs from chord")
        points = sample_linear(start, end)
        direction = unit(end - start)
        if max(angle_degrees(direction, vec(segment[side + "_tangent"])) for side in ("start", "end")) > 5.0:
            raise SampleError("REFERENCE_TANGENT_INCONSISTENT", kind)
        return points, "LINEAR_STRAIGHT_INTERPOLATION" if kind == "STRAIGHT" else "TWIST_CENTERLINE_CHORD_EQUALS_LENGTH", None
    arc = segment.get("exact_arc_parameters")
    if kind in ("VERTICAL_BEND", "HORIZONTAL_BEND") and isinstance(arc, dict):
        if arc.get("source") != "DIRECT_PLANNER_ARC_PARAMETERS":
            raise SampleError("REFERENCE_GEOMETRY_UNDERDETERMINED", "untrusted arc parameter source")
        axis = unit(vec(arc["axis"]))
        sign = math.copysign(1.0, float(arc["sweep_angle_rad"]))
        if (angle_degrees(sign * np.cross(axis, start - vec(arc["center"])), vec(segment["start_tangent"])) > 5.0 or
            angle_degrees(sign * np.cross(axis, end - vec(arc["center"])), vec(segment["end_tangent"])) > 5.0):
            raise SampleError("REFERENCE_ARC_PARAMETER_MISMATCH", "arc tangents disagree with export")
        points = sample_explicit_arc(start, arc["center"], arc["axis"], float(arc["radius"]),
                                     float(arc["sweep_angle_rad"]), SPACING_VOXEL / 2)
        if np.linalg.norm(points[-1] - end) > 1e-4 or abs(abs(float(arc["radius"]) * float(arc["sweep_angle_rad"])) - float(segment["length"])) > 1e-3:
            raise SampleError("REFERENCE_ARC_PARAMETER_MISMATCH")
        return points, "EXACT_PLANNER_ARC_PARAMETERS", None
    raise SampleError("REFERENCE_GEOMETRY_UNDERDETERMINED",
                      f"{kind}: missing planner curve samples or center/radius/axis/sweep; endpoint chord prohibited")


def enrich_row(source, source_path, source_line, planner_export_root=None):
    sample_id = source["sample_id"]
    source_root = source_path.parent
    raw_relative = Path(source["artifacts"]) / f"busbar_path_candidate_{int(source['candidate_id']):03d}.txt"
    raw_path = (source_root / raw_relative).resolve()
    if not raw_path.is_relative_to(source_root) or not raw_path.is_file():
        raise SampleError("RAW_FINE_EXPORT_NOT_FOUND")
    raw_hash = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    segments = source.get("fine_segments", [])
    if not segments:
        raise SampleError("REFERENCE_PATH_EMPTY")
    direct_segments = None
    sidecar_path = None
    if planner_export_root is not None:
        export_root = Path(planner_export_root).resolve()
        sidecar_path = (export_root / raw_relative.with_suffix(".reference.json")).resolve()
        if not sidecar_path.is_relative_to(export_root):
            raise SampleError("PLANNER_EXPORT_PATH_INVALID")
        if sidecar_path.is_file():
            export = json.loads(sidecar_path.read_text(encoding="utf-8"))
            if not isinstance(export, dict):
                raise SampleError("PLANNER_EXPORT_SCHEMA_INVALID", "top-level object")
            if (export.get("export_version") != "fine_planner_reference_v1" or
                    export.get("sample_id") != sample_id or
                    export.get("map_sha256") != source["map_sha256"] or
                    export.get("raw_fine_export_sha256") != raw_hash):
                raise SampleError("PLANNER_EXPORT_SOURCE_MISMATCH")
            direct_segments = export.get("segments")
            if not isinstance(direct_segments, list) or len(direct_segments) != len(segments):
                raise SampleError("PLANNER_EXPORT_GEOMETRY_INVALID", "segment count")
            for index, (direct, segment) in enumerate(zip(direct_segments, segments)):
                if not isinstance(direct, dict):
                    raise SampleError("PLANNER_EXPORT_SCHEMA_INVALID", "segment object")
                if direct.get("segment_index") != index or direct.get("segment_type") != segment["type"]:
                    raise SampleError("PLANNER_EXPORT_SOURCE_MISMATCH", "segment identity")
    points, pose_rows, metadata = [], [], []
    max_junction_position_error = 0.0
    max_junction_tangent_discontinuity = 0.0
    for index, segment in enumerate(segments):
        dense, method, direct_poses = _segment_points(segment, direct_segments[index] if direct_segments else None)
        if points:
            gap = float(np.linalg.norm(points[-1] - dense[0]))
            max_junction_position_error = max(max_junction_position_error, gap)
            if gap > 1e-4:
                raise SampleError("REFERENCE_DISCONTINUITY", f"segment {index} gap {gap:.6g}")
            discontinuity = angle_degrees(pose_rows[-1]["tangent"], unit(dense[1] - dense[0]))
            max_junction_tangent_discontinuity = max(max_junction_tangent_discontinuity, discontinuity)
            if discontinuity > 5.0:
                raise SampleError("REFERENCE_TANGENT_DISCONTINUITY", f"segment {index}: {discontinuity:.6g} deg")
        start_index = len(points) - 1 if points else 0
        normal_start, normal_end = vec(segment["start_normal"]), vec(segment["end_normal"])
        constant_orientation = (direct_poses is None and segment["type"] == "STRAIGHT" and
                                angle_degrees(normal_start, normal_end) < 0.1)
        for j, position in enumerate(dense):
            if index > 0 and j == 0:
                # Both source segments explicitly share this junction.
                previous = pose_rows[-1]
                next_normal = direct_poses[0].get("normal") if direct_poses is not None else normal_start
                next_valid = direct_poses[0].get("orientation_valid") is True if direct_poses is not None else True
                if (previous["orientation_valid"] and
                        (not next_valid or next_normal is None or
                         angle_degrees(previous["normal"], next_normal) > 0.1)):
                    previous["normal"] = previous["binormal"] = None
                    previous["orientation_valid"] = False
                continue
            if direct_poses is not None:
                direct_pose = direct_poses[j]
                tangent = unit(vec(direct_pose["tangent"]))
                orientation_valid = direct_pose.get("orientation_valid") is True
                normal_value = direct_pose.get("normal")
                binormal_value = direct_pose.get("binormal")
            else:
                tangent = unit(dense[min(j + 1, len(dense) - 1)] - dense[max(j - 1, 0)])
                orientation_valid = constant_orientation or j in (0, len(dense) - 1)
                normal = normal_start if constant_orientation or j == 0 else normal_end
                if orientation_valid:
                    t, u, v = orthogonal_frame(tangent, normal)
                    normal_value, binormal_value = u.tolist(), v.tolist()
                else:
                    normal_value = binormal_value = None
            points.append(position)
            pose_rows.append({"position": position.tolist(), "tangent": tangent.tolist(),
                              "normal": normal_value, "binormal": binormal_value,
                              "orientation_valid": orientation_valid,
                              "segment_index": index, "segment_type": segment["type"]})
        metadata.append({"segment_index": index, "segment_type": segment["type"],
                         "start_index": start_index, "end_index": len(points) - 1,
                         "geometry_source": "DIRECT_PLANNER_SAMPLES" if direct_poses is not None else "EXPORTED_SEGMENT_SUMMARY",
                         "geometry_reconstruction_method": method,
                         "orientation_source": ("DIRECT_PLANNER_POSES" if direct_poses is not None else
                                                "EXPORTED_ENDPOINT_CONSTANT" if constant_orientation else
                                                "EXPORTED_ENDPOINTS_ONLY")})
    positions = np.stack(points)
    edge_lengths = np.linalg.norm(np.diff(positions, axis=0), axis=1)
    if not np.all(edge_lengths > 1e-9):
        raise SampleError("REFERENCE_ZERO_LENGTH_STEP")
    s = np.r_[0.0, np.cumsum(edge_lengths)]
    total = float(s[-1])
    for pose, arc in zip(pose_rows, s):
        pose["s_voxel"] = float(arc)
        pose["s_normalized"] = float(arc / total)
    record = {
        "sample_id": sample_id, "reference_geometry_version": VERSION, "status": "VALID",
        "source": {"planner": "collaborator_fine_planner",
                   "source_kind": ("LEGACY_TXT_PLANNER_SIDECAR" if direct_segments is not None
                                   else "LEGACY_TXT_SUMMARY"),
                   "source_record": str(source_path.as_posix()), "source_line": source_line,
                   "source_sample_id": sample_id, "scene_id": source["case_id"],
                   "map_id": source["map_sha256"], "raw_fine_export": str(raw_path.as_posix()),
                   "raw_fine_export_sha256": raw_hash, "generator_version": GENERATOR_VERSION,
                   "planner_reference_export": str(sidecar_path.as_posix()) if direct_segments is not None else None,
                   "planner_reference_export_sha256": (hashlib.sha256(sidecar_path.read_bytes()).hexdigest()
                                                       if direct_segments is not None else None)},
        "reference_centerline": positions.tolist(), "reference_poses": pose_rows,
        "segments": metadata,
        "validation": {"point_count": len(positions), "length_voxel": total,
                       "max_spacing_voxel": float(max(edge_lengths)),
                       "max_junction_position_error_voxel": max_junction_position_error,
                       "max_junction_tangent_discontinuity_deg": max_junction_tangent_discontinuity,
                       "orientation_valid_nodes": sum(p["orientation_valid"] for p in pose_rows),
                       "sampling_rule": "max_spacing_0.25_voxel"}}
    validate_reference(record, source)
    return record


def enrich(input_path, output, planner_export_root=None):
    input_path = Path(input_path).resolve()
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"output exists: {output}")
    output.mkdir(parents=True)
    valid = rejected = total = dense_points = orientation_valid_nodes = 0
    bend_sources = bend_recovered = 0
    reasons = {}
    with input_path.open(encoding="utf-8") as source_file, \
         (output / "accepted_enriched.jsonl").open("x", encoding="utf-8") as target_file, \
         (output / "rejected.jsonl").open("x", encoding="utf-8") as rejected_file:
        for line, raw in enumerate(source_file, 1):
            if not raw.strip():
                continue
            total += 1
            row = json.loads(raw)
            has_bend = any("BEND" in segment.get("type", "") for segment in row.get("fine_segments", []))
            bend_sources += has_bend
            try:
                reference = enrich_row(row, input_path, line, planner_export_root)
                valid += 1
                bend_recovered += has_bend
                dense_points += reference["validation"]["point_count"]
                orientation_valid_nodes += reference["validation"]["orientation_valid_nodes"]
            except (SampleError, KeyError, TypeError, ValueError) as exc:
                reason = exc.code if isinstance(exc, SampleError) else "SOURCE_SCHEMA_INVALID"
                detail = str(exc)
                reference = {"sample_id": row.get("sample_id"),
                             "reference_geometry_version": VERSION, "status": "UNRECOVERABLE",
                             "reason": reason, "details": detail,
                             "source": {"source_record": str(input_path.as_posix()), "source_line": line,
                                        "map_id": row.get("map_sha256")}}
                rejected_file.write(json.dumps(reference, allow_nan=False) + "\n")
                rejected += 1
                reasons[reason] = reasons.get(reason, 0) + 1
            row[VERSION] = reference
            target_file.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + "\n")
    stats = {"source_samples": total, "enriched_valid": valid, "unrecoverable": rejected,
             "rejection_reasons": reasons, "sampling_max_spacing_voxel": SPACING_VOXEL,
             "dense_points": dense_points, "orientation_valid_nodes": orientation_valid_nodes,
             "bend_source_samples": bend_sources, "bend_recovered": bend_recovered}
    (output / "stats.json").write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--planner-export-root", type=Path,
                        help="Optional root of planner-generated *.reference.json sidecars")
    args = parser.parse_args()
    print(json.dumps(enrich(args.input, args.output, args.planner_export_root)))


if __name__ == "__main__":
    main()
