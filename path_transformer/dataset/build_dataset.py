"""Build A06 dataset from immutable integration-test accepted.jsonl."""

import argparse
import bisect
import hashlib
import json
import math
from pathlib import Path
import statistics

import numpy as np

from .geometry import SampleError, align_paths, frames, reference_polyline, signed_theta, target, vec
from .sdf_features import features as sdf_features, field_from_occupancy, read_sdf_binary, read_text_map
from .split_dataset import split_samples
from .validate_dataset import validate_directory, validate_sample
from path_transformer.reference.validate_reference_geometry import VERSION as REFERENCE_VERSION, validate_reference


SCHEMA_VERSION = "path_transformer_dataset_v1"
GENERATOR_VERSION = "A06-v1"
SOURCE = "collaborator_fine_planner"


def percentile(values, q):
    if not values:
        return None
    return float(np.percentile(values, q))


def summary(values, *, absolute=False):
    values = [abs(float(v)) if absolute else float(v) for v in values]
    if not values:
        return {"min": None, "mean": None, "median": None, "max": None, "std": None, "p95": None}
    return {"min": min(values), "mean": statistics.fmean(values),
            "median": statistics.median(values), "max": max(values),
            "std": statistics.pstdev(values), "p95": percentile(values, 95)}


def request_number(source, name):
    args = source.get("request", {}).get("args", [])
    try:
        value = float(args[args.index(name) + 1])
    except (ValueError, IndexError, TypeError) as exc:
        raise SampleError("MISSING_REQUEST_PARAMETER", name) from exc
    if not math.isfinite(value) or value <= 0:
        raise SampleError("MISSING_REQUEST_PARAMETER", name)
    return value


def safe_source_map(root, source):
    relative = Path(source["map_file"])
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise SampleError("MAP_PATH_INVALID")
    return path


def _reference_normal_at(reference, s):
    """Return a sourced section normal only within a valid pose interval."""
    poses = reference["reference_poses"]
    values = [pose["s_normalized"] for pose in poses]
    right = bisect.bisect_left(values, s)
    if right < len(poses) and abs(values[right] - s) < 1e-8:
        return poses[right]["normal"] if poses[right]["orientation_valid"] else None
    if right == 0 or right == len(poses):
        return None
    left = right - 1
    if not (poses[left]["orientation_valid"] and poses[right]["orientation_valid"] and
            poses[left]["segment_index"] == poses[right]["segment_index"]):
        return None
    fraction = (s - values[left]) / (values[right] - values[left])
    a, b = np.asarray(poses[left]["normal"]), np.asarray(poses[right]["normal"])
    interpolated = (1 - fraction) * a + fraction * b
    return interpolated.tolist() if np.linalg.norm(interpolated) > 1e-8 else None


def _unrepresentable_native_twist_at(reference, s):
    """A bounded per-node theta cannot supervise a complete larger turn."""
    if reference is None or reference.get("source", {}).get("source_kind") != "NATIVE_FINE_PLANNER":
        return False
    poses = reference["reference_poses"]
    for segment in reference["segments"]:
        if segment["segment_type"] != "STRAIGHT_TWIST" or \
                abs(segment["twist_angle_rad"]) <= math.pi / 2 + 1e-6:
            continue
        begin = poses[segment["start_index"]]["s_normalized"]
        end = poses[segment["end_index"]]["s_normalized"]
        if begin - 1e-8 <= s <= end + 1e-8:
            return True
    return False


def make_sample(source, source_path, source_line, field, environment_path):
    coarse = np.asarray(source["coarse_centerline"], dtype=np.float64)
    enriched = source.get(REFERENCE_VERSION)
    if enriched is not None:
        if enriched.get("status") != "VALID":
            raise SampleError(enriched.get("reason", "REFERENCE_GEOMETRY_UNDERDETERMINED"),
                              enriched.get("details", "enriched reference unavailable"))
        validate_reference(enriched, source)
        reference = np.asarray(enriched["reference_centerline"], dtype=np.float64)
    else:
        reference = reference_polyline(source["fine_segments"])
    aligned, s_values, coarse_length, reference_length = align_paths(coarse, reference)
    poses = source["coarse_poses"]
    if len(poses) != len(coarse):
        raise SampleError("FRAME_COUNT_MISMATCH")
    local_frames = frames(poses)
    # Native map resolution is 1 voxel; inference fixes endpoints exactly.
    endpoint_mismatch = max(float(np.linalg.norm(coarse[0] - aligned[0])),
                            float(np.linalg.norm(coarse[-1] - aligned[-1])))
    if endpoint_mismatch > 0.5:
        raise SampleError("ENDPOINT_MISMATCH", f"{endpoint_mismatch:.6g} voxel")
    request = source.get("request", {})
    start, goal = vec(request.get("start")), vec(request.get("goal"))
    if max(np.linalg.norm(coarse[0] - start), np.linalg.norm(coarse[-1] - goal)) > 0.5:
        raise SampleError("ENDPOINT_MISMATCH", "coarse endpoints differ from request")
    width = request_number(source, "--width")
    thickness = request_number(source, "--thickness")
    if not source.get("map_sha256") or not source.get("sample_id"):
        raise SampleError("PROVENANCE_MISSING")
    all_straight = enriched is None and all(segment["type"] == "STRAIGHT" for segment in source["fine_segments"])
    normal = None
    if all_straight:
        normals = [vec(segment[side + "_normal"]) for segment in source["fine_segments"]
                   for side in ("start", "end")]
        if all(float(np.dot(normals[0], n)) / (np.linalg.norm(normals[0]) * np.linalg.norm(n)) > 1 - 1e-5 for n in normals):
            normal = normals[0]
    nodes, targets, reconstruction_errors, longitudinal, spatial_offsets = [], [], [], [], []
    for i, (point, reference_point, frame) in enumerate(zip(coarse, aligned, local_frames)):
        t, u, v = frame
        du, dv, long_residual, error = target(point, reference_point, frame)
        if i in (0, len(coarse) - 1):
            du = dv = 0.0
            error = float(np.linalg.norm(point - reference_point))
        reference_normal = (_reference_normal_at(enriched, float(s_values[i]))
                            if enriched is not None else normal)
        theta = signed_theta(frame, reference_normal) if reference_normal is not None else None
        if theta is not None and (abs(theta) > math.pi / 2 + 1e-6 or
                                  _unrepresentable_native_twist_at(enriched, float(s_values[i]))):
            theta = None  # outside the target range or ambiguous cumulative turn
        sdf = sdf_features(field, point, u, v, width, thickness)
        nodes.append({"position": point.tolist(), "tangent": t.tolist(), "normal": u.tolist(),
                      "binormal": v.tolist(), "s_normalized": float(s_values[i]), "sdf": sdf})
        targets.append({"du": du, "dv": dv, "theta": theta, "theta_valid": theta is not None,
                        "longitudinal_residual": long_residual})
        reconstruction_errors.append(error)
        longitudinal.append(abs(long_residual))
        spatial_offsets.append(float(np.linalg.norm(reference_point - point)))
    # Reject when the unmodelled along-track discrepancy exceeds one map cell.
    # A discrepancy above half a cell remains visible as a quality warning.
    if max(reconstruction_errors) > 1.0:
        raise SampleError("LONGITUDINAL_RESIDUAL_EXCEEDS_ONE_VOXEL",
                          f"max reconstruction error {max(reconstruction_errors):.6g}")
    try:
        source_file = str(source_path.relative_to(Path.cwd()).as_posix())
    except ValueError:
        source_file = str(source_path.as_posix())
    result = {
        "schema_version": SCHEMA_VERSION, "sample_id": source["sample_id"],
        "provenance": {"source": (enriched["source"]["planner"] if enriched is not None else SOURCE), "source_file": source_file,
                       "source_line": source_line, "source_sample_id": source["sample_id"],
                       "scene_id": source.get("case_id"), "map_id": source["map_sha256"],
                       "map_file": source["map_file"],
                       "generator_version": GENERATOR_VERSION, "reference_type": "supervision_candidate"},
        "environment": {"map_sha256": source["map_sha256"], "sdf_artifact": environment_path,
                        "sdf_shape_zyx": list(field.shape), "origin": [0, 0, 0],
                        "voxel_size": [1, 1, 1], "layout": "z_y_x_flat", "dtype": "float32_le",
                        "sdf_valid": True,
                        "sdf_artifact_sha256": hashlib.sha256(field.astype("<f4").tobytes()).hexdigest(),
                        "empty_map_sentinel": bool(np.all(field == 1.0e30))},
        "busbar": {"width": width, "thickness": thickness, "units": "voxel"},
        "start": {"position": start.tolist()}, "goal": {"position": goal.tolist()},
        "coarse_path": nodes,
        "reference": {"source": (enriched["source"]["planner"] if enriched is not None else SOURCE), "role": "supervision_candidate",
                      "reference_geometry_version": REFERENCE_VERSION if enriched is not None else None,
                      "independent_geometry_validation": bool(source.get("independent_geometry_validation", False)),
                      "global_optimality_proven": bool(source.get("global_optimality_proven", False)),
                      "centerline": reference.tolist(), "aligned_centerline": aligned.tolist(),
                      "source_segment_types": ([segment["segment_type"] for segment in enriched["segments"]]
                                               if enriched is not None else
                                               [segment["type"] for segment in source["fine_segments"]])},
        "targets": targets,
        "validation": {"geometry_label_valid": True, "full_feature_valid": True,
                       "endpoint_mismatch_voxel": endpoint_mismatch,
                       "coarse_length_voxel": coarse_length, "reference_length_voxel": reference_length,
                       "alignment_offset_voxel": summary(spatial_offsets),
                       "reconstruction_error_voxel": summary(reconstruction_errors),
                       "longitudinal_residual_abs_voxel": summary(longitudinal),
                       "theta_valid_count": sum(t["theta_valid"] for t in targets),
                       "sdf_invalid_query_count": 0,
                       "quality_threshold_voxel": 1.0,
                       "quality_warnings": (["LONGITUDINAL_RESIDUAL_ABOVE_HALF_VOXEL"]
                                            if max(reconstruction_errors) > 0.5 else [])}}
    validate_sample(result)
    return result


def build(input_path, output, seed=42, dry_run=False, map_root=None):
    input_path = Path(input_path).resolve()
    root = Path(map_root).resolve() if map_root is not None else input_path.parent
    output = Path(output).resolve()
    if output.exists() and not dry_run:
        raise FileExistsError(f"output already exists: {output}")
    sources, rejected = [], []
    for line, raw in enumerate(input_path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        try:
            sources.append((line, json.loads(raw)))
        except json.JSONDecodeError as exc:
            rejected.append({"sample_id": None, "source_line": line,
                             "reason": "SOURCE_JSON_INVALID", "details": str(exc)})
    seen_ids, cache, samples = set(), {}, []
    for line, source in sources:
        sample_id = source.get("sample_id")
        try:
            if sample_id in seen_ids:
                raise SampleError("DUPLICATE_SAMPLE_ID")
            seen_ids.add(sample_id)
            map_hash = source["map_sha256"]
            if map_hash not in cache:
                occupancy = read_text_map(safe_source_map(root, source), map_hash)
                cache[map_hash] = field_from_occupancy(occupancy)
            field = cache[map_hash]
            sample = make_sample(source, input_path, line, field, f"environments/{map_hash}.bin")
            samples.append(sample)
        except (SampleError, KeyError, TypeError, ValueError) as exc:
            rejected.append({"sample_id": sample_id, "source_line": line,
                             "reason": exc.code if isinstance(exc, SampleError) else "SOURCE_SCHEMA_INVALID",
                             "details": str(exc)})
    if not samples:
        raise ValueError("no valid samples produced")
    splits, split_status = split_samples(samples, seed)
    all_targets = [t for sample in samples for t in sample["targets"]]
    all_nodes = [n for sample in samples for n in sample["coarse_path"]]
    nonzero_threshold_voxel = 1.0e-4
    nonzero_theta_threshold_rad = 1.0e-4
    sample_by_id = {sample["sample_id"]: sample for sample in samples}
    split_labels = {}
    for split_name, ids in splits.items():
        rows = [sample_by_id[sample_id] for sample_id in ids]
        split_labels[split_name] = {
            "samples": len(rows),
            "nonzero_du_samples": sum(any(abs(t["du"]) > nonzero_threshold_voxel for t in row["targets"]) for row in rows),
            "nonzero_dv_samples": sum(any(abs(t["dv"]) > nonzero_threshold_voxel for t in row["targets"]) for row in rows),
            "nonzero_theta_samples": sum(any(t["theta_valid"] and abs(t["theta"]) > nonzero_theta_threshold_rad for t in row["targets"]) for row in rows)}
    valid_theta = [t["theta"] for t in all_targets if t["theta_valid"]]
    stats = {"total_source_samples": len(sources) + sum(r["reason"] == "SOURCE_JSON_INVALID" for r in rejected),
             "accepted_source_samples": len(sources),
             "dataset_samples_generated": len(samples), "dataset_samples_rejected": len(rejected),
             "geometry_aligned_samples": len(samples), "du_dv_valid_samples": len(samples),
             "theta_valid_samples": sum(all(t["theta_valid"] for t in s["targets"]) for s in samples),
             "theta_valid_nodes": sum(t["theta_valid"] for t in all_targets),
             "theta_invalid_nodes": sum(not t["theta_valid"] for t in all_targets),
             "theta_nonzero_nodes": sum(t["theta_valid"] and abs(t["theta"]) > nonzero_theta_threshold_rad for t in all_targets),
             "theta_rad": summary(valid_theta), "theta_abs_rad": summary(valid_theta, absolute=True),
             "sdf_valid_samples": len(samples), "fully_valid_training_samples": len(samples),
             "nonzero_threshold_voxel": nonzero_threshold_voxel,
             "nonzero_correction_samples": sum(any(abs(t["du"]) > nonzero_threshold_voxel or
                                                    abs(t["dv"]) > nonzero_threshold_voxel for t in row["targets"]) for row in samples),
             "nonzero_correction_nodes": sum(abs(t["du"]) > nonzero_threshold_voxel or
                                             abs(t["dv"]) > nonzero_threshold_voxel for t in all_targets),
             "split_label_coverage": split_labels,
             "unique_source_maps": len({x["map_sha256"] for _, x in sources}),
             "unique_generated_maps": len({s["provenance"]["map_id"] for s in samples}),
             "node_count": summary([len(s["coarse_path"]) for s in samples]),
             "coarse_path_length_voxel": summary([s["validation"]["coarse_length_voxel"] for s in samples]),
             "du": summary([t["du"] for t in all_targets]), "du_abs": summary([t["du"] for t in all_targets], absolute=True),
             "dv": summary([t["dv"] for t in all_targets]), "dv_abs": summary([t["dv"] for t in all_targets], absolute=True),
             "longitudinal_residual_abs_voxel": summary([t["longitudinal_residual"] for t in all_targets], absolute=True),
             "reconstruction_error_voxel": summary([t["longitudinal_residual"] for t in all_targets], absolute=True),
             "alignment_offset_voxel": summary([s["validation"]["alignment_offset_voxel"]["max"] for s in samples]),
             "quality_warning_samples": sum(bool(s["validation"]["quality_warnings"]) for s in samples),
             "sdf_center": summary([n["sdf"]["center"] for n in all_nodes]),
             "sdf_unbounded_samples": sum(s["environment"]["empty_map_sentinel"] for s in samples),
             "sdf_invalid_query_count": 0,
             "rejection_reasons": {reason: sum(r["reason"] == reason for r in rejected)
                                   for reason in sorted({r["reason"] for r in rejected})}}
    manifest = {"schema_version": SCHEMA_VERSION, "source": str(input_path.as_posix()),
                "source_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
                "split": {"strategy": "group_by_map_sha256", "seed": seed, "status": split_status,
                          "insufficient_unique_groups": len(cache) < 30},
                "samples": [{"sample_id": s["sample_id"], "map_id": s["provenance"]["map_id"],
                             "path": "samples/" + hashlib.sha256(s["sample_id"].encode()).hexdigest()[:20] + ".json"}
                            for s in samples]}
    if not dry_run:
        output.mkdir(parents=True)
        for directory in ("samples", "splits", "environments"):
            (output / directory).mkdir()
        for record, sample in zip(manifest["samples"], samples):
            (output / record["path"]).write_text(json.dumps(sample, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        for map_hash in {s["provenance"]["map_id"] for s in samples}:
            cache[map_hash].astype("<f4").tofile(output / "environments" / f"{map_hash}.bin")
        for name, ids in splits.items():
            (output / "splits" / f"{name}.txt").write_text("".join(i + "\n" for i in ids), encoding="utf-8")
        (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        (output / "stats.json").write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
        (output / "rejected.jsonl").write_text("".join(json.dumps(row) + "\n" for row in rejected), encoding="utf-8")
        report = ["# A06 Dataset V1 statistics", "",
                  f"Source rows: {stats['total_source_samples']}",
                  f"Generated: {len(samples)}; rejected: {len(rejected)}",
                  f"Theta-valid samples: {stats['theta_valid_samples']}; SDF-valid samples: {stats['sdf_valid_samples']}",
                  f"Unbounded empty-map SDF samples: {stats['sdf_unbounded_samples']}",
                  f"Quality-warning samples: {stats['quality_warning_samples']}",
                  f"Groups: {stats['unique_generated_maps']}; split status: {split_status}; seed: {seed}",
                  f"Train/val/test samples: {len(splits['train'])}/{len(splits['val'])}/{len(splits['test'])}",
                  f"Max reconstruction error (voxel): {stats['reconstruction_error_voxel']['max']:.6g}",
                  f"P95 absolute longitudinal residual (voxel): {stats['longitudinal_residual_abs_voxel']['p95']:.6g}",
                  "", "## Rejections", ""]
        report.extend(f"- {row['sample_id']}: {row['reason']} ({row['details']})" for row in rejected)
        report += ["", "The split is provisional and the data is not sufficient for a production model.", ""]
        (output / "report.md").write_text("\n".join(report), encoding="utf-8")
        validate_directory(output)
    return stats, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--map-root", type=Path, help="Directory containing source map_file artifacts")
    args = parser.parse_args()
    stats, manifest = build(args.input, args.output, args.seed, args.dry_run, args.map_root)
    print(json.dumps({"generated": stats["dataset_samples_generated"],
                      "rejected": stats["dataset_samples_rejected"],
                      "split_status": manifest["split"]["status"]}))


if __name__ == "__main__":
    main()
