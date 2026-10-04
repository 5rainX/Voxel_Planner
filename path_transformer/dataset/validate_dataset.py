"""Validate A06 schema, geometry invariants, SDF artifacts and split leakage."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from .geometry import SampleError
from .sdf_features import query, read_sdf_binary


VERSION = "path_transformer_dataset_v1"


def require(condition, message):
    if not condition:
        raise SampleError("SCHEMA_VALIDATION_FAILED", message)


def finite_vector(value):
    a = np.asarray(value, dtype=float)
    require(a.shape == (3,) and np.isfinite(a).all(), "expected finite xyz vector")
    return a


def validate_sample(sample):
    require(sample.get("schema_version") == VERSION, "wrong schema version")
    require(isinstance(sample.get("sample_id"), str) and sample["sample_id"], "missing sample_id")
    provenance = sample["provenance"]
    for key in ("source", "source_file", "source_line", "source_sample_id", "scene_id", "map_id", "map_file", "generator_version", "reference_type"):
        require(provenance.get(key) is not None and provenance[key] != "", f"missing provenance.{key}")
    environment = sample["environment"]
    require(environment.get("sdf_valid") is True, "SDF must be valid in full sample")
    require(environment.get("map_sha256") == provenance["map_id"], "map identity mismatch")
    require(isinstance(environment.get("sdf_artifact_sha256"), str) and
            len(environment["sdf_artifact_sha256"]) == 64, "missing SDF artifact hash")
    require(len(environment["sdf_shape_zyx"]) == 3 and min(environment["sdf_shape_zyx"]) > 0, "bad SDF shape")
    require(environment.get("layout") == "z_y_x_flat" and environment.get("dtype") == "float32_le", "bad SDF metadata")
    require(sample["busbar"]["width"] > 0 and sample["busbar"]["thickness"] > 0, "bad busbar dimensions")
    coarse = sample["coarse_path"]
    targets = sample["targets"]
    aligned = sample["reference"]["aligned_centerline"]
    require(len(coarse) >= 2 and len(coarse) == len(targets) == len(aligned), "node/target/reference count mismatch")
    require(sample["reference"]["role"] == "supervision_candidate", "wrong reference role")
    finite_vector(sample["start"]["position"])
    finite_vector(sample["goal"]["position"])
    s_values = []
    previous_normal = None
    for i, (node, label, ref) in enumerate(zip(coarse, targets, aligned)):
        point = finite_vector(node["position"])
        t, u, v = (finite_vector(node[name]) for name in ("tangent", "normal", "binormal"))
        ref = finite_vector(ref)
        for axis in (t, u, v):
            require(abs(np.linalg.norm(axis) - 1) < 1e-5, "non-unit frame")
        require(max(abs(float(np.dot(t, u))), abs(float(np.dot(t, v))), abs(float(np.dot(u, v)))) < 1e-5,
                "non-orthogonal frame")
        require(np.linalg.norm(np.cross(t, u) - v) < 1e-5, "handedness mismatch")
        if previous_normal is not None:
            require(float(np.dot(previous_normal, u)) >= -1e-6, "frame normal flips")
        previous_normal = u
        s = float(node["s_normalized"])
        require(math.isfinite(s), "nonfinite s")
        s_values.append(s)
        for key in ("center", "u_pos", "u_neg", "v_pos", "v_neg"):
            require(math.isfinite(float(node["sdf"][key])), f"invalid sdf.{key}")
        for key in ("du", "dv", "longitudinal_residual"):
            require(math.isfinite(float(label[key])), f"invalid target.{key}")
        require((label["theta"] is None) == (label["theta_valid"] is False), "theta mask mismatch")
        if label["theta_valid"]:
            require(math.isfinite(float(label["theta"])) and abs(label["theta"]) <= math.pi,
                    "invalid theta")
        reconstructed = point + label["du"] * u + label["dv"] * v
        require(np.linalg.norm(reconstructed - ref) <= sample["validation"]["quality_threshold_voxel"] + 1e-6,
                "reconstruction exceeds threshold")
        if i in (0, len(coarse) - 1):
            require(label["du"] == label["dv"] == 0, "endpoint offset must be zero")
    require(abs(s_values[0]) < 1e-8 and abs(s_values[-1] - 1) < 1e-8, "s endpoints invalid")
    require(all(b > a for a, b in zip(s_values, s_values[1:])), "s not strictly increasing")
    require(sample["validation"]["geometry_label_valid"] is True and
            sample["validation"]["full_feature_valid"] is True, "validity flags invalid")
    return True


def validate_directory(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    require(manifest["schema_version"] == VERSION, "manifest version")
    ids, groups, assignments = set(), {}, {}
    for name in ("train", "val", "test"):
        file = directory / "splits" / f"{name}.txt"
        require(file.exists(), f"missing {name} split")
        for sample_id in file.read_text(encoding="utf-8").splitlines():
            require(sample_id not in assignments, "sample crosses splits")
            assignments[sample_id] = name
    for row in manifest["samples"]:
        require(row["sample_id"] not in ids, "duplicate sample_id")
        ids.add(row["sample_id"])
        sample = json.loads((directory / row["path"]).read_text(encoding="utf-8"))
        validate_sample(sample)
        require(sample["sample_id"] == row["sample_id"], "manifest sample mismatch")
        key = sample["provenance"]["map_id"]
        require(key == row["map_id"], "manifest map mismatch")
        group_split = groups.setdefault(key, assignments.get(sample["sample_id"]))
        require(group_split == assignments.get(sample["sample_id"]) and group_split is not None,
                "map leakage or unsplit sample")
        artifact = directory / sample["environment"]["sdf_artifact"]
        require(artifact.is_file(), "missing SDF artifact")
        require(hashlib.sha256(artifact.read_bytes()).hexdigest() == sample["environment"]["sdf_artifact_sha256"],
                "SDF artifact hash mismatch")
        field = read_sdf_binary(artifact, sample["environment"]["sdf_shape_zyx"])
        for node in sample["coarse_path"]:
            measured = query(field, node["position"])
            require(np.isclose(measured, node["sdf"]["center"], rtol=1e-6, atol=1e-5),
                    "SDF center does not match artifact")
    require(ids == set(assignments), "split IDs differ from manifest")
    stats = json.loads((directory / "stats.json").read_text(encoding="utf-8"))
    require(stats["dataset_samples_generated"] == len(ids), "stats sample count mismatch")
    return {"valid": len(ids), "invalid": 0,
            "warnings": [manifest["split"]["status"]] if manifest["split"]["status"] != "GROUPED" else []}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(validate_directory(args.directory)))


if __name__ == "__main__":
    main()
