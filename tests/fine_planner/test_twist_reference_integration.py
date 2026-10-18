"""Check native straight twist through Reference V1 and Dataset V1."""

import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from path_transformer.dataset.build_dataset import build
from path_transformer.dataset.geometry import SampleError
from path_transformer.dataset.validate_dataset import validate_directory
from path_transformer.reference.validate_reference_geometry import validate_file, validate_reference


def main(executable):
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        map_bytes = b"map: 16 16 16\n"
        (root / "map.txt").write_bytes(map_bytes)
        map_hash = hashlib.sha256(map_bytes).hexdigest()
        output = root / "twist_reference.json"
        subprocess.run([str(executable), str(output), map_hash, "twist"], check=True)
        reference = json.loads(output.read_text(encoding="utf-8"))
        assert validate_reference(reference)
        segment = reference["segments"][0]
        assert segment["segment_type"] == "STRAIGHT_TWIST"
        assert abs(segment["twist_angle_rad"] - math.pi / 2) < 1e-8
        assert abs(segment["twist_length_voxel"] - 4.0) < 1e-8
        assert segment["twist_profile"] == "LINEAR"
        poses = reference["reference_poses"]
        middle = poses[len(poses) // 2]
        assert abs(middle["normal"][1] + math.sqrt(0.5)) < 1e-6
        assert abs(middle["normal"][2] - math.sqrt(0.5)) < 1e-6
        assert all(pose["orientation_valid"] for pose in poses)
        damaged = json.loads(json.dumps(reference))
        damaged["segments"][0]["twist_angle_rad"] *= -1
        try:
            validate_reference(damaged)
        except SampleError:
            pass
        else:
            raise AssertionError("wrong signed twist was accepted")

        damaged = json.loads(json.dumps(reference))
        damaged["reference_poses"][1]["segment_type"] = "STRAIGHT"
        try:
            validate_reference(damaged)
        except SampleError:
            pass
        else:
            raise AssertionError("wrong twist pose segment type was accepted")

        row = {
            "sample_id": reference["sample_id"], "case_id": "mvp",
            "map_sha256": map_hash, "map_file": "map.txt",
            "coarse_centerline": [[2, 2, 2], [4, 2, 2], [6, 2, 2]],
            "coarse_poses": [[0, 0, 1, 1, 0, 0]] * 3,
            "request": {"start": [2, 2, 2], "goal": [6, 2, 2],
                        "args": ["--width", "0.5", "--thickness", "0.5"]},
            "independent_geometry_validation": False,
            "global_optimality_proven": False,
            "reference_geometry_v1": reference,
        }
        wrong_dimensions = json.loads(json.dumps(row))
        wrong_dimensions["request"]["args"][1] = "0.7"
        try:
            validate_reference(reference, wrong_dimensions)
        except SampleError:
            pass
        else:
            raise AssertionError("native section width mismatch was accepted")
        rows = [row]
        for case, expected_angle in (("negative90", -math.pi/2),
                                     ("positive360", 2*math.pi),
                                     ("smoothstep", 1.1)):
            case_output = root / f"{case}.json"
            subprocess.run([str(executable), str(case_output), map_hash,
                            "twist-case", case, "0.25"], check=True)
            case_reference = json.loads(case_output.read_text(encoding="utf-8"))
            assert validate_reference(case_reference)
            assert abs(case_reference["segments"][0]["twist_angle_rad"] - expected_angle) < 1e-12
            case_row = json.loads(json.dumps(row))
            case_row["sample_id"] = case_reference["sample_id"]
            case_row["reference_geometry_v1"] = case_reference
            if case == "smoothstep":
                case_row["request"]["args"] = ["--width", "0.6", "--thickness", "0.2"]
            rows.append(case_row)
        source = root / "native_twist.jsonl"
        source.write_text("".join(json.dumps(item, allow_nan=False) + "\n" for item in rows),
                          encoding="utf-8")
        assert validate_file(source) == {"total": 4, "valid": 4, "unrecoverable": 0}
        invalid_row = json.loads(json.dumps(row))
        invalid_row["sample_id"] = "mvp:twist-invalid"
        invalid_row["reference_geometry_v1"]["sample_id"] = invalid_row["sample_id"]
        invalid_row["reference_geometry_v1"]["reference_poses"][1]["segment_type"] = "STRAIGHT"
        with source.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(invalid_row, allow_nan=False) + "\n")
        stats, _ = build(source, root / "dataset")
        assert stats["dataset_samples_generated"] == 4
        assert stats["dataset_samples_rejected"] == 1
        samples = {sample["sample_id"]: sample for path in
                   (root / "dataset" / "samples").glob("*.json")
                   for sample in [json.loads(path.read_text())]}
        assert all(sample["reference"]["source_segment_types"] == ["STRAIGHT_TWIST"]
                   for sample in samples.values())
        assert any(target["theta_valid"] and target["theta"] > 0.1
                   for target in samples["mvp:twist"]["targets"])
        assert any(target["theta_valid"] and target["theta"] < -0.1
                   for target in samples["cross-validation:negative90"]["targets"])
        assert any(target["theta_valid"] and target["theta"] > 0.1
                   for target in samples["cross-validation:smoothstep"]["targets"])
        assert all(not target["theta_valid"] and target["theta"] is None
                   for target in samples["cross-validation:positive360"]["targets"])
        assert validate_directory(root / "dataset")["valid"] == 4


if __name__ == "__main__":
    main(Path(sys.argv[1]))
