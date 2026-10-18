"""Native C++ Reference Geometry V1 through the existing Python dataset flow."""

import hashlib
import json
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
        export = root / "native_reference.json"
        subprocess.run([str(executable), str(export), map_hash], check=True)
        reference = json.loads(export.read_text(encoding="utf-8"))
        assert validate_reference(reference)
        assert reference["source"]["source_kind"] == "NATIVE_FINE_PLANNER"
        assert "raw_fine_export" not in reference["source"]
        assert reference["segments"][1]["segment_type"] == "HORIZONTAL_BEND"
        assert reference["segments"][1]["radius_voxel"] == 2
        damaged = json.loads(json.dumps(reference))
        damaged["segments"][1]["radius_voxel"] = 3
        try:
            validate_reference(damaged)
        except SampleError:
            pass
        else:
            raise AssertionError("damaged native arc was accepted")
        damaged = json.loads(json.dumps(reference))
        damaged["segments"][1]["signed_sweep_rad"] *= -1
        try:
            validate_reference(damaged)
        except SampleError:
            pass
        else:
            raise AssertionError("wrong native sweep direction was accepted")

        row = {
            "sample_id": reference["sample_id"], "case_id": "mvp",
            "map_sha256": map_hash, "map_file": "map.txt",
            "coarse_centerline": [[2, 4, 4], [4, 4, 4], [6, 6, 4]],
            "coarse_poses": [[0, 0, 1, 1, 0, 0], [0, 0, 1, 1, 0, 0],
                             [0, 0, 1, 0, 1, 0]],
            "request": {"start": [2, 4, 4], "goal": [6, 6, 4],
                        "args": ["--width", "0.5", "--thickness", "0.5"]},
            "independent_geometry_validation": False,
            "global_optimality_proven": False,
            "reference_geometry_v1": reference,
        }
        source = root / "native_accepted.jsonl"
        source.write_text(json.dumps(row, allow_nan=False) + "\n", encoding="utf-8")
        assert validate_file(source) == {"total": 1, "valid": 1, "unrecoverable": 0}
        stats, manifest = build(source, root / "dataset")
        assert stats["dataset_samples_generated"] == 1
        assert stats["du_dv_valid_samples"] == 1
        assert stats["sdf_valid_samples"] == 1
        assert manifest["split"]["strategy"] == "group_by_map_sha256"
        assert validate_directory(root / "dataset")["valid"] == 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_reference_integration.py <test_fine_planner executable>")
    main(Path(sys.argv[1]))
