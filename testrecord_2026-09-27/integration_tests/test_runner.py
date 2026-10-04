import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from make_cases import create_cases
from run_tests import collect_result, make_command, read_segments, run_case, validate_map


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / "stdout.log").write_text("", encoding="utf-8")

    def tearDown(self):
        self.temp.cleanup()

    def test_timeout_cannot_become_no_path_or_label(self):
        (self.root / "voxel_planner_release_status.txt").write_text("status=NO_PATH\n")
        result, pairs = collect_result(self.root, None, timed_out=True)
        self.assertEqual(result["status"], "PROCESS_TIMEOUT")
        self.assertEqual(pairs, [])

    def test_overwritten_coarse_error_is_preserved(self):
        (self.root / "voxel_planner_release_status.txt").write_text("status=NO_PATH\n")
        (self.root / "voxel_planner_release_coarse_group_request_001_status.txt").write_text("status=ERROR\nmessage=search limit\n")
        result, pairs = collect_result(self.root, 1)
        self.assertEqual(result["status"], "COARSE_ERROR")
        self.assertFalse(pairs)

    def test_preprocess_error_overwritten_by_no_path(self):
        (self.root / "voxel_planner_release_status.txt").write_text("status=NO_PATH\n")
        (self.root / "stdout.log").write_text("Voxel_Planner_Release error: bad map\n")
        self.assertEqual(collect_result(self.root, 1)[0]["status"], "COARSE_ERROR")

    def test_no_outputs_never_counts_as_success(self):
        self.assertEqual(collect_result(self.root, 0)[0]["status"], "OUTPUT_ERROR")

    def test_process_crash_never_counts_as_no_path(self):
        self.assertEqual(collect_result(self.root, 3221225781)[0]["status"], "PROCESS_ERROR")

    def test_fine_timeout_is_not_no_path(self):
        (self.root / "planning_candidates.csv").write_text("fine_success,message\nfalse,PLANNER_TIMEOUT\n")
        self.assertEqual(collect_result(self.root, 1)[0]["status"], "FINE_TIMEOUT")

    def test_declared_success_without_exports_is_rejected(self):
        (self.root / "planning_candidates.csv").write_text(
            "candidate_id,coarse_tangent_group,route_index,fine_success\n1,1,1,true\n")
        with self.assertRaises(ValueError):
            collect_result(self.root, 0)

    def test_map_checks_each_coordinate_not_just_flat_index(self):
        path = self.root / "map.txt"
        path.write_text("map: 4 4 4\n5 0 0 0\n")
        with self.assertRaises(ValueError):
            validate_map(path)

    def test_truncated_fine_geometry_is_rejected(self):
        path = self.root / "path.txt"
        path.write_text("Segment 0:\n  Type: STRAIGHT\n  Start: (0, 0, 0)\n")
        with self.assertRaises(ValueError):
            read_segments(path)

    def test_command_calls_real_coarse_stage_without_silent_fallback(self):
        manifest = create_cases(self.root / "fixtures")
        case = json.loads(manifest.read_text())["cases"][0]
        command = make_command(Path("planner.exe"), case, 3, 1000, False)
        self.assertIn("--voxel-planner-release", command)
        self.assertIn("--endpoint-tangent-candidate", command)
        case["args"].append("--no-release-initial-path")
        with self.assertRaises(ValueError):
            make_command(Path("planner.exe"), case, 3, 1000, False)

    def test_missing_input_recorded_without_labels(self):
        case = {"id": "missing", "map": "missing.txt", "start": [0, 0, 0], "goal": [1, 1, 1]}
        result, pairs = run_case(Path("missing.exe"), case, self.root, self.root,
                                 SimpleNamespace(timeout=1, fine_timeout=1, coarse_limit=1, threads=1, export_obj=False))
        self.assertEqual(result["status"], "INPUT_ERROR")
        self.assertFalse(pairs)
        self.assertTrue((self.root / "missing/result.json").is_file())

    def test_generator_is_reproducible_and_does_not_overwrite(self):
        first = create_cases(self.root / "first", count=2, seed=42)
        second = create_cases(self.root / "second", count=2, seed=42)
        self.assertEqual(first.read_bytes(), second.read_bytes())
        for path in first.parent.glob("*.txt"):
            self.assertEqual(path.read_bytes(), (second.parent / path.name).read_bytes())
        with self.assertRaises(FileExistsError):
            create_cases(first.parent)


if __name__ == "__main__":
    unittest.main()
