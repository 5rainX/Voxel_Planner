import copy
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from path_transformer.dataset.build_dataset import make_sample
from path_transformer.dataset.geometry import SampleError, signed_theta
from path_transformer.reference.enrich_reference import enrich_row
from path_transformer.reference.geometry import sample_explicit_arc, sample_linear
from path_transformer.reference.validate_reference_geometry import validate_reference


def straight(start, end, normal=(0, 0, 1), kind="STRAIGHT", end_normal=None):
    start, end = np.asarray(start), np.asarray(end)
    direction = ((end - start) / np.linalg.norm(end - start)).tolist()
    return {"type": kind, "start": start.tolist(), "end": end.tolist(),
            "start_normal": list(normal), "end_normal": list(end_normal or normal),
            "start_tangent": direction, "end_tangent": direction,
            "length": float(np.linalg.norm(end - start)), "angle": 0.0, "cost": 1.0}


def source_for(segments, coarse=None):
    start, end = segments[0]["start"], segments[-1]["end"]
    if coarse is None:
        coarse = [start, end]
    poses = []
    for i, _ in enumerate(coarse):
        tangent = np.asarray(coarse[min(i+1, len(coarse)-1)]) - np.asarray(coarse[max(i-1, 0)])
        tangent = (tangent / np.linalg.norm(tangent)).tolist()
        poses.append([0, 0, 1] + tangent)
    return {"sample_id": "synthetic:1", "case_id": "synthetic", "map_sha256": "a" * 64,
            "map_file": "case/input.txt", "artifacts": "case", "candidate_id": 1, "fine_segments": segments,
            "coarse_centerline": coarse, "coarse_poses": poses,
            "request": {"start": start, "goal": end,
                        "args": ["--width", "0.2", "--thickness", "0.2"]},
            "independent_geometry_validation": False, "global_optimality_proven": False}


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source_path = self.root / "accepted.jsonl"
        self.source_path.write_text("fixture\n", encoding="utf-8")
        raw = self.root / "case" / "busbar_path_candidate_001.txt"
        raw.parent.mkdir()
        raw.write_text("synthetic fixture\n", encoding="utf-8")

    def tearDown(self):
        self.temp.cleanup()

    def enrich(self, source):
        return enrich_row(source, self.source_path, 1)

    def test_straight_dense_sampling(self):
        points = sample_linear([0, 0, 0], [1, 0, 0])
        self.assertEqual(len(points), 5)
        self.assertLessEqual(np.linalg.norm(np.diff(points, axis=0), axis=1).max(), 0.25)
        np.testing.assert_allclose(points[0], [0, 0, 0])
        np.testing.assert_allclose(points[-1], [1, 0, 0])

    def test_arc_dense_sampling(self):
        points = sample_explicit_arc([1, 0, 0], [0, 0, 0], [0, 0, 1], 1, math.pi / 2)
        self.assertLessEqual(np.linalg.norm(np.diff(points, axis=0), axis=1).max(), 0.25)
        np.testing.assert_allclose(np.linalg.norm(points, axis=1), 1, atol=1e-8)
        np.testing.assert_allclose(points[-1], [0, 1, 0], atol=1e-8)

    def test_segment_join_no_duplicate(self):
        source = source_for([straight([0, 1, 1], [1, 1, 1]),
                             straight([1, 1, 1], [2, 1, 1])])
        ref = self.enrich(source)
        self.assertEqual(len(ref["reference_centerline"]), 9)
        self.assertTrue(validate_reference(ref, source))
        self.assertEqual(ref["validation"]["max_junction_position_error_voxel"], 0)

    def test_segment_join_continuity(self):
        source = source_for([straight([0, 1, 1], [1, 1, 1]),
                             straight([1.1, 1, 1], [2.1, 1, 1])])
        with self.assertRaisesRegex(SampleError, "segment 1 gap"):
            self.enrich(source)

    def test_vertical_bend_requires_explicit_parameters(self):
        bend = {"type": "VERTICAL_BEND", "start": [2, 1, 1], "end": [1, 2, 1],
                "start_tangent": [0, 1, 0], "end_tangent": [-1, 0, 0],
                "start_normal": [0, 0, 1], "end_normal": [0, 0, 1],
                "length": math.pi / 2, "angle": math.pi / 2, "cost": 1}
        with self.assertRaises(SampleError):
            self.enrich(source_for([bend]))
        bend["exact_arc_parameters"] = {"source": "DIRECT_PLANNER_ARC_PARAMETERS",
                                         "center": [1, 1, 1], "axis": [0, 0, 1],
                                         "radius": 1, "sweep_angle_rad": math.pi / 2}
        source = source_for([bend], [[2, 1, 1], [1.5, 1.5, 1], [1, 2, 1]])
        ref = self.enrich(source)
        self.assertTrue(validate_reference(ref, source))
        self.assertEqual(ref["segments"][0]["geometry_reconstruction_method"], "EXACT_PLANNER_ARC_PARAMETERS")
        self.assertLess(ref["validation"]["max_spacing_voxel"], 0.25)

    def test_orientation_frame_orthogonal_and_matches_tangent(self):
        source = source_for([straight([0, 1, 1], [2, 1, 1])])
        ref = self.enrich(source)
        for pose in ref["reference_poses"]:
            t, u, v = (np.asarray(pose[k]) for k in ("tangent", "normal", "binormal"))
            self.assertAlmostEqual(np.dot(t, u), 0)
            self.assertAlmostEqual(np.dot(u, v), 0)
            np.testing.assert_allclose(np.cross(t, u), v)
            np.testing.assert_allclose(t, [1, 0, 0])

    def test_orientation_continuity_without_twist(self):
        source = source_for([straight([0, 1, 1], [1, 1, 1]),
                             straight([1, 1, 1], [2, 1, 1])])
        ref = self.enrich(source)
        normals = [p["normal"] for p in ref["reference_poses"]]
        self.assertTrue(all(np.dot(a, b) > 0.99 for a, b in zip(normals, normals[1:])))

    def test_known_twist_theta_and_mask(self):
        source = source_for([straight([0, 1, 1], [1, 1, 1]),
                             straight([1, 1, 1], [2, 1, 1], kind="TWIST", end_normal=(0, -1, 0)),
                             straight([2, 1, 1], [3, 1, 1], normal=(0, -1, 0))],
                            [[0, 1, 1], [1, 1, 1], [2, 1, 1], [3, 1, 1]])
        ref = self.enrich(source)
        self.assertTrue(any(not p["orientation_valid"] for p in ref["reference_poses"]))
        self.assertAlmostEqual(signed_theta((np.array([1,0,0]), np.array([0,0,1]), np.array([0,-1,0])),
                                            ref["reference_poses"][-1]["normal"]), math.pi/2)

    def test_angle_wraparound(self):
        frame = (np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1]))
        self.assertAlmostEqual(abs(signed_theta(frame, [0,-1,1e-8])), math.pi, places=6)

    def test_enriched_reference_preferred_and_vertical_bend_builds(self):
        bend = {"type": "VERTICAL_BEND", "start": [2, 1, 1], "end": [1, 2, 1],
                "start_tangent": [0, 1, 0], "end_tangent": [-1, 0, 0],
                "start_normal": [0, 0, 1], "end_normal": [0, 0, 1],
                "length": math.pi/2, "angle": math.pi/2, "cost": 1,
                "exact_arc_parameters": {"source": "DIRECT_PLANNER_ARC_PARAMETERS",
                                         "center": [1, 1, 1], "axis": [0, 0, 1],
                                         "radius": 1, "sweep_angle_rad": math.pi/2}}
        source = source_for([bend], [[2, 1, 1], [1.5, 1.5, 1], [1, 2, 1]])
        source["reference_geometry_v1"] = self.enrich(source)
        field = np.full((4, 4, 4), 2.0, dtype="<f4")
        sample = make_sample(source, self.source_path, 1, field, "environments/x.bin")
        self.assertEqual(sample["reference"]["reference_geometry_version"], "reference_geometry_v1")
        self.assertGreater(len(sample["reference"]["centerline"]), 2)

    def test_theta_valid_from_reference_pose(self):
        source = source_for([straight([0, 1, 1], [2, 1, 1], normal=(0, -1, 0))],
                            [[0, 1, 1], [1, 1, 1], [2, 1, 1]])
        source["reference_geometry_v1"] = self.enrich(source)
        sample = make_sample(source, self.source_path, 1, np.full((4,4,4), 2.0, dtype="<f4"), "environments/x.bin")
        self.assertTrue(all(t["theta_valid"] for t in sample["targets"]))
        self.assertTrue(all(abs(t["theta"] - math.pi/2) < 1e-6 for t in sample["targets"]))

    def test_theta_mask_when_orientation_missing(self):
        source = source_for([straight([0,1,1],[1,1,1],kind="TWIST",end_normal=(0,-1,0))],
                            [[0,1,1],[0.5,1,1],[1,1,1]])
        source["reference_geometry_v1"] = self.enrich(source)
        sample = make_sample(source,self.source_path,1,np.full((4,4,4),2.0,dtype="<f4"),"environments/x.bin")
        self.assertIsNone(sample["targets"][1]["theta"])
        self.assertFalse(sample["targets"][1]["theta_valid"])


if __name__ == "__main__":
    unittest.main()
