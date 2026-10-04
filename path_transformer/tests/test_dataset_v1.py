import copy
import hashlib
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from path_transformer.dataset.geometry import (
    SampleError, align_paths, frames, reference_polyline, signed_theta, target,
)
from path_transformer.dataset.sdf_features import (
    features, field_from_occupancy, query, read_sdf_binary, read_text_map,
)
from path_transformer.dataset.split_dataset import split_samples
from path_transformer.dataset.validate_dataset import validate_sample


def valid_sample():
    node = lambda x, s: {"position": [x, 1, 1], "tangent": [1, 0, 0],
                         "normal": [0, 1, 0], "binormal": [0, 0, 1],
                         "s_normalized": s,
                         "sdf": dict.fromkeys(("center", "u_pos", "u_neg", "v_pos", "v_neg"), 0.5)}
    label = {"du": 0.0, "dv": 0.0, "theta": None,
             "theta_valid": False, "longitudinal_residual": 0.0}
    return {
        "schema_version": "path_transformer_dataset_v1", "sample_id": "case:1",
        "provenance": {"source": "collaborator_fine_planner", "source_file": "source.jsonl",
                       "source_line": 1, "source_sample_id": "case:1", "scene_id": "case",
                       "map_id": "hash", "map_file": "map.txt", "generator_version": "A06-v1",
                       "reference_type": "supervision_candidate"},
        "environment": {"map_sha256": "hash", "sdf_artifact": "environments/hash.bin",
                        "sdf_artifact_sha256": "0" * 64, "empty_map_sentinel": False,
                        "sdf_shape_zyx": [3, 3, 3], "origin": [0, 0, 0],
                        "voxel_size": [1, 1, 1], "layout": "z_y_x_flat",
                        "dtype": "float32_le", "sdf_valid": True},
        "busbar": {"width": 1, "thickness": 1, "units": "voxel"},
        "start": {"position": [0, 1, 1]}, "goal": {"position": [1, 1, 1]},
        "coarse_path": [node(0, 0), node(1, 1)],
        "reference": {"source": "collaborator_fine_planner", "role": "supervision_candidate",
                      "independent_geometry_validation": False, "global_optimality_proven": False,
                      "centerline": [[0, 1, 1], [1, 1, 1]],
                      "aligned_centerline": [[0, 1, 1], [1, 1, 1]],
                      "source_segment_types": ["STRAIGHT"]},
        "targets": [copy.deepcopy(label), copy.deepcopy(label)],
        "validation": {"geometry_label_valid": True, "full_feature_valid": True,
                       "quality_threshold_voxel": 1.0}}


class AlignmentTests(unittest.TestCase):
    def test_identical_path_alignment(self):
        path = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]
        aligned, s, _, _ = align_paths(path, path)
        np.testing.assert_allclose(aligned, path)
        np.testing.assert_allclose(s, [0, 0.5, 1])

    def test_different_sampling_density(self):
        aligned, _, _, _ = align_paths([[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                                       [[0, 0, 0], [2, 0, 0]])
        np.testing.assert_allclose(aligned[1], [1, 0, 0])

    def test_uniform_u_offset(self):
        coarse = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]
        aligned, _, _, _ = align_paths(coarse, [[0, 2, 0], [2, 2, 0]])
        frame = frames([[0, 1, 0, 1, 0, 0]])[0]
        self.assertTrue(all(abs(target(np.asarray(p), r, frame)[0] - 2) < 1e-8
                            for p, r in zip(coarse, aligned)))

    def test_nonuniform_reference_sampling(self):
        aligned, _, _, _ = align_paths([[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                                       [[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0], [2, 0, 0]])
        np.testing.assert_allclose(aligned[1], [1, 0, 0])

    def test_repeated_points(self):
        with self.assertRaises(SampleError):
            align_paths([[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [1, 0, 0]])
        aligned, _, _, _ = align_paths([[0, 0, 0], [1, 0, 0]],
                                       [[0, 0, 0], [0, 0, 0], [1, 0, 0]])
        np.testing.assert_allclose(aligned[-1], [1, 0, 0])

    def test_bend_without_centerline_is_rejected(self):
        with self.assertRaisesRegex(SampleError, "VERTICAL_BEND"):
            reference_polyline([{"type": "VERTICAL_BEND", "start": [0, 0, 0],
                                 "end": [1, 0, 0], "length": 1.2}])


class GeometryTests(unittest.TestCase):
    def test_local_frame_orthogonal(self):
        t, u, v = frames([[0, 2, 1, 1, 1, 0]])[0]
        self.assertAlmostEqual(np.dot(t, u), 0)
        self.assertAlmostEqual(np.dot(t, v), 0)
        self.assertAlmostEqual(np.dot(u, v), 0)
        np.testing.assert_allclose(np.cross(t, u), v)

    def test_frame_continuity(self):
        f = frames([[0, 1, 0, 1, 0, 0], [0, -1, 0, 1, 0, 0]])
        self.assertGreater(np.dot(f[0][1], f[1][1]), 0.99)

    def test_du_dv_reconstruction(self):
        frame = frames([[0, 1, 0, 1, 0, 0]])[0]
        du, dv, longitudinal, error = target(np.array([0, 0, 0]), np.array([0, 2, 3]), frame)
        self.assertAlmostEqual(du, 2)
        self.assertAlmostEqual(dv, 3)
        self.assertAlmostEqual(longitudinal, 0)
        self.assertAlmostEqual(error, 0)

    def test_theta_signed_rotation_and_wraparound(self):
        frame = frames([[0, 1, 0, 1, 0, 0]])[0]
        self.assertAlmostEqual(signed_theta(frame, [0, 1, 0]), 0)
        self.assertAlmostEqual(signed_theta(frame, [0, 0, 1]), math.pi / 2)
        self.assertAlmostEqual(abs(signed_theta(frame, [0, -1, 1e-8])), math.pi, places=6)

    def test_theta_missing_mask(self):
        sample = valid_sample()
        self.assertFalse(sample["targets"][0]["theta_valid"])
        self.assertIsNone(sample["targets"][0]["theta"])
        self.assertTrue(validate_sample(sample))


class SdfTests(unittest.TestCase):
    def setUp(self):
        occupancy = np.zeros((3, 3, 3), dtype=bool)
        occupancy[1, 1, 1] = True
        self.field = field_from_occupancy(occupancy)

    def test_sdf_center_query(self):
        self.assertAlmostEqual(query(self.field, [1, 1, 1]), -0.5)
        self.assertAlmostEqual(query(self.field, [2, 1, 1]), 0.5)
        self.assertAlmostEqual(query(self.field, [0, 0, 0]), math.sqrt(3) - 0.5, places=6)

    def test_sdf_offset_queries(self):
        values = features(self.field, [1, 1, 1], [1, 0, 0], [0, 1, 0], 2, 2)
        self.assertAlmostEqual(values["center"], -0.5)
        self.assertAlmostEqual(values["u_pos"], 0.5)
        self.assertAlmostEqual(values["v_neg"], 0.5)

    def test_sdf_out_of_bounds(self):
        with self.assertRaises(SampleError):
            query(self.field, [-0.1, 1, 1])

    def test_sdf_binary_roundtrip_and_map_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "map.txt"
            raw = b"map: 3 3 3\n1 1 1 0\n"
            path.write_bytes(raw)
            occupancy = read_text_map(path, hashlib.sha256(raw).hexdigest())
            np.testing.assert_array_equal(occupancy, self.field < 0)
            binary = Path(directory) / "sdf.bin"
            self.field.astype("<f4").tofile(binary)
            np.testing.assert_allclose(read_sdf_binary(binary, (3, 3, 3)), self.field)
            with self.assertRaises(SampleError):
                read_text_map(path, "wrong-hash")

    def test_sdf_matches_brute_force_distances(self):
        occupancy = np.zeros((4, 4, 5), dtype=bool)
        occupancy[0, 1, 3] = True
        occupancy[3, 3, 0] = True
        field = field_from_occupancy(occupancy)
        obstacles = np.argwhere(occupancy)
        free = np.argwhere(~occupancy)
        for z in range(4):
            for y in range(4):
                for x in range(5):
                    point = np.array([z, y, x])
                    candidates = free if occupancy[z, y, x] else obstacles
                    expected = min(np.linalg.norm(point - candidate) for candidate in candidates) - 0.5
                    if occupancy[z, y, x]:
                        expected = -expected
                    self.assertAlmostEqual(field[z, y, x], expected, places=5)


class SchemaSplitTests(unittest.TestCase):
    def test_valid_sample_passes(self):
        self.assertTrue(validate_sample(valid_sample()))

    def test_missing_field_fails(self):
        sample = valid_sample()
        del sample["provenance"]["map_id"]
        with self.assertRaises((KeyError, SampleError)):
            validate_sample(sample)

    def test_nonfinite_value_fails(self):
        sample = valid_sample()
        sample["coarse_path"][0]["position"][0] = float("nan")
        with self.assertRaises(SampleError):
            validate_sample(sample)

    def test_target_length_mismatch_fails(self):
        sample = valid_sample()
        sample["targets"].pop()
        with self.assertRaises(SampleError):
            validate_sample(sample)

    def test_same_map_never_crosses_splits(self):
        samples = [{"sample_id": f"s{i}", "provenance": {"map_id": f"map{i // 2}"}}
                   for i in range(10)]
        splits, _ = split_samples(samples)
        lookup = {sample_id: name for name, ids in splits.items() for sample_id in ids}
        self.assertTrue(all(lookup[f"s{i}"] == lookup[f"s{i + 1}"] for i in range(0, 10, 2)))

    def test_split_is_deterministic(self):
        samples = [{"sample_id": str(i), "provenance": {"map_id": str(i // 2)}} for i in range(10)]
        self.assertEqual(split_samples(samples, seed=42), split_samples(list(reversed(samples)), seed=42))


if __name__ == "__main__":
    unittest.main()
