"""Run Path Transformer on Voxel_Planner integration-test JSONL."""

import argparse
import json
from pathlib import Path

import torch

from model import PathTransformer


def unit(vector):
    length = torch.linalg.vector_norm(vector)
    if not torch.isfinite(length) or length < 1e-8:
        raise ValueError("zero or non-finite pose vector")
    return vector / length


def features_and_frame(sample):
    points = torch.tensor(sample["coarse_centerline"], dtype=torch.float32)
    poses = torch.tensor(sample["coarse_poses"], dtype=torch.float32)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] < 2:
        raise ValueError("coarse_centerline must contain at least two xyz points")
    if poses.shape != (points.shape[0], 6):
        raise ValueError("coarse_poses must contain one [nx,ny,nz,tx,ty,tz] per point")
    if not torch.isfinite(points).all() or not torch.isfinite(poses).all():
        raise ValueError("path and poses must be finite")

    normals = []
    binormals = []
    tangents = []
    for pose in poses:
        tangent = unit(pose[3:])
        raw_normal = pose[:3] - torch.dot(pose[:3], tangent) * tangent
        normal = unit(raw_normal)
        tangents.append(tangent)
        normals.append(normal)
        binormals.append(unit(torch.linalg.cross(tangent, normal)))
    tangents = torch.stack(tangents)
    normals = torch.stack(normals)
    binormals = torch.stack(binormals)

    segment_lengths = torch.linalg.vector_norm(points[1:] - points[:-1], dim=1)
    arc = torch.cat((torch.zeros(1), torch.cumsum(segment_lengths, dim=0)))
    if arc[-1] <= 0:
        raise ValueError("coarse_centerline has zero length")
    scale = max(float(arc[-1]), 1.0)
    request_args = sample.get("request", {}).get("args", [])

    def option(name):
        try:
            return float(request_args[request_args.index(name) + 1]) / scale
        except (ValueError, IndexError):
            raise ValueError(f"missing required request option: {name}") from None

    count = points.shape[0]
    widths = torch.full((count, 1), option("--width"))
    thicknesses = torch.full((count, 1), option("--thickness"))
    # The current integration export has no SDF samples. Zero marks unavailable.
    clearances = torch.zeros((count, 1))
    features = torch.cat((
        (points - points[0]) / scale, tangents, normals,
        (arc / scale)[:, None], widths, thicknesses, clearances,
    ), dim=1)
    return features[None], points, normals, binormals


def predict(sample, model, status):
    features, points, normals, binormals = features_and_frame(sample)
    with torch.inference_mode():
        inferred = model(features)[0]
    corrections = inferred.clone()
    corrections[0] = 0.0
    corrections[-1] = 0.0
    corrected = (points + corrections[:, :1] * normals
                 + corrections[:, 1:2] * binormals)
    return {
        "sample_id": sample["sample_id"],
        "status": status,
        "units": "voxel",
        "clearance_feature_available": False,
        "coarse_centerline": points.tolist(),
        "predicted_centerline": corrected.tolist(),
        "local_corrections_u_v_theta": corrections.tolist(),
        "hard_constraints_verified": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="accepted.jsonl from the Voxel_Planner integration test")
    parser.add_argument("--output", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--checkpoint", type=Path,
                      help="trained checkpoint with architecture and state_dict")
    mode.add_argument("--smoke", action="store_true",
                      help="untrained interface test; outputs must not be used as plans")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; choose a new file")

    model = PathTransformer()
    if args.checkpoint:
        checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=True)
        if checkpoint.get("architecture") != "PathTransformer-v1" or not checkpoint.get("trained"):
            parser.error("checkpoint must be a trained PathTransformer-v1 state_dict")
        model.load_state_dict(checkpoint["state_dict"], strict=True)
        status = "CHECKPOINT_INFERENCE_UNVALIDATED"
    else:
        status = "UNTRAINED_SMOKE_ONLY"
    model.eval()

    rows = []
    with args.input.open("r", encoding="utf-8") as source:
        for line_number, line in enumerate(source, 1):
            if line.strip():
                try:
                    rows.append(predict(json.loads(line), model, status))
                except (KeyError, ValueError, TypeError) as error:
                    parser.error(f"input line {line_number}: {error}")
    if not rows:
        parser.error("input has no samples")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as target:
        for row in rows:
            target.write(json.dumps(row, ensure_ascii=False) + "\n")
    print(f"{status}: {len(rows)} samples -> {args.output}")


if __name__ == "__main__":
    main()
