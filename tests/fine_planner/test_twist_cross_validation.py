"""Independent quaternion and section-corner checks of native twist exports."""

import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from path_transformer.reference.validate_reference_geometry import validate_reference


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def distance(a, b):
    return norm([x - y for x, y in zip(a, b)])


def hamilton(a, b):
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return (aw*bw-ax*bx-ay*by-az*bz,
            aw*bx+ax*bw+ay*bz-az*by,
            aw*by-ax*bz+ay*bw+az*bx,
            aw*bz+ax*by-ay*bx+az*bw)


def quaternion_rotate(vector, axis, angle):
    half = angle / 2
    q = (math.cos(half), *(math.sin(half)*x for x in axis))
    conjugate = (q[0], -q[1], -q[2], -q[3])
    return hamilton(hamilton(q, (0, *vector)), conjugate)[1:]


def corner(center, normal, binormal, width, thickness, sign_n, sign_b):
    return tuple(center[i] + sign_n*width*normal[i]/2 +
                 sign_b*thickness*binormal[i]/2 for i in range(3))


def main(executable):
    cases = {
        "zero": (0.0, 0.5, 0.5),
        "positive90": (math.pi/2, 0.5, 0.5),
        "negative90": (-math.pi/2, 0.5, 0.5),
        "positive180": (math.pi, 0.5, 0.5),
        "positive360": (2*math.pi, 0.5, 0.5),
        "nonaxis": (1.2, 0.4, 0.9),
        "arbitraryframe": (-0.7, 0.8, 0.3),
        "short": (0.1, 0.2, 0.4),
        "smoothstep": (1.1, 0.6, 0.2),
    }
    maxima = {"position": 0.0, "orientation": 0.0, "corner": 0.0,
              "cumulative": 0.0, "orthogonality": 0.0, "density": 0.0}
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        map_hash = hashlib.sha256(b"cross-validation-empty-grid").hexdigest()
        for case, (angle, width, thickness) in cases.items():
            sampled = {}
            for spacing in (0.25, 0.125):
                output = root / f"{case}_{spacing}.json"
                subprocess.run([str(executable), str(output), map_hash,
                                "twist-case", case, str(spacing)], check=True)
                record = json.loads(output.read_text(encoding="utf-8"))
                assert validate_reference(record)
                segment = record["segments"][0]
                poses = record["reference_poses"]
                assert segment["segment_type"] == "STRAIGHT_TWIST"
                assert abs(segment["twist_angle_rad"] - angle) < 1e-12
                assert len(poses) >= 2
                start, end = segment["start"], segment["end"]
                length = distance(start, end)
                assert length > 0
                tangent = tuple((end[i] - start[i])/length for i in range(3))
                n0 = segment["twist_start_normal"]
                cumulative = 0.0
                previous_normal = None
                for index, pose in enumerate(poses):
                    u = index/(len(poses)-1)
                    profile_u = u*u*(3-2*u) if segment["twist_profile"] == "SMOOTHSTEP" else u
                    expected_angle = angle*profile_u
                    expected_center = tuple(start[i] + u*(end[i]-start[i]) for i in range(3))
                    expected_normal = quaternion_rotate(n0, tangent, expected_angle)
                    expected_binormal = cross(tangent, expected_normal)
                    maxima["position"] = max(maxima["position"],
                                             distance(pose["position"], expected_center))
                    maxima["orientation"] = max(maxima["orientation"],
                        distance(pose["tangent"], tangent),
                        distance(pose["normal"], expected_normal),
                        distance(pose["binormal"], expected_binormal))
                    frame = (pose["tangent"], pose["normal"], pose["binormal"])
                    maxima["orthogonality"] = max(maxima["orthogonality"],
                        *(abs(norm(vector)-1) for vector in frame),
                        *(abs(dot(frame[i], frame[j])) for i, j in ((0,1),(0,2),(1,2))),
                        distance(cross(frame[0], frame[1]), frame[2]))
                    for sn in (-1, 1):
                        for sb in (-1, 1):
                            maxima["corner"] = max(maxima["corner"], distance(
                                corner(pose["position"], pose["normal"], pose["binormal"],
                                       width, thickness, sn, sb),
                                corner(expected_center, expected_normal, expected_binormal,
                                       width, thickness, sn, sb)))
                    if previous_normal is not None:
                        increment = math.atan2(dot(tangent, cross(previous_normal, pose["normal"])),
                                               dot(previous_normal, pose["normal"]))
                        assert abs(increment) < math.pi/2  # No branch-cut ambiguity per sample.
                        cumulative += increment
                        assert increment*angle >= -1e-10
                    previous_normal = pose["normal"]
                maxima["cumulative"] = max(maxima["cumulative"], abs(cumulative-angle))
                assert distance(poses[-1]["normal"], segment["twist_end_normal"]) < 1e-10
                sampled[spacing] = poses
            # Compare the two densities at shared normalized positions.
            coarse, fine = sampled[0.25], sampled[0.125]
            assert (len(fine)-1) % (len(coarse)-1) == 0
            ratio = (len(fine)-1)//(len(coarse)-1)
            for i, pose in enumerate(coarse):
                finer = fine[i*ratio]
                maxima["density"] = max(maxima["density"],
                    distance(pose["position"], finer["position"]),
                    distance(pose["normal"], finer["normal"]))
        # Double precision operations with at most a few dozen samples and
        # JSON output at 17 significant digits; 1e-10 voxel is below input precision.
        for name, error in maxima.items():
            assert error < 1e-10, f"{name} error {error}"
    print(f"PASS {len(cases)} cases x 2 densities; " +
          ", ".join(f"max_{key}={value:.3e}" for key, value in maxima.items()))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
