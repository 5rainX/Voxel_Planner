# P2.1 source comparison and independent numerical validation

## Fixed source reviewed

The source review preceded the additional P2.1 cross-validation changes. The
commit IDs below came from read-only `git ls-remote` checks. No third-party code
was copied into this repository, built, or added as a dependency.

| Project and fixed revision | Source actually inspected | License | Relevant observation and reuse boundary |
| --- | --- | --- | --- |
| [Eigen](https://github.com/eigen-mirror/eigen/tree/3c3bf4532bf9e059e3999827b00b7452c2bcde7d), `3c3bf4532bf9e059e3999827b00b7452c2bcde7d` | [`AngleAxis::toRotationMatrix`](https://github.com/eigen-mirror/eigen/blob/3c3bf4532bf9e059e3999827b00b7452c2bcde7d/Eigen/src/Geometry/AngleAxis.h#L242-L276), [`QuaternionBase::_transformVector`, angle-axis assignment and `slerp`](https://github.com/eigen-mirror/eigen/blob/3c3bf4532bf9e059e3999827b00b7452c2bcde7d/Eigen/src/Geometry/Quaternion.h) | MPL-2.0 in the inspected source headers | Angle-axis rotation requires a unit axis. Quaternion vector rotation provides a separate numerical route. Quaternion interpolation selects an equivalent short path; endpoint orientation alone does not encode a signed full turn. We used the mathematics for an independent test, not Eigen as a runtime dependency. |
| [ParallelTransportFrame](https://github.com/vicrucann/ParallelTransportFrame/tree/3357b15725a3d953ae740d069d44a2595d6a2d81), `3357b15725a3d953ae740d069d44a2595d6a2d81` | [`PTFTube::extractTangents`, `firstFrame`, `nextFrame`](https://github.com/vicrucann/ParallelTransportFrame/blob/3357b15725a3d953ae740d069d44a2595d6a2d81/libPTFTube/PTFTube.cpp#L256-L347) | [MIT](https://github.com/vicrucann/ParallelTransportFrame/blob/3357b15725a3d953ae740d069d44a2595d6a2d81/LICENSE) | Rotates a frame to follow changing tangents; its collinear fallback chooses an arbitrary initial direction. On a straight centerline, tangent transport introduces no active twist. Its transport principle is useful for a future curved baseline frame, but neither its arbitrary initial frame nor its OpenSceneGraph dependency meets this explicit-section-orientation contract. |
| [OCCT](https://github.com/Open-Cascade-SAS/OCCT/tree/3d097a0328e71b826377d4814ab05ec3c3d23871), `3d097a0328e71b826377d4814ab05ec3c3d23871` | [`BRepOffsetAPI_MakePipe` constructor](https://github.com/Open-Cascade-SAS/OCCT/blob/3d097a0328e71b826377d4814ab05ec3c3d23871/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MakePipe.cxx), [`BRepOffsetAPI_MakePipeShell::SetMode` and profile setup](https://github.com/Open-Cascade-SAS/OCCT/blob/3d097a0328e71b826377d4814ab05ec3c3d23871/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MakePipeShell.cxx#L793-L845) | [LGPL-2.1 with OCCT exception](https://github.com/Open-Cascade-SAS/OCCT/blob/3d097a0328e71b826377d4814ab05ec3c3d23871/OCCT_LGPL_EXCEPTION.txt), or commercial license | A general CAD sweep accepts spine and profile and has frame-mode choices. It may be useful for later CAD output or swept-volume studies. Its API does not itself provide this planner's signed cumulative twist, process-rate limits, or SDF validation. No OCCT dependency was introduced. |

## Selected geometry

The repository uses a right-handed frame with unit tangent `t`, section normal
`n`, binormal `b = t × n`, radian angles, and positive rotation by the right-hand
rule around the start-to-end tangent. `n0` is the *initial section direction*,
not a world-axis angle. On a straight segment with `u` the normalized arc
length, `theta(u) = delta_theta × f(u)`, where `f(u)=u` for `LINEAR` and
`f(u)=3u²−2u³` for `SMOOTHSTEP`. Thus `n(u)=cos(theta)n0+sin(theta)b0` and
`b(u)=−sin(theta)n0+cos(theta)b0`. The center remains
`start + u × (end−start)`. The explicit signed `delta_theta` remains in the
segment and native Reference V1 export, even when the final frame equals the
initial frame after 360°. This is relative rotation about `t`; an absolute
world orientation is fixed by the supplied initial frame.

The native implementation uses a Rodrigues rotation and orthonormalizes the
resulting frame. The independent test computes a scalar Hamilton quaternion
product `q × (0,v) × conjugate(q)` directly. It does not call the native
rotation helper or reconstruct the answer from exported intermediate normals.
For each sample it compares the center, `t/n/b`, and four rectangle corners
`center ± (width/2)n ± (thickness/2)b`. It also accumulates *signed* adjacent
frame rotations using `atan2(t·(n_prev × n_next), n_prev·n_next)`; this checks
360° without reducing it to the endpoint orientation. The test checks frame
unit length, orthogonality, handedness, endpoint normal, direction of each
increment, and shared samples across two densities.

## Reproduction and measured result

From the repository root on the current Windows build configuration:

```powershell
cmake --build build -j 4
python tests/fine_planner/test_twist_cross_validation.py bin/test_fine_planner.exe
ctest --test-dir build --output-on-failure
python -m unittest discover -s path_transformer/tests -q
```

The independent test covers 0°, +90°, −90°, +180°, +360°, a non-axis-aligned
3D line, an arbitrary legal starting normal, a 0.25-voxel legal segment, and
smoothstep; widths and thicknesses vary. Each runs at 0.25- and 0.125-voxel
maximum sample spacing. The C++ test separately rejects zero-length twists.
Its measured maxima over all 18 exports were:

| Error | Maximum |
| --- | ---: |
| Center position | `0.000e+00` voxel |
| `t/n/b` vector | `5.579e-16` |
| Four corner positions | `8.882e-16` voxel |
| Signed cumulative angle | `2.220e-16` rad |
| Unit length, orthogonality and handedness | `2.719e-16` |
| Coincident samples at two densities | `0.000e+00` voxel |

The pass threshold is `1e-10` in the corresponding units. It allows double
precision arithmetic and 17-digit JSON serialization while remaining far
below the existing `1e-4`-voxel Reference geometry validation tolerance; it
was set before measuring the final maxima. During test development, a new
quarter-point assertion failed on the +180° case because the integer sample
index was treated as exactly `u=0.25`. The assertion now uses its actual
`index/(count−1)` parameter. This was a test error, not a geometry change.

## Decision and remaining scope

No Eigen, OpenSceneGraph, or OCCT dependency was introduced. The test contains
an independently written scalar quaternion calculation informed by Eigen's
documented/source mathematics, without copying its implementation. General
parallel transport and CAD sweeping do not replace explicit signed section
rotation. The current SDF check remains conservative at sampled centerline
points using a section half-diagonal; it is not a continuous full swept-volume
collision proof. No coupled bend–twist manufacturing model is claimed.
