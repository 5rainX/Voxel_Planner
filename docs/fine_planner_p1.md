# Fine Planner MVP: structure and Reference Geometry V1 contract

## Repository responsibilities

- `include/` holds stable C++ interfaces; `include/fine_planner/` defines the new geometry module.
- `src/` holds implementations; `src/fine_planner/` is independent of the legacy `Module*` directory layout.
- `tests/` holds maintained automated tests; the standalone legacy source tests and `tests/fine_planner/` are registered with CTest. Historical runs and artifacts stay in `testrecord_2026-09-27/`.
- `docs/` holds architecture and data-contract notes. Generated build and dataset outputs stay ignored.

Existing `voxel_planner`, `VoxelPlannerAPI.h`, `Module1_Main/`, `Module2_Morphology/`, `Module3_AStar/`, SDF and Dataset V1 remain in place. The new CMake `fine_planner` target links `voxel_planner`; it does not change the public coarse planner interface.

## Geometry interface and limits

`FinePlanner::plan(request, sdf)` accepts an ordered coarse path, endpoint tangent/normal constraints, positive width/thickness, safe margin, minimum bend radius and an optional ordered list of exact `SegmentSpec`s. All coordinates and radii are native voxel units; signed sweep is radians. With no segment specs, the module accepts only an exactly collinear, monotonically ordered straight coarse path. It rejects a nonstraight path instead of guessing a bend.

`SegmentSpec` supports `STRAIGHT`, `HORIZONTAL_BEND` and `VERTICAL_BEND`. A bend must include start/end, center, axis, radius and signed sweep. The type name is a caller process label; the explicit axis determines the actual plane. Samples include exact segment endpoints, analytical tangent and a section frame transported through the arc by Rodrigues rotation. Point spacing is at most 0.25 voxel (arcs are sampled more densely). Shared junctions are stored once and require position and G1 tangent continuity. Start and goal position, tangent and section normal are checked.

`GeometryValidator::validate(request, result, sdf)` is a separate module callable on an already constructed result. It checks ordered segment coverage, exact straight/arc equations, radius, sweep, frame orthonormality, endpoints, spacing and sampled SDF clearance. Clearance uses the enclosing cross-section half diagonal plus requested safe margin. This is a conservative **sampled** section check; it is not continuous full-size swept-volume proof. The export records `validation_scope=SAMPLED_SDF_CONSERVATIVE_SECTION` and `independent_full_size_collision_proven=false`.

`SqpProjector` is not treated as SQP or as a source of manufacturing features. It may later supply a position candidate, but exact process segments and independent validation remain necessary.

## Native export versus legacy sidecar

`writeReferenceGeometryV1` writes the existing `reference_geometry_v1` object with ordered `reference_centerline`, per-point `reference_poses`, process `segments`, normalized arc length, validation and provenance. Each arc segment carries its actual center, axis, radius and signed sweep. The new source has `source_kind=NATIVE_FINE_PLANNER`, a source sample/scene/map ID and **no** `raw_fine_export` or fabricated TXT hash. The legacy A06.2 sidecar retains its raw candidate TXT hash requirement. The Python validator distinguishes these source kinds, checks native ordered segment geometry and rejects invalid arc parameters.

For Dataset V1, a caller places this object in a new accepted-style JSONL row under `reference_geometry_v1`, alongside the real coarse path/poses, request dimensions, scene ID, original map path and its content SHA-256. `path_transformer.dataset.build_dataset` then uses the existing SDF, label and split pipeline. The CTest integration test demonstrates this without altering historical accepted records. A native producer must supply real map SHA-256; the older C++ DatasetGenerator's path/size hash is a different convention and is not used here.

## Build and tests

```powershell
cmake -S . -B build
cmake --build build -j 4
ctest --test-dir build --output-on-failure
python -m unittest discover -s path_transformer/tests -q
```

The original CTest baseline had two passing entries. P1-R1 registered and reproduced three previously unrun failing tests against the isolated pre-P1 HEAD. `test_macro_diversity` exposed a reuse penalty hardcoded to a one-voxel neighborhood despite `path_blocking_radius`; the penalty now follows the configured radius, while returned cost excludes search-only penalties. `test_se3_astar` expected a true arc length from the 26-neighbor coarse planner; its cost check now uses the selected lattice steps, and a separate route checks conditional pose export. `api_integration_test` required conditional poses on an open map; it now verifies the complete `centerline_poses` sequence and validates optional conditional descriptions. All 14 CTest entries pass after these targeted fixes.

The MVP does not yet fit standard process segments from arbitrary A* routes, construct straight twist segments, or prove continuous rectangular-section collision freedom. It exports locally validated supervision candidates, never global optimality claims. Dataset expansion and A07 training remain separate stages.
