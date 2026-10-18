# A06 — Path Transformer Training Dataset V1

This dataset records **supervision candidates** from the collaborator fine planner. It does not train a model. The existing `model.py` and `infer.py` retain their 13-input smoke interface; this dataset has a richer feature schema for A07 to consume deliberately.

## Build and validate

From the repository root, with Python, NumPy, and no PyTorch requirement:

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
python -m path_transformer.dataset.build_dataset --input testrecord_2026-09-27/integration_tests/results/accepted.jsonl --output dataset_outputs/a06_v1_verified
python -m path_transformer.dataset.validate_dataset dataset_outputs/a06_v1_verified
python -m unittest discover -s path_transformer/tests -v
```

`--dry-run` reads and validates without creating output. The output directory must not exist; the builder will not overwrite historical data. `dataset_outputs/` is Git-ignored. Outputs are `manifest.json`, per-sample JSON under `samples/`, shared little-endian float32 SDFs under `environments/`, group split ID lists, `stats.json`, and `rejected.jsonl`. Sample filenames are hashes of source IDs because source IDs contain `:` on Windows. The JSON schema is `schema_v1.json`; `validate_dataset.py` also checks cross-field geometry, artifact hashes, and split leakage.

## Source and provenance

The source is the immutable `testrecord_2026-09-27/integration_tests/results/accepted.jsonl`. `map_file` resolves relative to its containing directory. Before reading a map, the builder verifies that its bytes match the source `map_sha256`. The SHA-256 is the stable `map_id` and split group key; two scenes with identical map bytes share one group. `scene_id` comes from `case_id`; no ID is invented. Each sample retains source path, line, sample ID, map path/hash, and generator version. The manifest hashes the whole source JSONL. `reference.role` remains `supervision_candidate`, with the source's independent validation and optimality flags.

## Coordinate and label conventions

Coordinates, width, thickness, SDF distances, `du`, `dv`, and residuals use native **voxel units**. For each coarse point, use its `coarse_poses` tangent `t`; project its normal perpendicular to `t` and normalize to get `u`; set `v = normalize(t × u)`. This is the frame convention in `infer.py`. If a normal flips against the previous node, flip `u` and `v` together. Degenerate poses are rejected.

The fine segments are processed in source order. `STRAIGHT` and `TWIST` have a known straight centerline only when their declared length agrees with the endpoint chord. Adjacent endpoints must meet. `VERTICAL_BEND` and other curved segments are rejected as `REFERENCE_GEOMETRY_UNDERDETERMINED`: their endpoints and summary angle/length do not uniquely specify the actual centerline. No synthetic bend is inserted.

Compute cumulative arc length separately on the coarse and reconstructed reference polylines, normalize each to `[0,1]`, and linearly interpolate the reference at every coarse normalized arc position. The aligned reference and original reference polyline are both saved. For aligned position `R_i` and coarse point `P_i`, `du = (R_i-P_i)·u_i`, `dv = (R_i-P_i)·v_i`, and `longitudinal_residual = (R_i-P_i)·t_i`. Reconstruction is `P_i + du·u_i + dv·v_i`. Both endpoints receive zero `du/dv`, matching inference; endpoint mismatch above 0.5 voxel rejects the sample.

The hard reconstruction threshold is **1 voxel**, the native map resolution. A maximum error above **0.5 voxel**, its half-cell localization scale, produces `LONGITUDINAL_RESIDUAL_ABOVE_HALF_VOXEL` in `validation.quality_warnings`; above 1 voxel rejects the sample. `stats.json` reports reconstruction residual, spatial alignment offset, and longitudinal residual. Spatial alignment offset measures coarse-to-reference displacement after interpolation, not an independently known correspondence error.

For a sourced reference with valid section normals, `theta = atan2(u_ref·v, u_ref·u)` after projecting the reference normal perpendicular to the coarse tangent. Native Fine Planner twists with total signed angle within the present model's ±π/2 output range can yield valid per-node theta labels. Angles outside that range are masked. Every node covered by a native twist whose total angle exceeds ±π/2 is also masked: a bounded per-node angle cannot encode its cumulative turn, even when a 360° endpoint has the same frame as the start. Uncertain reference orientation likewise gives `theta=null` and `theta_valid=false`; no fake zero is filled. A07 must honor this per-node mask.

## SDF features

The repository's accepted records do not carry matching SDF binaries, but they do carry exact map paths and SHA-256 hashes. The builder parses those verified text occupancy maps using the same raw obstacle declaration convention as `VoxelIO.h`. It computes the 3D signed Euclidean distance field using the same X/Y/Z squared-distance-transform passes and half-voxel surface offset as `GlobalSdf` (`src/environment/SdfVolume.cpp`), then stores it once per unique accepted map. Its `z_y_x_flat` little-endian float32 layout matches `DatasetWriter.h`. `validate_dataset.py` checks each stored center query against the shared binary.

For each coarse node the five trilinear SDF queries are `P`, `P ± (width/2)u`, and `P ± (thickness/2)v`. Clearance is not added to these probe offsets. Query points outside the volume are rejected as `SDF_OUT_OF_BOUNDS`; no clamping or silent zero occurs. A map without obstacles uses the C++ `1e30` finite sentinel and is marked `environment.empty_map_sentinel=true`; A07 must mask or normalize this sentinel before feeding a model. These SDF features **do not constitute full-size collision verification**.

## Split and limitations

`split_dataset.py` assigns whole `map_sha256` groups to train, validation, or test with seed 42, recorded in the manifest. With fewer than 30 map groups the split status is `PROVISIONAL_SMALL_DATASET`. The current 13 source rows yield only a small set of reconstructible centerlines; it is **too small for a production-quality model**. Most accepted paths contain bends whose actual centerlines are unavailable, and the remaining straight examples have limited nonzero corrections. The source fine planner's acceptance is not proof of global optimality or independent geometry/manufacturing validation. Dataset V1 enables a supervised data interface; it does not prove learned path quality. A07 must explicitly choose feature normalization and handling of masked theta, SDF sentinels, and sample quality warnings before training.
