# A06.1 — Reference Geometry Enrichment

This module creates `reference_geometry_v1` alongside each immutable historical `accepted.jsonl` row. It does not alter the fine planner or the historical files.

## Source and information loss

`testrecord_2026-09-27/integration_tests/run_tests.py` invokes `hybrid_validation_cli.exe`, built by the integration CMake target from a **different repository's** `src/main.cpp`, `busbar_segment.cpp`, `segment_generator.cpp`, and `busbar_planner.cpp`. The historical `provenance.json` names its source as `E:\Learning\Subjects\路径生成\busbar_hybrid_state\busbar_hybrid_state`; that path and the executable are unavailable on this machine. The run produces candidate TXT, CSV and OBJ files. `read_segments()` parses the candidate TXT's type, endpoints, endpoint normals/tangents, angle, length and cost, then places only those summaries in `accepted.jsonl`. The committed OBJ says `Occupied-Voxel Shell Surface`; its mesh is a rasterized occupied volume, not an ordered centerline or per-point section frame.

The bend information may have existed inside the external planner, but the available source cannot prove this. It is lost **no later than its committed TXT/OBJ export**. The historical `VERTICAL_BEND` rows omit bend center, radius, plane axis, sweep direction and dense curve. Some exported bend `length` values are even less than the endpoint chord (for example `single_obstacle:1` has length 1.8326 and chord √5 ≈ 2.2361). A continuous arc cannot have that length, so the summary cannot be interpreted as exact arc parameters. The eight bend cases remain `UNRECOVERABLE` with a per-sample rejection record. The future external fine-planner export should provide either the actual ordered sampled centerline and section frame, or exact center/radius/axis/signed sweep plus the planner's orientation transport law. The adjacent `busbar_path_candidate_*.txt` source and SHA-256 remain in enriched provenance.

## Supported geometry and orientation

For a `STRAIGHT` or `TWIST` segment, the chord is the exact centerline only when the exported length equals the chord and the endpoint tangents agree with its direction. We sample it at **at most 0.25 voxel**, using the existing C++ `SqpProjectorOptions.validation_sample_voxel_fraction = 0.25` on unit voxels. Segment endpoints are kept; identical junction endpoints are stored once. Position gaps above 1e-4 voxel or tangent discontinuity above 5° reject the record. The output records maximum spacing, junction errors, normalized cumulative arc length, and segment/process ownership for every node.

For non-twisting straight segments whose exported endpoint normals agree, the exported frame is held constant and orthogonalized against the geometric tangent. A `TWIST` retains the exported start and end orientation, but its **interior orientation is invalid** because the export lacks the planner's twist angle profile and direction. No shortest-arc interpolation is claimed to reproduce an unknown physical twist. `orientation_valid` is per point; invalid nodes have `normal=null`, `binormal=null`. At joints, incompatible endpoint normals invalidate the shared node. A real 180° twist is never removed by a generic sign-flip heuristic. Dataset theta uses the existing A06 `signed_theta` definition relative to the coarse `t/u/v` frame, and remains masked where orientation is unavailable.

`geometry.py` includes a deterministic exact-circle sampler for a **future explicit** `exact_arc_parameters` source (`center`, `radius`, `axis`, signed `sweep_angle_rad`, `source=DIRECT_PLANNER_ARC_PARAMETERS`). It verifies radius, plane, endpoint and length. Historical bend records have no such parameters, so this path is never used to invent their geometry.

## A06.2 direct planner export intake

The historical external fine planner source and executable are not present here, so this repository cannot add an exporter to that planner or recover the eight historical bends without its real output. P1's native Fine Planner has a separate Reference Geometry V1 export and does not claim to reproduce those historical runs. This optional importer is only for external planner output tied to an existing TXT. Place its JSON sidecar under a separate export root at the same relative case path as the candidate TXT, replacing `.txt` with `.reference.json` (for example `case/busbar_path_candidate_001.reference.json`). Run enrichment with `--planner-export-root <directory>`. Do not edit the historical TXT or `accepted.jsonl`.

```json
{
  "export_version": "fine_planner_reference_v1",
  "sample_id": "case:1",
  "map_sha256": "<source map SHA-256>",
  "raw_fine_export_sha256": "<candidate TXT SHA-256>",
  "segments": [
    {
      "segment_index": 0,
      "segment_type": "STRAIGHT",
      "poses": [
        {"position": [0, 0, 0], "tangent": [1, 0, 0],
         "normal": [0, 0, 1], "binormal": [0, -1, 0],
         "orientation_valid": true},
        {"position": [0.25, 0, 0], "tangent": [1, 0, 0],
         "normal": [0, 0, 1], "binormal": [0, -1, 0],
         "orientation_valid": true}
      ]
    }
  ]
}
```

The example shows field names; coordinates must match the corresponding source segment. Export one segment entry per original `fine_segments` entry, with both endpoints and at least two ordered poses. Positions and frames must come from the planner, with sample spacing at most 0.25 voxel. If a point's physical orientation is unavailable, set `orientation_valid=false` and both `normal` and `binormal` to `null`; theta will be masked. The importer checks sample ID, map hash, raw TXT hash, segment identity, endpoints, spacing, tangent agreement and frame orthogonality before Dataset V1 uses it. These checks bind and validate the supplied data; they cannot prove that the external planner generated it. A sidecar with a mismatched binding is rejected, never replaced with guessed geometry. An absent sidecar leaves the conservative A06.1 reconstruction and bend rejection in place. The enriched record stores the sidecar path and SHA-256 with `source_kind=LEGACY_TXT_PLANNER_SIDECAR`; P1 native exports use `source_kind=NATIVE_FINE_PLANNER` and have no historical TXT hash.

## Commands

From the repository root:

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
python -m path_transformer.reference.enrich_reference --input testrecord_2026-09-27/integration_tests/results/accepted.jsonl --output dataset_outputs/a06_1_v1_release/reference_geometry
python -m path_transformer.reference.validate_reference_geometry dataset_outputs/a06_1_v1_release/reference_geometry/accepted_enriched.jsonl
python -m path_transformer.dataset.build_dataset --input dataset_outputs/a06_1_v1_release/reference_geometry/accepted_enriched.jsonl --map-root testrecord_2026-09-27/integration_tests/results --output dataset_outputs/a06_1_v1_release/dataset
python -m path_transformer.dataset.validate_dataset dataset_outputs/a06_1_v1_release/dataset
python -m path_transformer.reference.compare_datasets --a06 dataset_outputs/a06_v1_release --a06-1 dataset_outputs/a06_1_v1_release
```

The enrichment input is streamed line by line and its output directory must not already exist. The derived `accepted_enriched.jsonl` retains each source row plus a `reference_geometry_v1` field. Valid rows have dense geometry; unrecoverable rows have a status and reason. The A06 builder prioritizes that field over legacy `fine_segments`; it never falls back after an explicit enriched rejection. `--map-root` points to the original hashed map artifacts, since the enriched JSONL lives in a new directory. The old A06 input without `reference_geometry_v1` remains supported.

The output is a *candidate* supervision representation. `independent_geometry_validation=false` and `global_optimality_proven=false` remain unchanged. Dense geometry and SDF features do not certify manufacturing safety. The current few maps and label distribution do not support production training.
