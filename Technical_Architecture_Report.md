# Technical Architecture Report

## 1. System Architecture Overview

Voxel Planner is a layered, pose-aware motion-planning system. The public
facade in `include/VoxelPlannerAPI.h` hides the implementation modules behind
two operations: `loadMap()` creates a reusable `ProcessedMap`, and
`findPaths()` performs one bounded multi-route planning request.

The execution pipeline is:

```text
Map Ingestion
    -> VoxelGrid construction and occupancy-block indexing
Pre-processing (Pose Footprints)
    -> discrete SE(3) pose table
    -> lazy/static pose-mask population
SE(3) Coarse A* Search
    -> state lattice over (voxel, poseId)
    -> collision-validated straight, bend, and twist actions
Macro-Diversity Iteration
    -> accept/reject by overlap ratio
    -> soft penalties on accepted centerlines and neighbors
API Output
    -> endpoint restoration, public PathResult conversion, cost sorting
```

### Map ingestion and preprocessing

`VoxelIO` accepts the repository's text format (`header x y z` followed by
voxel records) and ASCII structured VTK volumes. Raw occupancy is stored in a
dense, X-fastest `VoxelGrid` array. An 8 x 8 x 8 occupancy-block index gives a
cheap axis-aligned rejection test before detailed footprint or sweep checks.

`loadMap()` validates the physical dimensions, generates the discrete pose
table, and attaches pose-footprint data to the grid. Pose masks are sparse and
lazily evaluated: open voxels that allow every pose use an all-poses bit mask;
voxels with a strict subset of valid poses use packed 64-bit words in a sparse
slot table. This keeps the common open-space representation compact while
retaining exact per-pose collision information near boundaries and obstacles.

### Search and output

`CoarseAStar` builds an action catalog once per request. It caches static
straight sweeps and lazily caches bend/twist sweeps. Each A* expansion first
performs inexpensive center/pose and bounding-box checks, then validates the
full swept offset set when required. The search is repeated with accumulated
soft penalties until `max_paths` feasible, sufficiently distinct paths are
collected, a retry limit is reached, or the computation budget is exhausted.

The facade shifts the raw endpoints outward when necessary to find collision-
free terminal poses, invokes the private search, then restores the continuous
26-connected terminal segments so callers receive paths beginning at the
requested start and ending at the requested goal. Conditional waypoints are
converted to the public `PoseDescription` representation and returned in
cost-sorted `PathResult` objects.

## 2. Core Data Structures

### Search state and action records

The fundamental lattice state is:

```text
state = (voxelIndex, poseId)
stateKey = voxelIndex * poseCount + poseId
```

`voxelIndex` identifies the anchor voxel in the dense grid. `poseId` indexes a
discrete tangent/normal orientation. A state record stores the best `g` cost,
parent state key, parent action, and open/closed status. The priority queue
orders nodes by `f`, then `h`, then `g`, with a deterministic state-key tie
break. A second dominance table groups pose orientations with identical
tangent and equivalent normal direction, avoiding redundant expansions at the
same voxel.

The action catalog contains:

- one straight action per pose, with lattice length and precomputed sweep;
- bend actions for each projected global axis, sign, and angle from one
  `angle_step_deg` through 90 degrees;
- two twist actions per pose (positive and negative roll step), where the
  tangent is preserved and only the cross-section normal changes;
- offset bounds for every sweep, used by the occupancy-block fast path.

### `ProcessedMap`

`ProcessedMap` is the public opaque handle. It owns a `std::shared_ptr` to a
private implementation (`ProcessedMap::Impl`) containing the loaded
`VoxelGrid`, generated `BusbarPose` table, and validated `PlannerConfig`.
Copies share immutable map context; route-specific penalty state is created
inside each `findPaths()` request and is never stored in the public handle.
This separation makes loading reusable and prevents one route request from
contaminating another.

### `VoxelClass`

The public voxel classification is:

```cpp
enum class VoxelClass {
    BLOCKED = 0,
    UNCONDITIONAL = 1,
    POSE_CONDITIONAL = 2
};
```

`BLOCKED` denotes raw obstacles, map-boundary failures, or voxels for which no
generated pose has a collision-free footprint. `UNCONDITIONAL` denotes an
open voxel whose complete pose table is valid. `POSE_CONDITIONAL` denotes an
open voxel with a non-empty subset of valid poses. Search still checks the
specific `poseId` even for an unconditional voxel because transitions also
carry a swept volume.

### `EndpointConstraint`

`EndpointConstraint` is an internal terminal contract used by the facade. Its
default start and end tangents are unit +X vectors, and its minimum beginning
and ending lengths are zero. The implementation validates unit tangents and
finite non-negative lengths, resolves safe terminal pose states, and applies a
15-voxel endpoint-immunity radius while validating terminal states. The raw
public endpoints remain the source of truth; safe shifted endpoints are an
internal planning detail restored before API return.

### Public result structures

`PathResult` contains a voxel centerline, cost, and pose descriptions. A
`PoseDescription` stores normalized normal and tangent vectors for conditional
waypoints. `PlanStatus::NO_PATH` is represented by an empty public path
vector; invalid arguments and resource failures are surfaced as exceptions by
the facade.

## 3. Swept-Volume Strategy (扫掠体碰撞检测)

Collision validation has two distinct geometric layers: a static pose
footprint at an anchor voxel and a continuous swept volume along an action.

### Static pose footprint: pure cross-section

For a pose with unit tangent `t` and unit normal `n`, the binormal is
`b = normalize(t x n)`. The rectangular busbar cross-section is rasterized in
the `(b, n)` plane. Width and thickness are sampled at one sample per voxel
dimension (with at least one sample), centered around zero, rounded to the
nearest voxel, and deduplicated in a set of integer offsets. This produces the
pose's pure cross-section only; it does not include any movement along `t`.

The resulting offset list is used to populate the lazy pose mask. A pose is
allowed at an anchor only when every offset is inside the map and free of a raw
obstacle. The mask is a feasibility cache, not a substitute for transition
validation.

### Straight movement and the 0.49-voxel sampling rule

Straight actions translate the cross-section along the pose tangent. The
continuous segment is sampled at:

```text
sampleCount = max(1, ceil(distance / 0.49))
sampleDistance(i) = distance * i / sampleCount,  i = 0..sampleCount
```

At every sample, the same oriented cross-section is rasterized at the sampled
center. The union of all rounded offsets is the straight swept volume. The
0.49-voxel maximum surface step is intentionally below one voxel: a moving
surface cannot advance across a voxel-sized gap without producing an
intermediate sample. The union is deduplicated before search-time use.

At search time, the sweep's integer bounding box is tested against the
8-voxel occupancy blocks first. Only if that box is not immediately clear are
the exact offsets enumerated. This preserves exact collision semantics while
making large open regions inexpensive.

### Bend and twist swept volumes

For a bend, the center follows a circular arc and the cross-section is
parallel-transported through the arc. The sample count is based on the arc
surface length, including the cross-section radius:

```text
samples = max(1, ceil(angle * (bendRadius + crossSectionRadius) / 0.49))
```

For a twist, the center remains fixed and the cross-section rotates around its
tangent. Its sampling count is proportional to angular travel times the
cross-section radius divided by 0.49. Static endpoint masks therefore answer
“can this pose sit here?”, while continuous sweeps answer “can the entire
physical motion between these states occur without collision?” Both checks are
required for a valid non-terminal transition.

## 4. Busbar Kinematics (巴片的旋转与弯曲)

### Discrete pose generation and roll

The tangent lattice contains all 26 non-zero integer directions
`(tx, ty, tz) in {-1, 0, 1}^3`. For each tangent, a deterministic orthogonal
reference axis creates a base normal and binormal. The normal is then rolled in
steps of `angle_step_deg` (default 15 degrees) through 360 degrees. The table
therefore contains `26 * (360 / angle_step_deg)` poses (624 at the default
step), each with integer tangent, unit normal, pose ID, and roll angle.

### Rodrigues axis-angle rotation

All bend and twist frame updates use Rodrigues rotation. For a unit axis `k`,
vector `v`, and angle `theta`:

```text
R(k, theta)v = v cos(theta)
             + (k x v) sin(theta)
             + k (k . v) (1 - cos(theta))
```

The implementation rotates tangents, normals, and radial vectors with this
formula, then normalizes frames and verifies orthogonality and endpoint
alignment within a small geometric tolerance.

### Flat bends and vertical bends

For each source pose, each global X/Y/Z axis is projected onto the plane
orthogonal to the current tangent. Both signs are considered. For every
discrete bend angle from one roll step through 90 degrees, the rotated tangent
and normal are matched to the closest generated target pose.

The bend radius is selected from the axis' alignment with the source frame:

```text
flatAlignment     = abs(axis . sourceNormal)
verticalAlignment = abs(axis . sourceBinormal)

radius = flatAlignment >= verticalAlignment
       ? flat_bend_factor     * busbar_thickness
       : vertical_bend_factor * busbar_width
```

The default factors are 1.5 for flat bends and 0.557143 for vertical bends.
The physical action length is `radius * angle`. Its lattice displacement is
the rounded circular-arc endpoint. Zero-displacement actions are discarded.
The bend sweep rotates the radial vector and transported frame at every 0.49
voxel surface step and validates the end tangent and end normal against the
target pose.

### Twists (roll actions)

Twists connect poses with the same integer tangent and adjacent roll values.
The twist angle is `+angle_step_deg` or `-angle_step_deg`; the center does not
translate. The end normal must equal the start normal rotated about the common
tangent. The twist sweep is the union of cross-sections during that roll, so a
wide busbar cannot rotate through an obstacle merely because both endpoint
footprints are individually clear.

### Transition validation

For every candidate action, the planner validates the target center and target
pose mask, then the action's exact sweep offsets. Terminal states may bypass
the normal obstacle rejection only inside the endpoint-immunity region; this
allows the facade to resolve safe outward terminal states in a narrow entry
channel without making the rest of the route permissive.

## 5. Macro-Diversity A* Search (多路径生成策略)

### Escaping the single-optimal-path trap

The planner does not run an independent unmodified A* search for every path.
It runs one search, records the accepted centerline, and adds its spatial cost
to a request-local `CostPenaltyMap`. Subsequent searches optimize physical
length plus those penalties. This preserves the best path first while making
reuse progressively unattractive, allowing the search to leave the first
optimal corridor and discover alternate corridors.

`max_paths` is bounded to ten per request. A candidate is compared with all
accepted paths by unique-voxel centerline overlap:

```text
overlapRatio = sharedUniqueCandidateVoxels /
               uniqueCandidateVoxels
```

Candidates above `PlannerConfig::max_overlap_ratio` (default 0.30) are
rejected, penalized, and retried. Up to ten diversity retries are attempted
for the current output slot; this prevents an impossible diversity demand from
running indefinitely. Accepted paths are finally sorted by reported cost.

### Soft-penalty strategy

For every voxel on an accepted (or rejected over-overlap) centerline, the
planner visits its 3 x 3 x 3 neighborhood. Outside endpoint immunity:

```text
centerline voxel reuse:  +100000.0
neighbor voxel reuse:    +1000.0
```

The penalty is added to the target-voxel component of a transition's `g`
cost. The centerline penalty is intentionally massive, while the neighbor
penalty creates a soft exclusion band rather than a hard wall. This lets a
later route share necessary infrastructure when no alternative exists, but
strongly favors a genuine geometric departure where one is available.

### Chokepoint sharing and endpoint immunity

The facade passes `endpointImmunityRadius = 15.0` voxels to both endpoint
validation and diversity penalties. Penalty updates skip every candidate voxel
whose Euclidean distance from the raw start or raw goal is at most that radius.
Consequently, multiple paths may share a narrow terminal channel, exit the
start region, diverge in open space, and converge again near the goal. The
immunity applies only to reuse penalties and terminal acceptance; ordinary
non-terminal transitions still undergo pose and swept-volume validation.

The backward field also treats obstacles inside this endpoint region as
traversable for its guidance computation. This keeps the heuristic from
declaring a legitimate terminal approach unreachable merely because the raw
endpoint is embedded in a tight or partially occupied terminal structure.

### Backward 3D BFS heuristic

Before the first A* iteration, a reverse breadth-first search starts at the
goal voxel and visits the six axis-aligned neighbors. It stores a compact
16-bit distance array when the maximum grid Manhattan diameter fits; otherwise
it uses 32-bit distances. Obstacles are excluded except for the start/goal and
the endpoint-immunity region. The field records both distance-to-goal and the
number of reachable voxels.

The SE(3) heuristic combines that field with orientation and travel-heading
terms:

```text
guidedDistance = backwardBfsDistance, or Euclidean distance if unreachable
orientationPenalty = (1 - dot(currentTangent, goalTangent))
                     * heuristicScale * 5
travelHeadingPenalty = (1 - bestHeadingAlignment)
                       * heuristicScale * 20
h = heuristicScale * guidedDistance
    + orientationPenalty
    + travelHeadingPenalty
```

For the heading term, the planner looks for a neighboring voxel whose BFS
distance is lower than the current distance and rewards alignment with that
descending direction. If no descending neighbor is available, it falls back
to the direct current-to-goal heading. Thus a dead-end state is not guided by
straight-line distance alone: the heuristic points expansions toward the
known traversable basin and supplies a useful fallback when the voxel is
outside the BFS-reachable component.

### Search limits and robustness

The unpenalized search has a 50,000-expansion budget; penalized searches have a
100,000-expansion budget. The planner records the closest generated goal
distance, backward distance, generated-state count, and open-list size when a
limit is reached, and returns paths already completed when possible. Neighbor
candidate preparation is parallelized for sufficiently large candidate sets,
while collision checks remain deterministic through cached action geometry and
stable queue tie-breaking.

## 6. Performance and Correctness Invariants

- Every non-terminal transition must pass both the target pose mask and the
  continuous swept-volume collision test.
- Every generated frame keeps tangent and normal orthogonal and normalized;
  bend endpoints must agree with the selected target pose within geometric
  tolerance.
- Integer sweep offsets are deduplicated and accompanied by min/max bounds for
  constant-time broad-phase rejection.
- Pose masks and action sweeps are reused across all macro-diversity retries;
  only the penalty map and A* search records change between iterations.
- Public paths are restored to the caller's exact raw endpoints and public
  conditional-pose indices are remapped after terminal extension.

These invariants preserve the distinction between geometric feasibility,
kinematic continuity, and route diversity while keeping the public API small
and stable.
