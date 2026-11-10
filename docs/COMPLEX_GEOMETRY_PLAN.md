# Complex-geometry support: branch plan

Status: B0–B1b, B3, B4, the B5 storage memo, the B7a checkpoint identity and the patch binding are landed (2026-11-07); B6a–c and B7b remain deferred | Updated: 2026-11-07 | Target: split complex-geometry support into independently mergeable feature branches

> Exit condition: a complex-geometry case runs from a repository geometry file with measured fill ratio,
> entity-bound boundary conditions, link-accurate walls, and an agreed memory budget; every non-spike
> branch below has landed with its gate. Delete this plan when the last branch merges.

## 1. Decision

Complex geometry is four separable concerns, not one missing mesh reader:

1. **Geometry / entities** — a domain-sized geometry record owns role, entity, material, and a narrow-band
   surface description. The domain is supplied by the caller, not by one STL bounding box.
2. **Active ownership** — a block owns an explicit active-cell set; decomposition and load balance use
   owned active cells, not AABB cells.
3. **Boundary physics** — a surface-link table carries distance/fraction, wall velocity, and entity;
   interpolated bounce-back replaces staircase reflection.
4. **Storage / execution** — dense AABB remains the default; active-cell storage is a policy. Brick-sparse
   is a measured choice, not an assumed win: the baseline is the current block decomposition, and B5
   kept dense.

The test geometry is our own closed square helical duct, not a private asset. It is tracked under
`tests/feature/geo_stl/stl/` and used by `tests/feature/geo_stl` for polarity and fill-ratio checks.

## 2. Results

The fixture is measured by `tests/feature/geo_stl`; owned fill is owned active / allocated.

| Resolution | Blocks | Allocated | Owned | Owned fill | Active fill |
| --- | --- | --- | --- | --- | --- |
| h = 1.0 | 11 | 12540 | 1294 | 0.103 | 0.143 |
| h = 0.5 | 22 | 51510 | 10588 | 0.206 | 0.262 |

- **B0**: fixture polarity and fill measured; active cells scale 7.5x from h=1.0 to h=0.5. The
  compaction path remains opt-in, not production storage.
- **B2a**: owned cells are counted from each block's base AABB; the derived pointer table is
  invalidated on every block replacement. `voxel_map` and `brick_seam` cover the counts and the
  seam/RCB paths.
- **B2b**: each plan carries its own generation stamp; a shared pair survives only when both
  endpoints are active (each side tests its own field). The sparse pattern prunes 426 entries
  (102 shared pairs + 324 MPI cells); a receiver-side test would have mis-pruned 26 live pairs.
  `compact_mpi` (12 checks) and `plan_stamp` (5 checks) cover the compact and stale sequences.
- **B1a**: `VoxelGeometry` plus its builder; `writeVoxelFlags` fills the field manager without a
  geometry-to-field include; `InitVoxelMap` derives the active map. On 27 x 27 x 36 the record has
  Domain 24950, Solid 1294, band 3947, and 2724 patch cells. `voxel_entity` (20299 checks) and
  `voxel_mask` (5 checks, 0 mismatches) cover the record and the mask equivalence.
- **B1b**: FLBGEO v1 is the geometry file; the round-trip is covered by `io/geometry_file`
  (105455 checks). casegen renders a `.geo` recipe from `[geometry.source.*]` and
  `tools/geometry_build` rebuilds the `.flbgeo`; the `casegen_geo_box3d` golden runs 100 steps.
- **B3**: `FillAwareDecomposition` cuts a block only when its whole split subtree lowers the
  dense footprint (base box plus overlap halo); the table is geometry-derived and independent of
  the rank count. `geo_bulk` goes from owned fill 0.406 (8 blocks) to 0.453 (16) and 0.511 (24);
  `geo_thick` stays at 0.658 (8 blocks) because no cut lowers its footprint. The helical duct
  reaches only 0.210 (h=0.5, 32 blocks) and 0.126 (h=1.0, 24 blocks), below the hard floors
  0.25/0.14. RCB on the same flag field is worse (0.074-0.081), so the thin limit is the dense
  AABB shape, not the cut order; B5 compared the storage layouts at this point. `block/decompose_coverage` (58 checks, serial and two
  ranks) requires one owner per active cell and a rank-identical layout.
  `BlockLayoutEntry`/`GetBlockLayout` record the table.
- **B3 fixtures**: `geo_bulk` (32^3 domain, 24^3 solid) has raw fluid 0.578 and a one-block
  owned fill 0.482; `geo_thick` has two 4-cell plates with entities 2/3; `geo_thin` has two
  1-cell plates with entities 2/3. All three go through `writeVoxelFlags` and `InitVoxelMap`.
- **B5**: the storage spike is recorded in `docs/STORAGE_COMPARE_PLAN.md` and reproducible from
  `benchmarks/geo_storage`. On the fill-aware table the existing compact iteration is 0.98-1.06x
  dense MLUPS; the active-slot model is 0.73-0.875x dense field bytes and excludes metadata, while
  brick is 0.51-0.77x but pays 1.4-3.3x seam traffic. Dense remains the default; B6a is deferred,
  with active-slot the first candidate if a case becomes memory-bound. The thin-duct floor stays
  red at 0.210; that is a block-shape limit, not a field-layout one.
- **B7a**: checkpoint format v2 stores the geometry identity next to the fields: grid, origin and
  extent, voxel size, the per-cell entity array, every boundary patch and the record digest.
  `writeWithGeometry` writes it; `readWithGeometry` and `checkCompatibleWithGeometry` validate it,
  and a digest mismatch is a hard error. A snapshot that carries an identity cannot be resumed
  through the record-less `read`/`checkCompatible`, so the old "trust the replay" path cannot
  silently load a different geometry. The block mapping now compares the global block count and
  each local block's presence instead of requiring one rank to own the whole table, so any rank
  assignment of the same block table resumes (2->1, 1->2, 2->2). `feature/checkpoint_geometry`
  covers the header, the mismatch and a bit-identical restart; `block/checkpoint_geometry_mpi`
  covers two ranks. Format v1 files are rejected.
- **Patch binding**: `flbbc::bindPatches` maps `PatchKind` to a task flag and replaces the domain bit
  (like the AABB-painted cases) on the boundary-layer cells the record maps to that patch; band
  cells that do not touch a non-Domain cell are skipped, and a kind without a mapping is counted,
  not silently dropped.  A per-entity override splits two patches of the same kind, and `skipFlag`
  keeps e.g. a periodic face intact.  `[patch.*]` sections (entity/kind/box) let
  `tools/geometry_build` carry patches in `.flbgeo`.  `feature/patch_bc` builds that file with two
  Wall patches (246 and 116 band cells), checks that the override gives each its own flag and that
  exactly the boundary cells are painted, and that an unmapped kind paints nothing.
- **Example migration**: `examples/cavity3d` reads `cavity3d.flbgeo` (box domain plus an internal
  solid obstacle) instead of building the box in code; `[Mesh]` is validated against the file, the
  obstacle's Wall patch drives 12364 bounce-back cells, the moving lid stays AABB-driven because a
  domain-edge wall has no surface band, and checkpoint/restart uses the geometry identity.  A
  100-step run on the 100^3 grid and a restart from step 100 were exercised.
- **Moving-wall example**: `examples/couette3d` is a channel periodic in x/y with two solid slabs;
  both patches are Wall and the per-entity override gives the top slab the moving-wall task and the
  bottom the static one.  The centreline profile after 2000 steps follows the Couette line (lattice
  wall speed 0.002, values rise 0 -> 0.002 across the gap).  The x/y periodic face strip carries no
  wall task because the periodic copy refreshes it, so the printed profile is the centreline, not a
  plane average.
- **Gates**: `make test` PASS 52 / FAIL 0 / SKIP 0 (CPU, MPI and CUDA); `make test-casegen` exits 0
  with 14 golden cases (PASS 14 / FAIL 0 / SKIP 0) and the generator unit tests;
  `tools/check_headers.sh` exits 0 (the one baseline `utils/cuda_device.h` entry is unchanged).

Patches and entity/material values round-trip but do not drive tasks or boundary conditions yet.
Existing `examples/` still build masks through the helper flag path; example migration is deferred.

## 3. Target architecture

**Geometry.** `VoxelGeometry` is immutable and domain-sized: role (domain/solid/outside), entity,
material, and boundary patches. Patches overlay cells and never change occupancy. A narrow band
stores signed distance and normal near the surface. The builder lives in the IO layer because it
consumes `StlReader`; the record itself has no IO dependency. It accepts several STLs and
primitives and validates closure. Minimum wall thickness is not checked yet. Patches drive the
task/BC flags through `flbbc::bindPatches` (kind -> flag bit, unmapped kinds are counted), entities
drive the link wall velocity, and `[patch.*]` sections make the `.geo` recipe carry both.

**File format.** FLBGEO v1 is the decided primary format: one versioned, self-contained binary
carrying grid/domain, role/entity/material, patches, the narrow band, the build recipe, and a
digest. A case can be rebuilt and audited from one file. `.vti` is an optional visualization
export, not an interchange format: the repository has no VTK ImageData reader, and patches, the
recipe, and the digest would need a sidecar.

**Ownership.** `ownedN` (code `OwnedNum`) is the unique owned active count per block. `activeN`
(code `VoxNum`) is owned plus the halo slots whose owner cell is active, not every halo slot in
the block AABB. Every active interior cell belongs to exactly one block; each block's stored set
is closed under the stencil. The mask is built from geometry occupancy, not runtime fluid/gas
state. A geometry generation bumps only on a geometry/mask rebuild; runtime free-surface
reclassification does not bump it. Free surface requires overlap at least two.

**Boundary.** A surface link records direction, distance fraction, wall velocity, entity, and
whether the direction is blocked. Links are block-local and resolve a halo slot across a seam.
Interpolated bounce-back is the default for curved walls; sub-cell walls are either rejected at
build time or must be honoured by streaming and communication.

**Storage.** Dense AABB is the default and bit-identical for a fixed decomposition. Active-cell
storage uses a slot map. Brick-sparse needs an upstream halo or boundary-plane exchange; the
measured seam fraction, not only `(B+2)^3/B^3`, decides it.

**Communication and IO.** Plans are derived data stamped per plan with the geometry generation; a
mismatch aborts and the caller rebuilds. Sender and receiver apply the same pair mask; a size
mismatch is a hard error. A block layout record makes decomposition reproducible. Checkpoint
format v2 carries the geometry identity (grid, origin/extent, per-cell entity, patches and a record
digest) next to the fields and validates it on restart; the active mask and policy tag are B7b.
Only a same-block-table restart is promised, but any rank assignment of that table is allowed.

## 4. Branches

| Branch | Scope | Depends on | Exit | Gate |
| --- | --- | --- | --- | --- |
| `feat/geo-observe` (B0) | fill ratio, active/allocated counts, effective MLUPS, peak RSS, synthetic STL fixture, safe opt-in compaction | — | fixture polarity/fill measured; fixture and compact-path tests pass | `make test-cpu`; `tools/check_headers.sh`; landed |
| `feat/voxel-entity` (B1a) | in-memory geometry, builder, entities, materials, boundary patches, surface band | B0 | two entities/materials coexist; patches select cells; 2D keeps the existing path | `make test-cpu`; `feature/voxel_entity`, `feature/voxel_mask`; landed |
| `feat/voxel-io` (B1b) | FLBGEO read/write, casegen spec, `tools/geometry_build`, golden case; example migration landed in `examples/cavity3d` | B1a | file round-trip; one case driven from the file | `make test-casegen`; `make test-cpu`; landed |
| `feat/active-index` (B2a) | owned/active counts, inverse map, owned metrics, derived-state invalidation | B0 | owned/active metrics and pointer-table invalidation covered by tests | `make test-cpu`; `feature/brick_seam`; landed |
| `feat/active-comm` (B2b) | joint shared pair mask, per-plan generation stamp, `InitVoxelMap` invalidation, 2-rank compact smoke test | B2a | plan stamp aborts on mismatch; shared pairs require both endpoints active; compact MPI test passes | block MPI suite; `tests/block/compact_mpi`; `tests/block/plan_stamp`; landed |
| `feat/geo-decompose` (B3) | fill-aware split (subtree lookahead), coverage assertion, layout record, owned-count budget | B0, B2a, B2b | bulk/thick footprint within budget; block cap respected; coverage holds; thin hard floor **not met**, handed to B5 | block tests; `block/decompose_coverage`; landed |
| `feat/link-bc` (B4) | surface links, interpolated bounce-back, thin/corner, moving wall (constant and rigid rotation) | B1a | curved-wall accuracy; thin-wall mass conservation; seam fixture on one and two ranks; CUDA scope stated | `make test`; `feature/link_bc*`, `block/link_bc_seam`; landed |
| `spike/storage-compare` (B5) | time-boxed comparison of the current block decomposition, dense, active-slot, and brick layouts; seam copy vs row-base redirect; thin and thick fixtures | B1a, B2b, B3 | memo with resident bytes, MLUPS, seam bytes/step, communication, and complexity; one policy chosen, including "dense remains the default, B6a deferred" | `docs/STORAGE_COMPARE_PLAN.md`; `benchmarks/geo_storage`; memo landed |
| `feat/active-storage-cpu` (B6a) | policy chosen by B5: fields by slots/bricks, stream, collision, cell access; no free surface | B5, B3 | memory target met; no mixed dense/sparse; CPU checksums | `make test-cpu`; CUDA compile where available |
| `feat/active-storage-physics` (B6b) | free surface, CA, statistics, writers, checkpoint consumers | B6a | free-surface and CA fixtures pass; memory report | `make test-cpu`; checkpoint tests |
| `feat/active-cuda` (B6c) | device mirror and kernels for the chosen storage | B6b | full test suite and FP16 run green; CPU/GPU checksum | CUDA gate |
| `feat/geo-io` (B7a/B7b) | B7a: checkpoint geometry identity (grid/origin, per-cell entity, patches, digest), format v2, any rank assignment of the same block table; B7b: active mask and policy tag | B7a: B1b, B3; B7b: B6a | B7a: digest mismatch rejected, 2->1 / 1->2 / 2->2 restarts covered; B7b: mask round-trips | B7a: `feature/checkpoint_geometry`, `block/checkpoint_geometry_mpi`, checkpoint suite; landed |

B4 was merged from `feat/link-bc` and B3 from `feat/geo-decompose`. B4 needs B1a only: the link
table reads the halo flags that `writeVoxelFlags` fills from the record, not the B2b seam plan. B5
is a spike: its memo and benchmark are on `dev`, no storage code is merged. B7 does not promise
cross-decomposition restart; that is a separate cell-coordinate remap. Every MPI gate also needs
one explicit `mpirun` run in the runner's form; CUDA is not a gate before B6c.

B4 landed with three recorded gaps: the moving-wall term is 9.5% L2 in the flat Couette fixture and
9.7% on the curved inner half (the `q = 0.75` one-sided link, a uniform slope deficit of 0.921);
925 of 7160 duct links fall back to the half-way position (a further 671 with `q > 1` are clamped to
the link end, down from 1596); and the apply path is host only. Patch binding closed the fourth gap: a
patch kind now paints the task flag, and `[patch.*]` carries the patch through the `.geo` recipe. Curved-wall accuracy is 7.5% L2 against Taylor-Couette, and the two-rank seam
reproduces the single-block link count and mass.

B1a was merged from `feat/voxel-entity`. The StlReader repair, the mask adapter, B1b,
`geometry_build`, and the casegen changes were committed directly on `dev`; the `feat/voxel-io` row
names the planned scope, not an existing branch. B3's row keeps its thin exit intentionally red:
that was the first input to B5.

## 5. Sequencing

1. **B0, B2a, B2b, B1a, B1b, B3, B4** have landed: fixture, owned index, per-plan stamping, shared
   pair mask, voxel record/builder, the FLBGEO path through casegen, the fill-aware dense
   decomposition, and the surface-link boundary. B3 met bulk and thick; the thin-duct hard floor is
   not met.
2. **B5 is complete**: `docs/STORAGE_COMPARE_PLAN.md` and `benchmarks/geo_storage`. Dense remains
   the default and B6a is deferred. **B7a has landed**; patch binding is on `feat/patch-bc`.
3. **B6a, B6b, B6c** are deferred until a case is measured memory-bound or the thin fixture must
   meet its 0.25 floor; then active-slot is the first candidate. **B7b** follows B6a.
4. B4's open refinements (moving-wall `q = 0.75`, patch-to-BC binding, duct fallback) are follow-ups
   on `dev`, not new branches.

## 6. Rejected alternatives

- **Unstructured mesh.** Rejected: LBM needs a regular lattice; constant offsets and rotation are the
  performance model. A different solver is out of scope.
- **The current compaction map as the end state.** Rejected: it compacts iteration and part of
  communication, not storage, streaming, free surface, CUDA, or IO. It remains an opt-in stepping stone.
- **Solid-cell elision for a fluid-domain STL.** Rejected for this fixture: the STL encloses the fluid
  domain. The idea remains valid for embedded solids and is out of scope.
- **Widening the flag byte only.** Rejected: more bits do not add entity semantics, surface distance, or
  persistence.
- **A full sparse neighbour table as the first storage step.** Rejected as a first step: its index memory
  rivals the populations and worsens locality. B5 measured brick seam traffic and deferred it.
- **Brick-sparse without an upstream halo.** Rejected: rotation would move values across a brick face
  without the upstream data.
- **Cell-level compact active set with one global cyclic shift.** Rejected: the stencil translations
  generate the full lattice; a full rectangular brick is the finest unit that keeps the rotation.
- **Per-cell pointer/index redirect as the first storage scheme.** Rejected: slower than a structured
  seam copy; the row-base redirect is at best ~1–2% faster and needs a read-before-write phase.
- **Unbounded block splitting.** Rejected: overlap cost and block count grow; B3 stops at a measured
  fill/overhead optimum.
- **`.vti` as the primary geometry format, or HDF5 / JSON plus raw arrays.** Rejected: the repository
  has no VTK ImageData reader, and image data would still need a sidecar for patches, recipe, and the
  versioned digest; HDF5 adds a dependency, and split manifest/array files lose the self-contained
  record and digest. `.vti` stays an optional export.

## 7. B3 budget and open questions

Starting budget, measured on the helical duct above. Owned fill is owned active / allocated;
allocated cells are the primary footprint metric because peak RSS is dominated by fixed overhead
at this size.

| Budget | Current (h=0.5 / h=1.0 unless stated) | B3 target | Hard limit | B3 measured |
| --- | --- | --- | --- | --- |
| Owned fill, thin duct | 0.206 / 0.103 | 0.28–0.32 / 0.16–0.20 | ≥ 0.25 / ≥ 0.14 | 0.210 (32 blocks) / 0.126 (24) — **not met** |
| Owned fill, bulk fixture | 0.482 (one block); 0.365–0.406 (8–16 blocks) | 0.45–0.60 | ≥ 0.40 | 0.453 (16) / 0.511 (24) — met |
| Allocated cells | 51510 / 12540 | 0.64–0.73× / 0.52–0.65× | ≤ 0.82× / ≤ 0.74× | 50459 (0.98×) / 10261 (0.82×) — thin not met |
| Peak RSS | 9.66 MB (`geo_stl`, process-wide) | ≤ 0.90× | ≤ 1.00× | B5 recorded field bytes 6.8-8.6 MB; RSS was not measured for a storage alternative, so the row stays open |
| Blocks, one-rank fixture | 22 / 11 | 16–24 / 8–16 per rank | ≤ 32; ≥ 4 for MPI | 32 / 24 — cap respected, thin target not met; bulk 16 / 24 |

The allocated-cell ranges are the owned-fill targets expressed against the current allocation; the
hard limit is the fill floor. Uniform baseline (one rank, 32^3, h=1, overlap 1): `geo_bulk` at
8/16/27/64 blocks has owned fill 0.406/0.365/0.345/0.296; `geo_thick` has
0.658/0.593/0.560/0.480. B3's split reaches 0.453 at 16 blocks and 0.511 at 24 against the
8-block 0.406, so the bulk target is met; `geo_thick` stays at 0.658 because no cut lowers its
footprint. The thin duct stays at 0.210/0.126, and RCB on the same flag field is worse
(0.074-0.081); the dense AABB shape, not the cut order, is the limit. B3 therefore does not relax
the thin gate. B5 compared the storage layouts and kept dense (`docs/STORAGE_COMPARE_PLAN.md`), so
the thin limit is still the block shape, not the field layout. Peak RSS is the `geo_stl` process
value (9.66 MB); it is not a B3 storage result and B5 did not measure an alternative.

- CUDA scope and whether GPU complex geometry is required with the CPU feature.
- Thin-wall fixture: `tests/feature/geo_thin` (two 1-cell plates) is now tracked; bulk/thick and
  multi-entity remain in `tests/feature/geo_bulk` and `tests/feature/geo_thick`.
- B3 implements a split pass with a four-level subtree lookahead. A merge pass was not needed for
  bulk/thick; B5 did not choose a storage change, so merge stays a B3 follow-up only if the thin
  gate is reopened.
- B3's `MaxBlocks` and `TargetOwnedFill` are helper parameters; the `.ini`/casegen keys are not
  wired yet.
- 2D parity: the existing 2D path stays unchanged until a later branch.
