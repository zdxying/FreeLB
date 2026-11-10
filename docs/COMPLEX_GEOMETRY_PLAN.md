# Complex-geometry support: branch plan

Status: B0 + B2a + B2b + B1a + B1b landed on `dev` (2026-11-06); B3 next | Updated: 2026-11-06 | Target: split complex-geometry support into independently mergeable feature branches

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
   is a measured choice, not an assumed win: the baseline is the current block decomposition, and B5 may
   keep dense.

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
- **Gates**: `make test-cpu` PASS 36 / FAIL 0 / SKIP 4; `make test-casegen` 185 unit tests and
  14 golden cases (PASS 14 / FAIL 0 / SKIP 0); `tools/check_headers.sh` exits 0.

Patches and entity/material values round-trip but do not drive tasks or boundary conditions yet.
Existing `examples/` still build masks through the helper flag path; example migration is deferred.

## 3. Target architecture

**Geometry.** `VoxelGeometry` is immutable and domain-sized: role (domain/solid/outside), entity,
material, and boundary patches. Patches overlay cells and never change occupancy. A narrow band
stores signed distance and normal near the surface. The builder lives in the IO layer because it
consumes `StlReader`; the record itself has no IO dependency. It accepts several STLs and
primitives and validates closure. Minimum wall thickness is not checked yet; patches and entities
are not bound to tasks or boundary conditions until B4.

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
mismatch is a hard error. A block layout record makes decomposition reproducible. Checkpoint adds
the geometry digest and active mask; only a same-block-table restart is promised.

## 4. Branches

| Branch | Scope | Depends on | Exit | Gate |
| --- | --- | --- | --- | --- |
| `feat/geo-observe` (B0) | fill ratio, active/allocated counts, effective MLUPS, peak RSS, synthetic STL fixture, safe opt-in compaction | — | fixture polarity/fill measured; fixture and compact-path tests pass | `make test-cpu`; `tools/check_headers.sh`; landed |
| `feat/voxel-entity` (B1a) | in-memory geometry, builder, entities, materials, boundary patches, surface band | B0 | two entities/materials coexist; patches select cells; 2D keeps the existing path | `make test-cpu`; `feature/voxel_entity`, `feature/voxel_mask`; landed |
| `feat/voxel-io` (B1b) | FLBGEO read/write, casegen spec, `tools/geometry_build`, golden case; example migration deferred | B1a | file round-trip; one case driven from the file | `make test-casegen`; `make test-cpu`; landed |
| `feat/active-index` (B2a) | owned/active counts, inverse map, owned metrics, derived-state invalidation | B0 | owned/active metrics and pointer-table invalidation covered by tests | `make test-cpu`; `feature/brick_seam`; landed |
| `feat/active-comm` (B2b) | joint shared pair mask, per-plan generation stamp, `InitVoxelMap` invalidation, 2-rank compact smoke test | B2a | plan stamp aborts on mismatch; shared pairs require both endpoints active; compact MPI test passes | block MPI suite; `tests/block/compact_mpi`; `tests/block/plan_stamp`; landed |
| `feat/geo-decompose` (B3) | fill-aware split/merge, owned-count load balance, coverage assertion, layout record | B0, B2a, B2b | dense footprint within budget; block cap respected; coverage holds | block tests, checksums, MPI scaling |
| `feat/link-bc` (B4) | surface links, interpolated bounce-back, thin/corner, moving wall | B1a; B2b for seam-crossing links | curved-wall accuracy; thin-wall mass conservation; seam fixture; CUDA scope stated | `make test`; new link-BC baseline |
| `spike/storage-compare` (B5) | time-boxed comparison of the current block decomposition, dense, active-slot, and brick layouts; seam copy vs row-base redirect; thin and thick fixtures | B1a, B2b, B3 | memo with resident bytes, MLUPS, seam bytes/step, communication, and complexity; one policy chosen, including "dense remains the default, B6a deferred" | review only |
| `feat/active-storage-cpu` (B6a) | policy chosen by B5: fields by slots/bricks, stream, collision, cell access; no free surface | B5, B3 | memory target met; no mixed dense/sparse; CPU checksums | `make test-cpu`; CUDA compile where available |
| `feat/active-storage-physics` (B6b) | free surface, CA, statistics, writers, checkpoint consumers | B6a | free-surface and CA fixtures pass; memory report | `make test-cpu`; checkpoint tests |
| `feat/active-cuda` (B6c) | device mirror and kernels for the chosen storage | B6b | full test suite and FP16 run green; CPU/GPU checksum | CUDA gate |
| `feat/geo-io` (B7a/B7b) | B7a: checkpoint origin/entity/patch/digest, same block-table restart; B7b: active mask and policy tag | B7a: B1b, B3; B7b: B6a | different rank count restarts; digest mismatch rejected; mask round-trips | checkpoint tests; format-version approval |

B5 is not merged. B7 does not promise cross-decomposition restart; that is a separate cell-coordinate
remap. Every MPI gate also needs one explicit `mpirun` run in the runner's form; CUDA is not a gate
before B6c.

B1a was merged from `feat/voxel-entity`. The StlReader repair, the mask adapter, B1b,
`geometry_build`, and the casegen changes were committed directly on `dev`; the `feat/voxel-io` row
names the planned scope, not an existing branch.

## 5. Sequencing

1. **B0, B2a, B2b, B1a, B1b** have landed: fixture, owned index, per-plan stamping, shared pair mask,
   voxel record/builder, and the FLBGEO path through casegen.
2. **B3** is next. **B4** follows B1a and uses B2b for seam-crossing links. **B7a** can follow B1b/B3.
3. **B5** after B3, with the memory budget and brick parameters fixed beforehand.
4. **B6a, B6b, B6c** after the chosen policy; **B7b** follows B6a.

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
  rivals the populations and worsens locality. Compare it in B5.
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

| Budget | Current (h=0.5 / h=1.0) | B3 target | Hard limit |
| --- | --- | --- | --- |
| Owned fill, thin duct | 0.206 / 0.103 | 0.28–0.32 / 0.16–0.20 | ≥ 0.25 / ≥ 0.14 |
| Owned fill, bulk fixture | TBD | 0.45–0.60 | ≥ 0.40 |
| Allocated cells | 51510 / 12540 | 0.64–0.73× / 0.52–0.65× | ≤ 0.82× / ≤ 0.74× |
| Peak RSS | ~8.8 MB | ≤ 0.90× | ≤ 1.00× |
| Blocks, one-rank fixture | 22 / 11 | 16–24 / 8–16 per rank | ≤ 32; ≥ 4 for MPI |

The allocated-cell ranges are the owned-fill targets expressed against the current allocation; the
hard limit is the fill floor. The bulk target waits for the second fixture; optimizing only the thin
duct would overfit. If B3 cannot reach the hard floor, B5 compares active-slot/brick storage
instead of relaxing the gate.

- CUDA scope and whether GPU complex geometry is required with the CPU feature.
- A second fixture: thin walls/multi-entity contact, plus a thick/bulk case for B5.
- 2D parity: the existing 2D path stays unchanged until a later branch.
