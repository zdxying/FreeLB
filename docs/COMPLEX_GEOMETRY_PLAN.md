# Complex-geometry support: branch plan

Status: B0 + B2a + B2b landed (2026-11-05); B1a/B3 next | Updated: 2026-11-05 | Target: split complex-geometry support into independently mergeable feature branches

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

### B0

- Fixture: 1032 triangles, closed and watertight, bbox 24.04 x 24.04 x 33.74, analytic volume/bbox 6.84%.
- Cell size 1.0: active 1789 / allocated 12540, block-local fill 0.143.
- Cell size 0.5: active 13501 / allocated 51510, block-local fill 0.262; active cells scale 7.5x.
- Fill rises with resolution because allocation is block-local and halo-inclusive; active-cell scaling is
  the resolution check, not equal fill.
- `tests/feature/geo_stl`: 21 checks; `tests/block/voxel_map`: 433 checks at B0.
- The compaction path remains opt-in and is not production storage.

### B2a

- `VoxelMap` now has AABB -> slot, slot -> AABB, and owned/halo identity.
- Owned cells are counted from each block's base AABB, not from the overlap-inclusive block AABB.
- Cell size 1.0: active 1789, owned 1294 / allocated 12540, owned fill 0.103.
- Cell size 0.5: active 13501, owned 10588 / allocated 51510, owned fill 0.206; owned cells scale 8.2x.
- The derived pointer table is invalidated by `CreateBlocks`, `Optimize`, `LoadOptimization`, and merge
  paths; the 2D helper has the same guard.
- `tests/block/voxel_map`: 644 checks; `tests/feature/brick_seam`: 510 checks.

### B2b

- `InitVoxelMap` bumps a geometry generation; lattice and field communication plans record the generation
  they were built against and abort on a mismatch.
- Shared communication pairs are kept only when both the send and receive cells are active.
- `tests/block/compact_mpi` runs on 2 ranks with `InitVoxelMap` before the plan: 1296 allocated,
  1152 active, 72 communicator pairs removed, plan stamp accepted, compact count matches.
- Existing MPI tests compile and pass against the stamped plans.
- `make test-cpu`: PASS 31, FAIL 0, SKIP 4 (the CUDA targets).

The active index, plan stamping, and the compact MPI case are complete.

## 3. Target architecture

**Geometry.** `VoxelGeometry` is immutable and domain-sized: role (domain/solid/outside), entity, material,
and boundary patches. Boundary patches are overlays on entity surfaces; they never change occupancy. A
narrow band stores signed distance and normal near the surface. The builder accepts several STLs and
primitives, transforms them into one grid, and validates closure and minimum wall thickness. Readers live
in the IO layer; geometry does not include IO in new code.

**Ownership.** `ownedN` is the unique owned active count per block; `activeN` is owned plus halo slots.
Every active interior cell belongs to exactly one block; each block's stored set is closed under the
stencil. The mask is built from geometry occupancy, not from runtime fluid/gas state. All fields of a
block share one slot map. A geometry generation bumps only on a geometry/mask rebuild; runtime
free-surface reclassification does not bump it. Free surface requires overlap at least two.

**Boundary.** A surface link records direction, distance fraction, wall velocity, entity, and whether the
direction is blocked. Links are block-local and resolve a halo slot across a seam. Interpolated
bounce-back is the default for curved walls; sub-cell walls are either rejected at build time or must be
honoured by streaming and communication.

**Storage.** Dense AABB is the default and bit-identical for a fixed decomposition. Active-cell storage
uses a slot map. Brick-sparse needs an upstream halo or boundary-plane exchange; the measured seam
fraction, not only `(B+2)^3/B^3`, decides it. A full rectangular brick is the finest unit that keeps the
per-brick rotation.

**Communication and IO.** Plans are derived data stamped with the geometry generation; a mismatch aborts
and the caller rebuilds. Sender and receiver apply the same pair mask; a size mismatch is a hard error.
A geometry file stores grid, role/entity/material, patches, build recipe, and the narrow band. A block
layout record makes decomposition reproducible. Checkpoint adds geometry digest and active mask; only a
same-block-table restart is promised.

## 4. Branches

| Branch | Scope | Depends on | Exit | Gate |
| --- | --- | --- | --- | --- |
| `feat/geo-observe` (B0) | fill ratio, active/allocated counts, effective MLUPS, peak RSS, synthetic STL fixture, safe opt-in compaction | — | fixture polarity/fill measured; fixture and compact-path tests pass | `make test-cpu`; `check_headers`; landed |
| `feat/voxel-entity` (B1a) | in-memory geometry, builder, entities, materials, boundary patches, surface band | B0 | two entities/materials coexist; patches select cells; 2D keeps the existing path | `make test-cpu`; new feature test |
| `feat/voxel-io` (B1b) | geometry file read/write, casegen spec, example migration | B1a | file round-trip; one case driven from the file | IO and casegen tests |
| `feat/active-index` (B2a) | owned/active counts, inverse map, owned metrics, derived-state invalidation | B0 | owned/active metrics and pointer-table invalidation covered by tests | `make test-cpu`; `feature/brick_seam`; landed |
| `feat/active-comm` (B2b) | joint shared pair mask, generation-stamped plans, `InitVoxelMap` invalidation, 2-rank compact smoke test | B2a | plan stamp aborts on mismatch; shared pairs require both endpoints active; compact MPI test passes | block MPI suite; `tests/block/compact_mpi`; landed |
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

## 5. Sequencing

1. **B0, B2a, B2b** have landed: fixture, owned index, plan stamping, and shared pair mask.
2. **B1a and B3** next. B1b follows B1a; B4 follows B1a and uses B2b for seams.
3. **B5** after B3, with the memory budget and brick parameters fixed beforehand.
4. **B6a, B6b, B6c** after the chosen policy; **B7a** can follow B1b/B3, **B7b** follows B6a.

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

## 7. Open questions

- Primary geometry format versus an optional `.vti` export.
- CUDA scope and whether GPU complex geometry is required with the CPU feature.
- Target fill ratio, peak memory (relative to the current block decomposition), and block-count cap.
- A second fixture: thin walls/multi-entity contact, plus a thick/bulk case for B5.
- 2D parity: the existing 2D path stays unchanged until a later branch.
