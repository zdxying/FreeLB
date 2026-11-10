# Storage comparison (B5): decision memo

Status: B5 decision landed on `dev` (2026-11-07); dense remains the default, B6a deferred | Updated: 2026-11-07 | Target: choose the storage policy after B3

> Exit condition: a later case is measured as memory-bound, or the thin fixture must meet the
> 0.25 owned-fill floor. Reopen B6a with this memo's baseline; otherwise this memo is the record.

## 1. Decision

**Dense AABB remains the default. B6a (active-slot storage) is deferred.**

Active-slot is the first candidate if a production case becomes memory-bound. Brick-sparse is not
chosen for the first storage branch: its seam traffic is 1.4-3.3x dense on these fixtures, and the
upstream halo exchange it needs is the hard part, not the storage itself.

The thin duct's 0.25 owned-fill floor is still not met (B3 reached 0.210); that is a block-shape
limit, not a storage-layout limit, so a storage rewrite would not fix it by itself.

## 2. Evidence

Command (one rank, serial):

```
make -C benchmarks/geo_storage
./benchmarks/geo_storage/geo_storage.exe
```

Machine: Intel i7-11800H (16 threads), g++ 13.3.0, `-O3 -march=native`, 2 warm-up steps and 50
timed steps. `dense` and `compact` are the same fill-aware table (B3, `max_blocks = 24`);
`compact` runs the existing `VoxelMap` loop pruning after `InitVoxelMap`, while the field storage
is still dense. `slot` and `brick` are first-order byte models, not solvers.

The `thin` fixture here is the two-plate `geo_thin`, not the helical duct; the duct's 0.210 fill is
the B3 number quoted in §1. The model columns exclude the flag field, boundary metadata and block
objects, so the true slot ratio is worse than the quoted one. No storage alternative was measured
for process RSS: the evidence is field bytes and MLUPS, and the budget's RSS row stays open.

| Fixture | Mode | Blocks | Allocated | Owned | Fill | Seam bytes/step | Field bytes | MLUPS |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bulk | dense | 24 | 37056 | 18944 | 0.511 | 0.90 MB | 6.82 MB | 17.85 |
| bulk | compact | 24 | 37056 | 18944 | 0.511 | 0.90 MB | 6.82 MB | 18.83 |
| thick | dense | 8 | 46656 | 30720 | 0.658 | 0.99 MB | 8.58 MB | 21.46 |
| thick | compact | 8 | 46656 | 30720 | 0.658 | 0.91 MB | 8.58 MB | 21.49 |
| thin | dense | 8 | 46656 | 31616 | 0.678 | 0.99 MB | 8.58 MB | 22.19 |
| thin | compact | 8 | 46656 | 31616 | 0.678 | 0.96 MB | 8.58 MB | 21.75 |

Storage models, relative to the dense field bytes:

| Fixture | Dense | Active-slot | Slot ratio | Brick 4^3 | Brick4 ratio | Brick4 seam | Brick 8^3 | Brick8 ratio | Brick8 seam |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bulk | 6.82 MB | 4.95 MB | 0.73x | 3.49 MB | 0.51x | 1.43 MB (1.59x) | 5.28 MB | 0.77x | 1.05 MB (1.17x) |
| thick | 8.58 MB | 7.27 MB | 0.85x | 5.65 MB | 0.66x | 2.92 MB (2.94x) | 6.03 MB | 0.70x | 1.40 MB (1.41x) |
| thin | 8.58 MB | 7.51 MB | 0.875x | 6.03 MB | 0.70x | 3.27 MB (3.29x) | 6.03 MB | 0.70x | 1.40 MB (1.41x) |

The active-slot model stores the fields per active slot (`VoxNum`, halo included), plus the
forward `int32` AABB-to-slot map and the reverse `size_t`/owned arrays. It excludes the flag
field, boundary metadata, and the block objects. The brick model aligns bricks to the grid, keeps
every brick that contains an active cell, and counts a seam face between two kept bricks; it
assumes the upstream halo exists. Seam bytes/step is the first-layer `Comm`/`MPIComm` cell count
times `q * sizeof(T)`; the all-layer communicator and the packed plan are not included.

## 3. Why dense

- The compact iteration is not faster: 0.98-1.06x dense MLUPS, within run-to-run noise. It prunes
  loops but does not change storage, so it cannot justify B6a on its own.
- Active-slot saves 12-27% of the field bytes in the model. On the deciding thin fixture it is
  0.875x, only marginally inside the B5 target of 0.90x, and the estimate excludes metadata.
- Brick storage saves 30-49%, but pays 1.4-3.3x seam traffic. At 4^3 the seam is 1.6-3.3x; the
  8^3 variant is cheaper to communicate but loses fill on the bulk shell (0.66).
- Dense is already tested, bit-identical for a fixed decomposition, and has no upstream-halo cost.
  The B3 thin limit is the AABB block shape, not the field layout, so a storage branch would not
  move the thin fill number.

## 4. Rejected alternatives

- **Brick-sparse as the first storage step.** Rejected for now: the measured seam traffic is
  1.4-3.3x dense and the upstream halo is the missing piece. The 4^3 brick only looks good on the
  storage column.
- **Active-slot as a drop-in.** Rejected: it is a field-layout rewrite (B6a), and the compact-loop
  measurement shows no speed headroom that would pay for it. Keep it as the first candidate when a
  case is memory-bound.
- **Changing the B3 thin gate to make B5 look successful.** Rejected: B3's thin result is 0.210
  against a 0.25 floor; the gate stays red until a storage/decomposition change actually moves it.

## 5. Next

- B4 landed; its refinements are tracked in `docs/COMPLEX_GEOMETRY_PLAN.md`.
- B7a (checkpoint layout/origin/digest) is next and can follow B1b/B3.
- The benchmark and this memo are on `dev`; the table reproduced there on 2026-11-07 with MLUPS
  inside run-to-run noise. Rerun `benchmarks/geo_storage` before reopening B6a.
- Reopen B6a only with a memory-bound production case or a thin-case requirement that dense cannot
  meet; start from `benchmarks/geo_storage` and this table.
