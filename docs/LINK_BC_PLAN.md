# Link boundary conditions (B4): stage 2

Status: merged into `dev` from `feat/link-bc` (2026-11-07): surface links + static IBB + one-rank
and 2-rank seams + fractional-wall, moving-wall and curved-wall (Taylor-Couette) accuracy; full
entity/BC binding, host-only CUDA and moving-wall refinement still open | Updated: 2026-11-07 |
Target: replace staircase reflection with a surface-link boundary

> Exit condition: the B4 accuracy/seam/CUDA-scope baseline is recorded (done at the merge), and the
> open items in §3 are closed or reassigned. The measured values live in the tests; delete this
> plan once §3 is empty.

## 1. Decision

- The link table is **derived**, block-local data built from `VoxelGeometry` and the flag field.
  It reads the void flag in the block halo, so a wall that lies in a neighbouring block is still
  owned by the cell that sees it; no new field and no extra communication are introduced.
- Static interpolated bounce-back uses the narrow-band signed distance/normal. For a link with
  fraction `q`:
  - `q < 1/2`: two-point Bouzidi formula using the fluid node behind the boundary cell;
  - `q >= 1/2` (and the thin/fallback case): one-sided formula using the opposite population.
  `q` is clamped to `1/2` when the node behind is outside the block or void, or when the band
  normal does not define a crossing. `Fallback` and `Thin` stay on the link for diagnostics.
  A boundary cell with no band sample (a domain-edge wall) gets a half-way link instead of being
  skipped. A cell on a periodic face is skipped by the link table; the periodic manager refreshes
  it, so the face is not treated as a wall.
- The manager is host-only and runs after streaming: it reads the pre-stream populations
  (`getPrevious`) and writes the opposite population of the wall-adjacent cell, so the next
  collision or communication sees the wall value. The CUDA path keeps the existing staircase
  bounce-back; B4 does not add a device link table.
- Moving walls use a per-entity wall-velocity map: the builder copies the entity velocity into
  the link, and the apply path adds the moving-wall term to the same Bouzidi branches.  A second
  table of `WallMotion{center, axis, omega}` covers rigid rotation, where the speed depends on the
  wall position (`u_w = omega x (x_wall - center)`); it is evaluated once per link because a rigid
  wall does not change in time.  Patches still do not drive tasks or boundary conditions beyond
  these entity tables.

## 2. Results

Helical duct (`tests/feature/link_bc`, 20 steps, one rank):

| Item | Value |
| --- | --- |
| Links | 7160 |
| Fractional (`q != 0.5`) | 5190 |
| Fallback (`q = 0.5`, normal did not cross) | 1596 |
| Thin (node behind outside/void) | 374 |
| Links with an entity | 7160 |
| Mass, one block | 1294 -> 1294 (drift 0) |
| Mass, 8 blocks | 1294 -> 1294 (drift 0) |
| Links, 8 blocks | 7160 (same as one block) |

The 8-block run is the seam fixture: the link table sees the same owned links through the halo,
and `NormalFullCommunicate()` keeps the multi-block mass unchanged. The duct is curved and thin;
this is the thin-wall mass check.

Two-rank seam (`tests/block/link_bc_seam`, `mpirun -np 2`, 2 x 2 x 2 blocks, four blocks per
rank): 8 blocks, 7160 links and 1294 fluid cells in total, mass 1294 -> 1294 over 20 steps, and
both ranks own links (checked with an `MPI_MIN`). The link table is built from block halos whose
flags `writeVoxelFlags` fills from global positions, so the MPI halo needs no flag communication;
`InitVoxelMap` runs before the lattice is built. Global counts are reduced by the test itself
rather than `getTotalOwnedCellNum`, whose non-root result is only fixed by the unmerged B3 branch.

Fractional-wall accuracy (`tests/feature/link_bc_accuracy`): a 16 x 16 x 8 channel periodic in
x/y, with STL walls at q=0.25 (bottom) and q=0.75 (top) and a constant body force along x.
The centreline u_x(z) compared with the normalised parabola has L2 = 0.0118 (threshold 0.06);
the test also checks the two wall fractions and entities. This is an IBB wall-position baseline,
not a curved-surface analytic.

Moving-wall Couette (`tests/feature/link_bc_moving`): bottom entity 2 is static, top entity 3
moves along x with Uw = 0.01.  The measured u_x(z) against the linear Couette profile has
L2 = 0.0948 (threshold 0.10); the top q=0.75 one-sided branch is the weak part of the moving term.

Curved-wall accuracy (`tests/feature/link_bc_curved`): Taylor-Couette between two concentric
cylinders, inner r = 5 rotating as a rigid body (entity 2), outer r = 10 at rest (entity 3),
periodic along z so the analytic profile has no end effects.  The two meshes are written as ASCII
STL at run time; the caps sit outside the domain.  Against `v_theta(r) = A r + B/r` at the
mid-plane: L2 = 0.0746 relative to the wall speed (threshold 0.10), split as 0.0967 on the moving
inner half and 0.0561 on the static outer half.  6000 and 12000 steps agree to 1e-4, so the
profile is converged.  This is the B4 curved-wall baseline and it also exercises
position-dependent wall velocity on a curved surface.

Gates: `make test` PASS 45 / FAIL 0 / SKIP 0 on `dev` after the merge (includes `feature/link_bc`,
`feature/link_bc_accuracy`, `feature/link_bc_moving`, `feature/link_bc_curved`,
`block/link_bc_seam` and the CUDA targets on the RTX 3060); `tools/check_headers.sh` exits 0 (the
one baseline `utils/cuda_device.h` entry remains).

## 3. Open items

- **Moving-wall accuracy.** The per-entity tables and wall term are implemented; the flat Couette
  L2 is 9.5% and the curved inner half is 9.7%, with the q=0.75 one-sided branch the main error.
  Refine the diagonal-link/wall-term interpolation before using it for production moving walls.
- **Entity/patch to BC binding.** Entity -> wall velocity (constant and rigid rotation) is wired;
  patches still do not drive tasks or other boundary conditions.
- **Fallback rate.** 1596/7160 duct links fall back to the half-way position because the nearest
  band normal does not cross the link. The channel and Couette tests have periodic faces skipped
  and no fallback in the measured profile; the duct fallback remains a visible risk for curved
  walls with sharp features.
- **CUDA scope.** Host only. The device mirror would need a per-cell link table and a device
  apply path; B4 does not promise it.
