# Link boundary conditions (B4): stage 2

Status: merged into `dev` from `feat/link-bc` (2026-11-07): surface links + static IBB + one-rank
and 2-rank seams + fractional-wall, moving-wall and curved-wall (Taylor-Couette) accuracy; full
entity/BC binding, host-only CUDA and moving-wall refinement still open | Updated: 2026-11-08 |
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
  skipped. A cell on a periodic face keeps the links whose direction stays inside the domain and
  drops the ones that leave it along a periodic axis (`periodicAxes`); the periodic copy owns the
  latter. Skipping the whole cell used to drop the wall links on the face strip, because both ends
  of a periodic pair were skipped and the copy cannot supply a link that was never built.
- The manager is host-only and runs after streaming: it reads the pre-stream populations
  (`getPrevious`) and writes the opposite population of the wall-adjacent cell, so the next
  collision or communication sees the wall value. The CUDA path keeps the existing staircase
  bounce-back; B4 does not add a device link table.
- Moving walls use a per-entity wall-velocity map: the builder copies the entity velocity into
  the link, and the apply path adds the moving-wall term to the same Bouzidi branches.  A second
  table of `WallMotion{center, axis, omega}` covers rigid rotation, where the speed depends on the
  wall position (`u_w = omega x (x_wall - center)`); it is evaluated once per link because a rigid
  wall does not change in time.  Patches are bound to the task flags by `flbbc::bindPatches` (see
  the global plan), not by the link manager.

## 2. Results

Helical duct (`tests/feature/link_bc`, 20 steps, one rank):

| Item | Value |
| --- | --- |
| Links | 7160 |
| Fractional (`q != 0.5`) | 5190 |
| Fallback (`q = 0.5`, normal did not cross) | 925 |
| `q > 1` clamped to 1 (normal crosses beyond the link end) | 671 |
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
The centreline u_x(z) compared with the normalised parabola has L2 = 0.0115 (threshold 0.06);
the test also checks the two wall fractions and entities. This is an IBB wall-position baseline,
not a curved-surface analytic.

Moving-wall Couette (`tests/feature/link_bc_moving`): bottom entity 2 is static, top entity 3
moves along x with Uw = 0.01.  The measured u_x(z) against the linear Couette profile has
L2 = 0.0600 at 6000 steps (threshold 0.10, slope ratio 0.933); the run converges rather than
decaying: 48000 steps 0.0518, 150000 steps 0.0336 with slope 0.999, 300000 steps 0.0365.  The
long-run mass is bounded (3240 -> 3099 at 150k, -4.4%; a slow copy leak remains, see §3).

Curved-wall accuracy (`tests/feature/link_bc_curved`): Taylor-Couette between two concentric
cylinders, inner r = 5 rotating as a rigid body (entity 2), outer r = 10 at rest (entity 3),
periodic along z so the analytic profile has no end effects.  The two meshes are written as ASCII
STL at run time; the caps sit outside the domain.  Against `v_theta(r) = A r + B/r` at the
mid-plane: L2 = 0.0605 at 12000 steps (threshold 0.10), split as 0.0839 on the moving
inner half and 0.0386 on the static outer half; 6000 and 12000 steps agree to 0.5%.  This is the
B4 curved-wall baseline and it also exercises position-dependent wall velocity on a curved surface.

Periodic-face wall links: the link table keeps a periodic cell's wall links whose direction stays
inside the domain (`periodicAxes`) and drops only the ones that leave it along a periodic axis.
Previously both ends of a periodic pair were skipped, so the wall had no link on the whole face
strip and the copy could not supply one that was never built.  The link counts rise (accuracy and
moving 1960 -> 2432, curved 2144 -> 2976).  With the periodic copy fixed below, the plane,
moving and curved L2 are 0.0106 / 0.0092 (0.017 at 150k) / 0.0276, against 0.0118 / 0.0948 /
0.0746 before either change.

IBB mass correction: the interpolated reflection changes the boundary cell's mass.  On a complete
flat wall the symmetric links cancel; a strip truncated by a periodic face leaks, and the moving
wall term adds to it.  Apply now removes the per-cell change with the zero-momentum mode
(`f_k -= dM w_k`, `sum w = 1`, `sum w c = 0`), so momentum is untouched.  The moving run then
converges instead of decaying (before the correction: 96000/300000 steps 0.19/0.33).  With the
periodic copy fixed below, its mass creeps up 1.7% by 150k; the change is carried by the
wall-normal directions of the moving top wall, not by the exchange.

Periodic exchange: `Stream()` pulls `f_i(x) <- f_i(x - c_i)`, so a copy must write the cell one
voxel outside the face, in the block halo; writing the face cell itself landed the image one cell
inside the domain.  `Setup` pairs each recv through the send index it records
(`LatticeComm::SendCommIndex`) instead of through the block id, which returned the first send when
one block held an x and a y pair.  A corner column belongs to both pairs, so its diagonal
populations come from a corner pair with a combined mask; `getPeriodicPopDir` takes the
intersection over the masked axes, because `getCommPopDir`'s union would overwrite the
axis-aligned directions with the diagonal partner's values.  `feature/periodic_copy` asserts the
halo image and the mass/momentum conservation of a fully periodic 2D box, and a revert to the
face-cell write fails 103 of its 104 checks.  The stream stage's per-direction mass balance is
exact (1e-13) with the copy in place; without the corner pairs it is off by 5e-4 per step, which
the moving wall's opposite drift used to mask.

Gates: `make test` PASS 53 / FAIL 0 / SKIP 0 on the branch (includes `feature/link_bc`,
`feature/link_bc_accuracy`, `feature/link_bc_moving`, `feature/link_bc_curved`,
`feature/periodic_copy`, `block/link_bc_seam` and the CUDA targets on the RTX 3060);
`tools/check_headers.sh` exits 0 (the one baseline `utils/cuda_device.h` entry remains);
`make test-casegen` PASS 14 / FAIL 0.

## 3. Open items

- **Moving-wall accuracy.** The flat Couette L2 is 0.9% at the 6000-step baseline and 1.7% at 150k
  steps (slope 0.994); the curved 6k L2 is 2.8% (inner 4.2%, outer 1.2%).  The startup transient is
  gone; what remains is the moving-wall mass term below.
- **Moving-wall mass at the periodic strips.** With the exchange fixed, the drift at 150k is +0.05%
  in the static body-force channel and +1.7% in the moving Couette case.  The copy's own balance is
  exact (per-direction mass change of the stream stage: 1e-13; `feature/periodic_copy`), and the
  moving run's change is carried by the wall-normal directions at the top wall, so the remaining
  source is the moving-wall term where a periodic face truncates the link set, not the exchange.
  The per-cell IBB correction removes the per-cell sum but leaves the stream-side wall exchange.
- **Periodic exchange usage.** A caller whose domain is periodic along two axes must add the four
  corner pairs next to the two face pairs (`Setup` with a combined mask such as `XN|YN`), or the
  corner halos keep whatever they held.  The link_bc tests and `feature/periodic_copy` do so;
  `examples/couette3d` is periodic in x/y and still uses the face pairs only, so its corners stay
  stale until it is ported.  The manager now rejects the two wiring mistakes at setup: a face box
  that is not the block's boundary layer (the copy would overwrite a real neighbour cell) and a
  recv/send pair whose cell lists differ (the copy pairs them by position).
- **Fallback rate.** 925/7160 duct links still fall back to the half-way position because the band
  normal does not cross the link.  A further 671 links with `q > 1` are now clamped to `q = 1`: the
  void neighbour proves the wall is before the link end, so the clamp is geometric rather than a
  guess.  This reduced the fallback count from 1596, kept the duct mass exact and moved the
  curved-wall L2 from 0.07457 to 0.07453 (unchanged to two digits).  The 925 no-cross links are the
  remaining risk.
- **CUDA scope.** Host only. The device mirror would need a per-cell link table and a device
  apply path; B4 does not promise it.
