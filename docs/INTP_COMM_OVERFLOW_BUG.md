# Known communication defects (Intp out-of-bounds, MPI send/recv layout, cazsAref2d startup, same-level partial not refreshed)

Status: the exit condition is the §17.4 fix, the §17.7 handover to `examples/cazsAref2d` maintenance, and closing §18.8.3/§18.8.4, at which point this document is deleted; the usable-range criterion for same-level partial is folded into §18.9 (voided once 3D partial/full bit-identical measurement lands or the partial path is deleted); the local Intp out-of-bounds (§1) and the MPI Intp/Aver send/recv layout (§16/§17) are fixed, the rest are not | Updated: 2026-11-09 | Target: keep only criteria, invariants, guard mappings and unfixed defects

## 0. Conclusions that still hold

| Conclusion | Status | Details |
| --- | --- | --- |
| The Intp template's `+1` row falls outside the block; the local path writes past the heap | **Fixed**: index clamped into the block | §1, `tests/block/comm_index` |
| The MPI Intp send and receive sides enumerate different coarse-cell sets, so the lengths disagree | **Fixed**: `intp::BuildStencil` is the single enumeration | §16, `intp-count` in `amr_comm` |
| Aver's two sides intersect in different blocks' coordinate frames (`recvs=10 sends=72`) | **Fixed**: `aver::BuildAverStencil` | §17, `cazsAref2d` OMP build |
| `MPI_ERR_TRUNCATE` after the `cazsAref2d` MPI output step | **Not fixed**, not localised | §17.4 |
| Same-level partial refreshes only the direction set from `getCommPopDir()` | **Not fixed**: the 2D magnitude depends on main-loop order, 3D untested, forbidden under AMR; criterion and measurements in §18.9 | §18 |
| On AMR geometry partial's `DirSends/DirRecvs` are out of bounds (counts vary with geometry) | **Not fixed** | §18.8.3 |
| `ApplyCellDynamicsWithHiddenComm` has no `count` overload, so sub-cycled AMR steps cannot use it | **Not fixed**: first step `ndiff=116901` against the level-filtered plain path | §18.8.4 |
| Hidden steps do no halo exchange in non-MPI builds | Known limitation, not a defect | §18.8.4 |

## 0.0 Defect list and guards

| Defect | Scope | Status | Guard |
| --- | --- | --- | --- |
| Intp index out of bounds (local) | 2D/3D geometry | **Fixed** (§1) | `tests/block/comm_index` |
| `cavref2d` heap corruption | local Intp | **Fixed** (§1) | same; `amr_comm` regression 575 items |
| MPI Intp send/recv layout inconsistent | MPI path | **Fixed** (§16) | `intp-count` in `amr_comm`; `comm_index` static check |
| `cazsAref2d` never inits MPI | example | **Fixed** (§15.1) | no automated test |
| `cazsblock2d` never inits MPI | example | **Fixed** (§17.5) | no automated test |
| `cazsAref2d` Aver send/recv length mismatch | example + `SharedComm` | **Fixed** (§17.1) | `aver-count` in `amr_comm` |
| `cazsAref2d` truncate after the MPI output step | MPI path | not fixed (§17.4) | none |
| `cazsAref2d` intermittent failure after refinement | runtime refinement path | **Not fixed**; owner: `examples/cazsAref2d` (§17.7) | none |
| Same-level partial communication misses halo refresh | 2D/3D, reproducible on a single rank | not fixed (§18, criterion in §18.9) | `tests/block/normal_comm_correctness` (reports only by default) |

## 0.1 Commit index

| Document section | Commit | What changed | Guard |
| --- | --- | --- | --- |
| §1 | `c4b2a7e` | the Intp template's `+1` clamped into the block (2D/3D) | `tests/block/comm_index` |
| §2 | `6723677` | `flbcomm::checkCellIndex()`: report and `exit(1)` on out-of-bounds | the out-of-bounds path itself has no test |
| §14 / §16 | `05c9adf` | `intp::BuildStencil` becomes the single enumeration for cross-level Intp | `intp-count` in `amr_comm` |
| §15.1 | `a438be2` | `cazsAref2d` gains `mpi().init()` | no automated test |
| §15.2 / §17 | `05c9adf` | `aver::BuildAverStencil` unifies the Aver send/recv enumeration | `aver-count` in `amr_comm`; `cazsAref2d` OMP build |
| §17.5 | `45283e0` | MPI wiring and shutdown for `cazsblock2d` | no automated test |

§17.4, §17.7, §18.8.3 and §18.8.4 are not in the table: they are unfixed items with no corresponding commit.

## 1. Local Intp out of bounds (fixed)

`AddtoSharedIntpComm()` (2D `block_geometry2d.hh`, 3D isomorphic) takes a 2×2 (3D 2×2×2) neighbourhood on the fine block for every coarse cell; the template's `+1` row unconditionally steps one cell towards +y/+z; when the intersection touches the fine block's last row/column it lands outside the block, and `popIntpCommunicate()` writes `CellPop` by index without checking.

**Fix/invariant**: the construction side clamps out-of-range coordinates into the block (full, including halo), so the boundary template degenerates to taking the edge cell; it is not skipped on the consumer side (writing one fine cell less is a numerical change, and cannot be distinguished from "the geometry does not require this cell anyway"). Local, MPI send and receive, 2D and 3D all come out of `intp::BuildStencil`. The guard is `tests/block/comm_index`: 2D/3D all geometries `badSrc=0 / badDst=0`.

## 2. `flbcomm::checkCellIndex()`: report and `exit(1)`

`src/parallel/comm_plan.h`. An out-of-bounds index means the geometry side enumerated something wrong; silently skipping it would hide the bug in the numbers, and could not be distinguished from "the geometry does not need this cell anyway". The 11 call sites are in `block_lattice.hh` (`unpackCrossRecv`, `averCommunicate`, `popIntpCommunicate`) and `field_struct.hh` (`MPIAverCommunicate`, `MPIIntpCommunicate`, `AllMPIAverCommunicate`, `AllMPIIntpCommunicate`, two overload sets each); the out-of-bounds path itself has no automated test.

## 3. Scope: only Intp bypasses the geometry library

Of the three classes — same-level `Comm.Comms`, `AverComm`, `IntpComm` — the same-level and Aver sides both use the same AABB intersection, with boundaries guaranteed by the geometry; only Intp hand-writes `+1` to make up 4/8 source points. The failure condition is the fine block's base being clipped by a domain boundary and the intersection landing on the fine block's last row/column. A probe must use `geo.getBlock(findBlockIndex(id))`, not `geo.getBlock(id)` (which takes a block by container index and gives false positives when block order and id disagree).

## 14. MPI Intp exchange: the old enumeration disagreed between the two sides (superseded by §16)

The old implementation enumerated the receive side by "fine block self ∩ coarse block base" (65 coarse cells) and the send side by "fine block base ∩ coarse block extended by 2 layers" (66; `ext2` denotes the old code's `getExtBlock(2)` extended block, now deleted), so the `Cells` elements meant different things on the two sides (coarse vs fine) and the lengths naturally disagreed; there were also block pairs with a list on one side only, where the send side's `startC - 1` underflowed into a huge `size_t`. Three alignment variants were tried and rejected (drop the send-side offset, switch the send side to `nbaseblock`, switch the send side to the receive side's `getCellIdx`); all have been deleted — do not repeat them. The semantic decision is in §16.2.

## 15. `cazsAref2d` startup blockage (fixed)

### 15.1 `LoadBalancing()` segfault: `mpi().init()` never called

`main()` takes no arguments and never calls `mpi().init()`, so `MpiManager::numTasks` keeps an uninitialised value (0 in the debugger), `LoadBalancing(int ProcessNum = mpi().getSize())` therefore does `BlockIndex.resize(0)`, and the following `BlockIndex[0].push_back(...)` accesses an empty vector. The fix is `main(int argc, char* argv[])` + `mpi().init(&argc, &argv)` (`cazsArefblock2d.cpp:193`), consistent with `cavref2d`/`dambreak2d`; `cazsblock2d`'s identical problem is fixed together in §17.5. No automated test.

### 15.2 `AddtoSharedAverComm` send/recv length mismatch

During geometry construction `setRecvSendIdx` reports `recvs=10 sends=72 ratio=7.2`: the two `getCellIdx` calls are used in different blocks' coordinate frames (the coarsened virtual block `nCblk`'s frame vs the `self`-scale AABB), and 18 vs 10 has no 1:4 relation; the pathological source is the large shape difference between the fine block 22×62 and the coarse block 42×42. The fix is in §17.1.

### 15.3 One wrong approach (do not repeat)

Deriving `VSends` from `Recvs` (`loc / 2` mapped into `nCblk`): it runs 200 steps and prints "Calculation Complete!", but all 10 mappings fall outside `nCblk`. "It runs" is not "it is correct"; index validity must be verified separately (`tests/block/comm_index`).

## 16. MPI Intp exchange (fixed)

### 16.1 Reproduction

Run `amr_comm`'s `intp-count` from `tests/block` on two ranks: a single rank passes, two ranks crash or hang (ASan reports an out-of-bounds in `unpackCrossRecv`; at the time `MPIInterpolateCommunicate` also had an `MPI_Barrier`, so without ASan it stops there waiting for a message — that barrier was removed in `4430d78`). The probe shows, besides the length mismatch, block pairs with a list on one side only (e.g. `3→2`, `4→2` receive 0).

### 16.2 Semantic decision: the template may have no source outside the coarse block boundary

Coarse sources are always clamped to the coarse block's full range (including halo), degenerating at the boundary to taking the edge cell. This matches the semantics the local path already accepts — the single-rank reference is this semantics, and MPI must agree with it. The extension scheme (the template keeping a full 2×2 support at the boundary) was rejected: it would require the local path to go back to an extended enumeration and read across blocks, while "clamping" is already proven measurable and explainable by `comm_index`.

### 16.3 One enumeration definition, three layouts

`intp::BuildStencil<D>()` (`:104`) in `src/geometry/intp_stencil.h` is the single official definition of "which coarse cells take part in Intp", and local and MPI, 2D and 3D all come out of it:

- support domain = `intersection(fine full block, coarse base block)`, counted in coarse cells;
- each coarse cell expands `NumStencils = 2^D` stencils in the row order of `getIntpWeight<D>()` (2D 4 / 3D 8), each reading the clamped 2×2 (2×2×2) coarse neighbourhood and writing one fine cell;
- `SrcNum = intp::SrcNum<D>()` (`:59`) is the coarse source count per stencil: 2D 4, 3D 8; `send == recv × SrcNum` is guaranteed by construction;
- the clamping boundary is the **full block (including halo)**.

6 call sites are wired (2D/3D × local/`InitMPIIntpComm` send/receive), with only three layouts:

| Call site | Layout |
| --- | --- |
| `AddtoSharedIntpComm` (local) | per stencil `[SrcNum coarse sources, 1 fine target]` |
| `InitMPIIntpComm` send side | per coarse cell `NumStencils × SrcNum` coarse sources |
| `InitMPIIntpComm` receive side | 1 fine target per stencil |

Call sites: `AddtoSharedIntpComm` / `InitMPIIntpComm` (one set each for 2D/3D). The old `getExtBlock(2)`, `startC - 1` and `getCoasenedBlock` enumerations are all deleted. Guard: `unpackCrossRecv` (receive unpack) calls `flbcomm::checkCellIndex`; send-side indices are generated clamped by `BuildStencil`.

### 16.4 Verification

| Verification | Result |
| --- | --- |
| `BuildStencil` vs the hand-written enumeration before the fix (`amr_comm` geometry, bit-identical) | 4 block pairs × 5120 entries equal — local indices unchanged |
| `tests/block/comm_index` (11 local geometries + MPI Intp static check) | 44 items all pass |
| `tests/block/amr_comm` including `intp-count` | 575 items all green (this path never ran before the fix) |
| Two-rank identity `send == recv*SrcNum` | holds for all 4 block pairs (including `3→2`, `4→2`, which were 0 before the fix) |

Limitation: `amr_comm` is 2D only; 3D Intp has only `comm_index`'s static check, with no numerical gate.

### 16.5 Invariants (for later maintenance)

1. Any change to `intp::BuildStencil` is a **numerical change**: the local single-rank reference and two-rank MPI change together, and `amr_comm` cannot tell "changed correctly together" from "changed incorrectly together"; take a per-configuration checksum per AGENTS.md §4.
2. A new cross-level communication path must reuse `intp::BuildStencil` and must not hand-write an enumeration; `comm_index`'s MPI check audits every block pair against "recv == send × SrcNum + indices in range".
3. The same class of problem (Aver's `SharedComm` length mismatch) is fixed the same way — already done, see §17.

## 17. Unified Aver send/receive + parallelism for two CA cases (fixed)

### 17.1 Fix for the local Aver length mismatch (§15.2)

`aver::BuildAverStencil<D>(CBase, FBase, coarseOverlap, fineOverlap, ...)` in `src/geometry/aver_stencil.h` is isomorphic to §16:

- support domain = `intersection(coarse FULL block, fine base block)`, counted in coarse cells; the coarse side includes the overlap ring facing the fine block (that ring is what the exchange writes);
- each coarse cell's fine sources = the 2×2 (2×2×2) neighbourhood of that coarse cell in the coarsened fine frame (anchor `FBase.min - 0.5*Fvox`, spacing `2*Fvox`), `loc*2 + off` then **clamped to the refined grid boundary** (`aver_stencil.h:158-167`);
- `SrcNum = aver::SrcNum<D>()` (2D 4, 3D 8); the three layouts (local `[SrcNum fine sources,1 coarse target]×N`, MPI Sends, MPI Recvs) come out of the same function, and `send == recv * SrcNum` is guaranteed by construction.

6 sites are wired: `AddtoSharedAverComm` (2D/3D local) and `InitMPIAverComm` send/receive (2D/3D), call sites `block_geometry2d.hh:792,809,914`, `block_geometry3d.hh:898,915,1020`; the old two-frame enumeration ("coarsened block + refined frame": `getCoasenedBlock` + `getRefinedCellIdx`) is deleted. On regular geometry (`amr_comm`) the old and new local enumerations are bit-identical (4 pairs, 320 entries). Under the new enumeration `cazsAref2d`'s pathological geometry satisfies the identity for all 20 block pairs (`bad=0`), and OMP can also get past geometry construction into the main loop — but "running 50000 steps to completion" is **intermittent** (§17.7).

**Invariant**: the support domain stops at the fine base, but the neighbourhood sampling `loc*2 + off` is clamped into the refined grid `nCmesh*2 - 1` (`aver_stencil.h:158-167`), and that grid covers the base plus one cell, so **a fine source slot landing outside the base (in the fine halo) is a normal result**: measured, each coarse block has 2–65 fine source slots in the halo. This is not a defect, but any reasoning of the form "fine sources are always inside the base" is wrong.

### 17.2 overlap is not a level rule

- `lv != 0 ? 2 : 1` appears only in the direct-construction path and in `Init(GeoHelper)` (`block_geometry2d.hh:128-131,186`); on the helper path (every case in this repository) all blocks have overlap = the helper's `_Overlap` = 1. Passing a level rule would break `amr_comm`'s aver/intp (17028 out of bounds); do not do it.
- overlap must be **this block's** `Block2D::getOverlap()`; the two sides of the Aver and Intp MPI branches each pass their own block's real value. That is why Intp's `BuildStencil` has a dual-overlap signature (one for fine, one for coarse).
- The same "level rule vs real value" mismatch appears again in the partial `Dir*` indices of §18.8.3, and that time there was no rollback — it is an unfixed defect.

### 17.3 Refined neighbourhood clamping (the root cause of 17028, fixed)

Aver's refined neighbourhood `loc*2 + off` ran past the refined grid when the support domain touched the edge (17028 was an out-of-range index); it is now clamped to `nCmesh*2 - 1`. For interior stencils this is a no-op (a bit-identical comparison probe confirms regular geometry is unchanged).

### 17.4 Remaining `cazsAref2d` MPI problem (not fixed)

After the fix, two-rank MPI gets past the old dead point and through the output step, with Solid% of the same order as OMP, but the field-level Aver/Intp exchange after the output step (`BlockFieldManager::runCrossPlan`) reports `MPI_ERR_TRUNCATE`. Probes verify that all 20 Aver/Intp block pairs satisfy the send/recv count identity; the truncate comes from one link between the two sides' active subsets and the posted requests under the count filter (`activeDepth`), not yet localised. This is the last known problem of MPI + AMR + this case.

### 17.5 MPI parallelism for `cazsblock2d` (fixed in passing)

`main()` gains `mpi().init` (as in §15.1), and geometry construction moves from direct `BlockGeometry2D(Ni, Nj, Thread_Num, ...)` (no rank distribution; under MPI every rank's block set is empty → runs on nothing) to Helper + `AdaptiveOptimization + LoadBalancing`. Measured: OMP 2000 steps normal; MPI np=2 completes all 2000 steps. Known difference: a different CA capture sequence at block seams causes a State difference in 16/91204 cells (Conc differs by 0.19 at the same position); this is the serial assumption of the CA algorithm on AMR block seams, not a communication error.

### 17.6 Verification summary

| Verification | Result |
| --- | --- |
| `aver::BuildAverStencil` vs the old local enumeration (regular geometry, bit-identical) | 4 pairs, 320 entries equal |
| `cazsAref2d` geometry, 20 block pairs send/recv identity (probe, np=2) | bad=0 |
| `amr_comm` aver-count / intp-count (np=2 vs single rank) | bit-identical |
| `cazsAref2d` OMP 50000 steps | completes intermittently (§17.7); the old code could not even build the geometry |
| `cazsblock2d` OMP 2000 steps / MPI np=2 2000 steps | both complete fully (the 16-cell seam difference is recorded) |

Limitation: 3D Aver/Intp still has only `comm_index`'s static check, with no numerical gate; `cazsAref2d`'s OMP run only proves it can start — for the numerical problems see §17.7; the MPI truncate is not fixed (§17.4).

### 17.7 Two pre-existing defects in runtime refinement (not fixed; owner: `examples/cazsAref2d`)

The owner is the AMR runtime refinement path of `examples/cazsAref2d`, unrelated to communication; `cazsAref2d` therefore cannot serve as a gate:

- **16 threads, nondeterministic**: with `Thread_Num=16` and `TotalStep=50000`, only some runs complete, and every failure is after refinement (step 1000): `BlockZhuStefanescu2D::UpdateDeltaFs()` reports `deltaf > 1` (`zhu_stefanescu2d.hh:125`) or `Ceq > 1`, the two blocks' output interleaves on stderr (`Apply_SimpleCapture` is `#pragma omp parallel for` over `BlockZS`, `zhu_stefanescu2d.h:423-424`), and two runs of the same binary take different branches — a race. With `RefineCheckStep` set to 1e9 (dynamic refinement off, the initial two-level geometry kept) 3/3 pass: the main loop and the initial cross-level communication are clean; the problem is "runtime refinement + re-initialisation".
- **≤8 threads, always fails and deterministically**: 1/2/4 threads, 6 runs, all exit at the same place with `[ERROR] AABB<T, D>: min > max`, from `DynamicBlockLatticeHelper2D::PopFieldInit()` (`block_lattice.hh:2081`) — it indexes the old block list with **new** block numbers: `for (iblock < BlockGeo.getBlockNum())` followed by `BlockGeoHelper.getAllOldBasicBlock(iblock)` (`:2104-2105`); the probe measured `n_old = 12`, `n_new = 22` on the out-of-bounds run, and taking the AABB of the `BasicBlock` read out of bounds then calls `exit(1)`. At 16 threads the two refinements happen to have the same block count (32/32), which sidesteps it.

## 18. Same-level partial communication: only inflow directions refreshed (not fixed)

`NormalCommunicate()` (`block_lattice.hh:1816`) refreshes, for every halo cell, only the direction set D given by `getCommPopDir()` (`src/parallel/communicator.h:64`); `normalFullCommunicate()` (`:1835`) refreshes all q. The two differ only in the halo's **unrefreshed slots D̄** (all q minus D). The defect is unrelated to AMR, reproduces at `np=1`, and exists in both 2D and 3D; whether it is acceptable depends on main-loop order, dimensionality and ghost readers — the criterion, measurements, rejected subsets and the cost of full are in §18.9, and only the defect itself and the open items are kept here.

Affected call sites: `cazsArefblock2d.cpp:398,406`, `openbd2d.cpp:189`, `cavity2d.cpp:253`, `cavity3d.cpp:257` (hidden defaults to Partial, with full used at the output step and the final checkpoint: `cavity2d.cpp:266,284`, `cavity3d.cpp:270,288`); `couette2d.cpp:215` and `cavref2d.cpp:228` use full and are unaffected; `cazsblock2d`'s field path (`cazsblock2d.cpp:424`) is also unaffected.

Guard: `tests/block/normal_comm_correctness` compares the inner checksums of partial and full under the same decomposition; by default it only reports, and `--strict` is the pass condition. AGENTS.md §3 records it as the acceptance condition for this known gap. It runs single-process only, and `tests/feature/hidden_comm` only shows hidden agrees with plain (both wrong the same way), so **there is no evidence for partial vs full on the MPI side**.

### 18.8 Open items for AMR and hidden

#### 18.8.3 partial's `Dir*` indices out of bounds on AMR geometry (not fixed)

On refined geometry, partial's `DirSends` / `DirRecvs` indices are built in the helper's `overlap = (level != 0) ? 2 : 1` frame (`block_geometry2d.hh:131`, same rule at `:186`, `:586`, `:985` in the same file), while the running block uses `overlap = _overlap = helper.getOverlap()` = 1 (`:128-131` in the same file, `useHelperOlap = true`). The defect is per-configuration identical on the merge baseline `ff01afd`; it was not introduced by this branch.

Measured on the 8 geometries of `tests/block/comm_index` at np=4 (`N` = block cell count; the counts vary with geometry; this table is authoritative — the counts from earlier uncommitted probes are not reproducible and are void):

| Configuration | send ≥ N | recv ≥ N | packed cells outside `_collideFirst` |
| --- | --- | --- | --- |
| lv2 inset | 39 | 70 | 148 |
| lv2 at corner | 6 | 6 | 55 |
| lv3 inset | 41 | 69 | 145 |
| uniform / ring refined | 0 | 0 | 0 |

Invariant: the out-of-range indices are read by `packNormalSend` (`block_lattice.hh:215`) and written by `unpackNormalRecv` (`:231`); the packed cells are not in `_collideFirst`, i.e. **pack first, collide after**, so what is sent is the pre-collision value; on the same geometry `MPIComm`'s count is 0, so that path is unaffected. Rebuilding the running geometry with `useHelperOlap = false` makes the two sides' sets agree, pinning the cause on the overlap frame mismatch — the same class as §14 (Intp send/recv layout) and §15.2 (Aver send/recv length): the two sides do not enumerate the same set in the same frame.

#### 18.8.4 Boundaries of the hidden path (not fixed)

- `ApplyCellDynamicsWithHiddenComm` has only one overload (`block_lattice.h:502`), with no `count` / level filter parameter, so **sub-cycled AMR steps cannot use it**: measured against the level-filtered plain path, the first step already gives `ndiff = 116901, maxdiff = 0.0667`.
- In non-MPI builds a hidden step does no halo exchange at all; the `#else` at `block_lattice.hh:1555` falls straight back to `ApplyCellDynamics`. Know this when writing single-process comparisons.
- Under `_VOX_ENABLED` both half-passes are wrong: they use the raw `id` and `flagarr[id]` (`block_lattice.hh:202-204`, `:294-296`), whereas `ApplyCellDynamics` uses `map[id]` / `getVoxNum()` (`block_lattice.hh:133-150`).
- The comment at `block_lattice.hh:73` says `getPrevious()` reads the neighbour cell. It reads this cell's pre-stream value (`CyclicArray::getPrevious`, `field.h:716-726`), i.e. a **self-read**; do not treat `getPrevious()` in a BC as the source of "reading a ghost".

### 18.9 Usable-range criterion for partial

partial is not an approximation of full but a trick that relies on error propagation cancelling; the criterion in this section is voided once 3D partial/full bit-identical measurement lands, or the partial path is deleted.

**Why only the full q closes**: D is taken by `getCommPopDir()` as the deduplicated union of `c(i)[axis] > 0` per axis (the sign follows the neighbour direction): D2Q9 face 3 / corner 5, D3Q19 face 5 / edge 9 / corner 12. `opp(D) != D`, while a pure bounce-back ghost needs `opp(D) ⊆ D`; and BGK mixes all q into every output through rho/u. No proper subset closes, so the necessary and sufficient condition for bit-identical equality can only be: the output of every operation that reads ghost D̄ (collision, BC, statistics, hidden packing) is insensitive to it, and nothing else moves it.

**Usable conditions** (check each for a new case; hitting one forbids partial):

- single block (`Comm.Comms` empty): both versions are no-ops.
- 2D, same level, one layer of ghost, overlap=1, no AMR, and the skeleton must be `Collide → Stream → Communicate`; measured bit-identical for 2×2 over 200 steps and 3×3 over 40 steps.
- main-loop order sets the magnitude; consistency is not a property of partial: with the repository's current `Stream → Collide → Communicate` the two diverge by step 3–5, and `Collide → Communicate → Stream` is of the same magnitude.
- forbidden if the ghost is read after being collided: `anti_pressure` reads the neighbour VELOCITY, the residual statistics walk the whole array including ghosts, AMR's Intp reads the coarse block halo's pop+rho+u, and CUDA has no cell list and can only collide over all indices.
- hidden comm can only be consistent with `CommMode::Full`: partial's D̄ lands exactly in hidden's overlap window.
- 3D has no evidence of partial and full being bit-identical — forbidden; AMR forbidden (it additionally has the `Dir*` out-of-bounds of §18.8.3, unrelated to partial/full).

**Rejected subsets and the cost of full**: patching corner cells, patching only edges/faces, and enlarging the direction set all fail to close; the "0 difference" of dropping q0 covers one configuration only with an unexplained cause, and must not be used as a bandwidth argument; one full pass at initialisation does not change the leakage (it happens at step 2 and reaches base at step 4). Cost of full (np=2, single node, 16 cores, median of 3–5 runs): 2D 120²/10×10 +5.4% (payload ×2.80), 3D 48³/4×4×4 +37.3% (×3.23); 2D is dominated by latency and packing, 3D by bandwidth and packing.

**Unverified**: there is no recorded run of 3D partial and full being bit-identical; "seam ghost self-synchronisation" under `Collide → Stream → Communicate` is only an observation over 2×2 and 5 steps, not a theorem; for the 1×N domain boundary residual only the position is known, not the trigger; the A/B/C order comparison numbers come from uncommitted probes.
