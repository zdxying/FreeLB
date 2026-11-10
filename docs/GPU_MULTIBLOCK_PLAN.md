# GPU multi-block support evaluation plan

Status: partial: M0 partially landed (M0-2 only) | Updated: 2026-11-09 | Target: prepare multi-block capability for AMR (cross-resolution) and multi-GPU, not to speed up a uniform grid.
Machine: RTX 3060 Laptop / sm_86 / 30 SM / 6 GB / CUDA 12.5 (WSL2) | Exit condition: the Go/No-Go gate in §7 gives a verdict; if No-Go, keep only the remaining M0 guard and the §8.2 active-index compaction. **§7 currently has no recorded verdict.**

## 0. Summary

- **Current state**: the GPU path runs only `BlockLats[0]`; with `block_num > 1` it silently computes only 1/N of the cells — a correctness hazard, not a performance issue, and it should be closed first.
- **Necessity**: only **AMR (cross-resolution)** and **multi-GPU** genuinely require multi-block. The two common justifications, VRAM and sparse geometry, **do not hold** once the arithmetic is done (see §2); do not use them as the basis for starting the work.
- **Performance cost**: the dominant term is kernel launch, not bandwidth and not wavefront quantisation. Measured marginal cost on this card under WSL is ~6.6 µs/launch; 2 launches per block per step, so 8 blocks cost about +47% time (throughput −32%) and 16 blocks about +103% (throughput −51%). **CUDA Graph can pull the N≤32 cost back to ≤15%**, and that is the key variable deciding whether to do this at all.
- **Recommendation**: do the remaining M0 guard first (zero risk); defer M1 until the two measurements in §8 and the active-index compaction are done, and start it only when the §7 Go condition holds; AMR/cross-resolution is a separate phase.

## 1. Current state (code facts)

| Location | Fact |
|---|---|
| `src/data_struct/block_lattice.hh:991-1018` | Comment `for now only one block is supported`; `CuDevStream` and both `CuDevApplyCellDynamics` overloads hard-code `BlockLats[0]` (`:995`, `:1004`, `:1015`); the multi-block loop is already written but commented out |
| `README.md` | `Cuda for multi-block structure is not supported for now` |
| `benchmarks/cavity3d_cu/cavity3d.ini` | `block_num = 1` (the GPU benchmark only measures a single block) |
| `src/data_struct/block_lattice.hh:401-402` | `CuDevStream` does one `<<<LatSet::q, 1>>>` launch per block (19 CTAs × 1 thread), only to rotate pointers |
| `src/data_struct/cuda_field.h` (`cudev::CyclicArray`) | In-block streaming is a **pointer rotation**, O(1) and independent of cell count; it does not apply across blocks or across resolutions |
| `src/data_struct/block_lattice.hh` (`NormalCommunicate`, `MPI*Communicate`) | Inter-block communication **all goes through host memory**; there is no device-side version |
| `src/data_struct/block_lattice_base.h:129-130`, `block_lattice.h:94-96` | One set of `dev_Delta_Index` + `dev_Fields` + `devOmega/dev_Omega/dev_fOmega` per block (3 single-element `cudaMalloc`, `block_lattice.h:115-120`) + `dev_BlockLat` |
| `src/data_struct/block_lattice.h:292`, `:306` | Both `CuDevApplyCellDynamicsKernel` variants already have `if (idx >= n) return` (`:298`, `:312`); the `Reg` variant is folded into the same kernel, the POP policy chosen by the `POPPOLICY` template parameter of `cudev::Cell` |
| `src/utils/field_checksum.h` | `PopChecksumDevice` accepts a single block only; the call site still uses `getBlockLat(0)` (`benchmarks/cavity3d_cu/cavity3d.cu`) |
| `src/data_struct/block_lattice.h:81-85` | Omega / _Omega / fOmega are already **per block** — the structure for a different τ at a different resolution is already in place |
| `src/geometry/block_geometry3d.h` (`getMaxLevel()`) | AMR multi-resolution blocks are first-class (on the CPU side) |

**An existing defect (closely tied to multi-block, fixed in M0)**

**Silent error**: with `block_num > 1` the GPU computes only block 0, the POP of every other block stays at its initial value, and the program still prints `Calculation Complete!` and MLUPS. This is the same class of problem as the sm_89 incident recorded in `benchmarks/cavity3d_cu/cavity3d.cu` (a failed launch is not checked). The M0-1 error guard has not landed yet (§4).

## 2. Necessity classification

### 2.1 Holds — cannot be done without multi-block

| Need | Why multi-block is required |
|---|---|
| **AMR / cross-resolution** | Different refinement levels = different Δx, different τ, different Δt. This information is per-block (`getMaxLevel()`, per-block Omega). A single block can only have one resolution, so **the GPU path and AMR are currently mutually exclusive**, and AMR is FreeLB's headline capability (the first two images in the README). |
| **Multi-GPU** | With one process/device per card, the block is the natural unit of distribution and load balancing. |

### 2.2 Does not hold — common justifications that fail the arithmetic (stated honestly)

| Common justification | Result of the arithmetic |
|---|---|
| **"A single block does not fit, VRAM is too small"** | Does not hold. Each cell is about 93 B (POP 19×4 + rho 4 + u 12 + flag 1); 6 GB of VRAM holds about 6.5e7 cells ≈ **400³**. A single block can fill the whole card's VRAM; VRAM is not a reason to split. |
| **"Complex geometry has lots of void, splitting can skip the idle threads"** | **Right direction, but it does not belong to multi-block.** The right fix for eliding void is **active-index compaction**: the CPU side already has the `GenericvectorManager` index table and the `ApplyCellDynamics(blockids)` overload, but the `BulkTaskIds` task-table build is commented out in the benchmark and the GPU does not consume a task table either (see §8.2). One block can do it too and the launch count stays 1; it is an optimisation **orthogonal** to multi-block and should be its own item. |
| **"Multi-block speeds up a uniform grid"** | The opposite. On a uniform dense grid on one card, multi-block adds only three pure costs: launch, halo and overlap (see §3). |

### 2.3 Conclusion

> **The target scenarios for multi-block are AMR and multi-GPU, not speeding up a uniform grid.**
> If the work is justified by "speed-up", by M1 you will find what you built is about 47%–103% slower than today (throughput −32%–−51%), and be stuck. Anchor the goal on "capability" first.

## 3. Performance cost analysis

### 3.1 Measurements (new in this round)

Probe: `benchmarks/multiblock_probe/gpu_multiblock_probe.cu` (an empty kernel measures pure launch; a DRAM-bound copy measures the real cost of splitting a fixed total workload into N launches; both are also run against CUDA Graph)

```
nvcc -O3 -arch=native benchmarks/multiblock_probe/gpu_multiblock_probe.cu -o /tmp/probe2 && /tmp/probe2   # repository root; wrap with wsl.exe under WSL2
```

**[A] Empty kernel, µs/step (500 reps + 20 warmup)**

| N launches | Per-launch | CUDA Graph | ns/launch | Graph speed-up |
|---|---|---|---|---|
| 1 | 7.68 | 7.73 | 7676 | 0.99x |
| 2 | 14.74 | 7.18 | 7370 | 2.05x |
| 4 | 31.01 | 7.50 | 7753 | 4.13x |
| 8 | 61.39 | 9.64 | 7674 | 6.37x |
| 16 | 122.54 | 16.16 | 7659 | 7.58x |
| 32 | 239.30 | 30.21 | 7478 | 7.92x |
| 64 | 487.54 | 60.11 | 7618 | 8.11x |
| 128 | 848.23 | 108.13 | 6627 | 7.84x |

Marginal cost: **6.6 µs/launch launched one by one; 0.79 µs/launch with Graph (8.4x difference)**.

> ⚠️ The absolute 7.6 µs/launch is high; WSL2's GPU passthrough inflates launch latency, native Linux is expected at 2–3 µs. **Take the absolute value from a native-Linux re-run, but the structural conclusion "marginal cost is linear + Graph compresses 8x" does not depend on the absolute value.**

**[B] DRAM-bound, fixed total workload (4M float, 32 MB/step, 308 GB/s at full speed) split into N launches**

| N launches | µs/step | GB/s | vs N=1 | Graph µs | Graph GB/s |
|---|---|---|---|---|---|
| 1 | 108.85 | 308.3 | 1.00x | 108.26 | 310.0 |
| 2 | 111.97 | 299.7 | 1.03x | 110.00 | 305.0 |
| 4 | 121.99 | 275.1 | 1.12x | 114.94 | 291.9 |
| 8 | 131.98 | 254.2 | 1.21x | 113.35 | 296.0 |
| 16 | 321.03 | 104.5 | **2.95x** | 121.99 | 276.9 |
| 32 | 575.43 | 58.3 | **5.29x** | 140.31 | 239.1 |
| 64 | 1334.98 | 25.1 | 12.26x | 175.78 | 190.9 |
| 128 | 986.70 | 34.0 | 9.06x | 245.52 | 136.7 |

How to read it:

- For N ≤ 8 the cost is mild (within 1.2x); **at N ≥ 16 it degrades sharply**; the degradation is **larger** than the pure launch sum of [A] (16 launches: [A] predicts +115 µs, measured +212 µs), which says that besides CPU-side enqueueing there is a GPU-side per-kernel cost (kernel teardown/front end); **the mechanism awaits nsys confirmation**, so do not extrapolate directly from [A];
- At N ≥ 32 the spread between two re-runs is about ±20%; numbers in that range are order-of-magnitude only;
- **Graph pulls N=16 from 2.95x back to 1.12x and N=32 from 5.29x back to 1.29x** — the only means in this project that turns multi-block from a loss into a wash.

### 3.2 Cost items

| Cost item | Magnitude | Varies with | Removable |
|---|---|---|---|
| **kernel launch** | measured marginal 6.6 µs/launch (WSL); 2 per block per step | ∝ block count | **CUDA Graph (→0.79 µs)**; then merge the N stream launches into 1 (grid = N×q), 2N → N+1 |
| Wavefront quantisation (tail) | `ceil(n/128)` CTAs rounded up to a multiple of the 150 resident CTAs; 125k cells → about 7% | significant for small blocks | a lower bound on block size |
| Insufficient occupancy | resident capacity = 30 SM × 5 CTA (96 regs × 128 thr = 12288 regs/CTA, 65536/12288 = 5) = **150 CTA = 19200 cells**; below that the SMs are underfed | small blocks | same as above |
| **overlap redundancy** | `(B+2)³/B³ − 1`: 50³ → 12.5%, 25³ → 26%, 16³ → 42% | worse the smaller the block | block-size lower bound; overlap is already a parameter |
| **halo exchange** | traffic O(N²ᐟ³) vs compute O(N), small in magnitude; **but going through D2H→CPU comm→H2D costs ~9 MB/step + two synchronisations ≈ ms, which destroys the entire gain** | surface-to-volume ratio | **a device-side halo kernel is mandatory**; the host `normalCommunicate` cannot be used directly |
| RegPop register residency | 96 regs / LOCAL:0 (D3Q19 FP32, measured); runtime dispatch by blockIdx breaks residency | implementation-dependent | one kernel per block template instantiation, with the block index as a **compile-time** argument, no runtime dispatch |
| VRAM/object overhead | about 5–6 `cudaMalloc` + H2D per block | ∝ block count | one large allocation plus pointer offsets; an initialisation-time cost, acceptable |

### 3.3 Scaled to the real solver (100³, FP32 2211 MLUPS ⇒ 452 µs/step)

Adding the **measured additional cost** of [B] (`t(N) − t(1)`, independent of workload) linearly; 2 launches per block per step, so b blocks correspond to N=2b:

| Blocks | Launches | Added µs | Projected µs/step | Time vs current | With Graph |
|---|---|---|---|---|---|
| 1 | 2 | +3 | 455 | 1.01x | 1.00x |
| 2 | 4 | +13 | 465 | 1.03x | 1.01x |
| 4 | 8 | +23 | 475 | 1.05x | 1.01x |
| 8 | 16 | +212 | 664 | **1.47x** | **1.03x** |
| 16 | 32 | +467 | 919 | **2.03x** | **1.07x** |
| 32 | 64 | +1226 | 1678 | 3.71x | 1.15x |

Both columns are relative to the current 452 µs; the [B] data at N ≥ 32 has about ±20% spread, so the last row is order-of-magnitude only.

**Implications**:

- A naive multi-block implementation stops paying off at **8 blocks or more**; **with Graph, the cost stays within 15% up to 32 blocks**.
- The conclusion flips once the grid is enlarged: at 400³ a step is ~29 ms, and 64 launches at 0.4 ms is only 1.4%. **So the essence of the cost is "cells per block", not "block count"** — the criterion should be `cells_per_block ≥ 5e4 ~ 1e5`, not a cap on the block count.

### 3.4 Resulting design constraint

> **Recommended lower bound per block ≈ 5×10⁴ cells (about 37³).**
> The lower bound comes from three things at once: the occupancy floor of 19200 cells, overlap redundancy starting to exceed 12% below 50³, and the launch-cost share requiring a few hundred µs of work per block.

## 4. Implementation plan

### M0 Safety guard (partially landed)

| Item | Content | File | Status |
|---|---|---|---|
| M0-1 | `BlockLatticeManager::CuDev*` **errors out** when `BlockLats.size() > 1`, no longer silently computing block 0 only | `block_lattice.hh:991-1018` | **Not done (open correctness item)** |
| M0-2 | both cell-dynamics kernels get an `if (idx >= n)` bounds check | `block_lattice.h:292`, `:306` | **Landed** |
| M0-3 | checksum covers all blocks (`getBlockLat(0)` → loop over all) | `benchmarks/cavity3d_cu/cavity3d.cu`, `field_checksum.h` | **Not done (open correctness item)** |

Rollback: M0-1/M0-3 only add a guard and a loop; no behaviour change (the single-block path checksum must stay bit-identical).

### M1 Single-card same-resolution multi-block (deferred)

**§8 takes priority**: do not start M1 until the two measurements in §8.7 and the active-index compaction in §8.2 are complete and the Go condition in §7 holds; if §8.7 decides "stop after §8.2", M1 is not done.

| Step | Content | File | Estimated size |
|---|---|---|---|
| S1 | restore the commented-out multi-block loop (`CuDevStream` / `CuDevApplyCellDynamics`, three entry points) | `block_lattice.hh` | ~15 lines |
| S2 | merge the stream launches: one kernel, grid = Nblock × q, replacing one `<<<q,1>>>` per block | `block_lattice.h` / `.hh` | ~20 lines |
| S3 | device-side halo kernel: H2D the `CommDirection` / send-recv index tables, write gather/scatter kernels; **same resolution only at first (plain copy, no interpolation)** | new `cuda_comm.h` | ~150 lines |
| S4 | one large VRAM allocation per block, replacing 5–6 small `cudaMalloc` | `block_lattice.h` | ~30 lines |
| S5 | benchmark: parameterise `block_num` + MLUPS sweep | `cavity3d_cu` | ~30 lines |

**Verification**: see §5. Rollback: bit-identical to the current state when `block_num = 1`.

### M2 Cost removal (about 1–2 days, depends on the M1 measurements)

- S6: CUDA Graph, one `cudaGraphLaunch` per step (reusable while the block count is unchanged; a change in block count needs a rebuild);
- S7: nsys to locate the real source of the sharp degradation in [B] (CPU enqueue vs GPU front end), and decide from that whether S8 is worth doing;
- S8 (optional): merge collision and stream into the same kernel, launch count N+1 → 1.

**Gate**: M2 stops at "≤15% cost with ≤32 blocks"; do not chase more.

### M3 Cross-resolution / AMR (phase 2, separate item, not implemented in this plan)

Three hard problems; do only a feasibility assessment first:

1. **Interface interpolation**: halo between different Δx needs space-time interpolation (`popIntp` is host-only today);
2. **Local time stepping**: the fine level takes 2 steps per coarse Δt while the coarse level takes 1 — this conflicts with the existing main-loop structure of "two kernels over the whole field per step", and the main loop has to be rewritten;
3. **Rotation streaming breaks**: `CyclicArray::rotate` relies on a uniform in-block shift, and at a coarse/fine interface it must degrade to a real gather/scatter.

### M4 Multi-GPU (long term)

One rank per card, block-to-device mapping plus reuse of the existing MPI communication. Depends on the M3 conclusion.

## 5. Verification plan

**Correctness (three levels, the same criteria as FP16_PLAN)**

1. **Regression**: with `block_num = 1`, the GPU output checksum is **bit-identical** to before M0;
2. **Cross-implementation**: GPU multi-block vs CPU multi-block (same `block_num`), 1000 steps, 100³, float — target relative difference in the sum ~1e-7 (single-block GPU vs CPU is known to be 2.3e-7);
3. **Physical invariant (strongest criterion)**: run the same case with 1 / 2 / 4 / 8 blocks and compare the velocity profile and total mass after convergence — **splitting must not change the physical result**; any quantity that drifts with the block count means the halo is wrong.

**Performance**: 4. MLUPS curves for `block_num ∈ {1,2,4,8,16,32,64}` × grid `{100³, 200³}`, plotted as "cost vs cells per block" rather than vs block count; 5. nsys capture of launch gaps, SM occupancy and halo-kernel share; 6. before/after Graph comparison (the comparison method already exists in §3.1 [B]).

**Stability**: 8 blocks, a 3000-step long run with no NaN and mass conserved; follow the way the FP16 long-run failure boundary is recorded.

## 6. Risks and rollback

| Risk | Description | Handling |
|---|---|---|
| **Silent error** (highest) | `block_num>1` computes one block only while the program prints completion as usual — an incident of this class has already happened once | M0-1 unconditionally first |
| Communication through the host | Once halo is done through D2H/H2D as a shortcut, it costs ms per step, the gain is gone and it is hard to notice | S3 is a hard gate; no performance evaluation until it is finished |
| CUDA Graph vs dynamism | A change in block count/grid needs a graph rebuild; AMR dynamic refinement triggers it often | M2 is limited to "static block structure"; dynamic refinement is left to M3 |
| **Over-engineering** | If M1 measures a loss at ≥8 blocks that Graph cannot rescue, **stop at M0 + active-index compaction**, keep multi-block only for AMR, the "must have" scenario, and never as the default path | §7 gate |

## 7. Decision gate (Go / No-Go)

**Current state: no recorded verdict.** None of the conditions below has been met or rejected by any measurement or review; M1 does not start before the two measurements in §8.7 exist.

**Go** — any one of: **AMR on GPU** is needed (FreeLB's core capability, currently not supported by the GPU path at all); **multi-GPU** is needed; in a complex-geometry case the void fraction is high enough to justify a geometry-fitted block decomposition (**note: try active-index compaction first, that path keeps the launch count unchanged**).

**No-Go**: only "to speed up a uniform grid" — do not do it. FP16 has already measured 1.86x (2211 → 4113 MLUPS), an order of magnitude more reliable than multi-block; do that first.

**M2 pass line**: with ≤32 blocks and ≥5×10⁴ cells per block, the multi-block cost relative to a single block is ≤15%. If that is not met, position multi-block as an "AMR-only path" and state in the doc that it is a net loss on a uniform grid.

## 8. The void cost and the VRAM trade-off under complex geometry (decides whether M1 starts)

This addresses the specific dilemma of "one card, no multi-GPU in the near term, but the bounding box of a complex geometry is far larger than the solid". The conclusion is
**do the active-index compaction of §8.2 first, and defer multi-block (M1) — one account has not been settled yet, and it is not a VRAM account.**

### 8.1 Void cells are currently paid for in full, and the currency is bandwidth

The default path (`--reg`, i.e. `cudev::Cell<..., cudev::RegPop>`) unconditionally moves the q distribution functions into registers in the RegPop constructor in `src/data_struct/cell.h` (`PopCache::load`), and the `flag` dispatch happens after that; `cell.flush()` likewise writes back unconditionally. So a void cell also generates 19 loads + 19 stores = 152 B/cell of traffic. Compare `--base` (`DirectPop`): `operator[]` returns a global-memory reference only, evaluated lazily, and a void cell triggers no access.

**Implication**: void cells under `--base` generate no traffic; under `--reg` (the default) they generate the full 152 B/cell; the 100³ cavity is 100% filled, so this cost is **completely invisible** in the existing benchmark. Counting useful cells, **effective MLUPS = measured MLUPS × φ** (φ = volume ratio of solid to bounding box): at φ = 0.3, 2211 MLUPS is really only about 663 effective MLUPS — an invisible loss of the 3x order, an order of magnitude larger than the multi-block throughput loss (−32% at 8 blocks).

### 8.2 The right order: take the step with a pure gain first

**Active-index compaction + deferred fetch** (about 1 day, orthogonal to multi-block)

- the kernel no longer starts threads from `0..N−1` but from the active-cell index table. The CPU side already has `GenericvectorManager` and the `ApplyCellDynamics(blockids)` overload, but the `BulkTaskIds` task-table build is commented out in the benchmark, and the GPU path does not consume a task table yet — the interface layer exists, the task table is not enabled;
- **rotation streaming makes this unusually clean**: in pull mode each thread accesses only the 19 slots of its own cell (`data_d[i]`), it does not gather neighbours. So a compacted index **needs no neighbour table, no change to `CyclicArray`, no halo**;
- gain: traffic per step ∝ N_active; cost: still 1 launch, and memory locality depends on the contiguous run length of active cells along the fastest axis; **it does not save VRAM**. This step wins in both outcomes, "do multi-block" and "do not", so it cannot lose.

### 8.3 The VRAM account — and an explicit threshold

Each cell is about 93 B (same arithmetic as §2.2), 6 GB ⇒ about 6.5e7 cells ≈ 400³. Criterion: when `N_bbox × 93 B ≤ 0.7 × available VRAM`, do only §8.2 and do not touch multi-block (keeping the §7 No-Go); only when VRAM is a real bottleneck is tiled allocation even discussed. The reason is plain: **wasting VRAM you already own costs nothing.** Tiling is only worth it when "it does not fit" or "the resolution must be lowered" — and the accuracy lost by lowering the resolution is far larger than the 20–40% of halo.

### 8.4 If VRAM really is the wall: the tile size tension is real

Tile edge m, solid volume V, surface area A (in cells): VRAM needed ≈ V·(1 + 6/m) + A·m, per-step time penalty ≈ 6/m (halo redundant computation) + launch cost, VRAM optimum m* = √(6V/A) (a sphere R=100 → m ≈ 14; a thin plate 8 cells thick → m ≈ 5).

**The VRAM optimum lands on small tiles (halo penalty 40%–100%), while the time penalty demands large tiles.** This tension is real and cannot be removed. But the launch term in the tension is not essential: do not write "one launch per block", write **one kernel walking the tile table** (grid-stride over tiles), with the launch count constant at O(1)/step; then the 6.6 µs/launch measured in §3.1 is no longer a function of the tile count. The remaining 6/m cannot be dodged, so the rule for choosing m becomes: **take the largest m the VRAM budget allows**, not the VRAM-optimal m*.

### 8.5 The third way was rejected by measurement: sparse commit of virtual memory

The idea: keep the dense index unchanged (constant offsets, rotate streaming, coalesced access all unchanged) and commit physical pages only for non-empty VA ranges — in theory saving VRAM at zero performance cost. It was evaluated with a temporary probe using `cuMemAddressReserve` + `cuMemCreate` + `cuMemMap` + `cuMemSetAccess` (not committed; deleted with that experiment).

Measured (this card): granularity 2 MB; FP32 POP gives 524288 **contiguous** cells per chunk (cube edge 80.6, 3.27 xy planes at nx=400), doubled for FP16. In linear (i,j,k) order, as long as the geometry spreads along the slowest axis, not a single chunk can be saved; finer granularity requires Morton order, which breaks constant offsets → rotation streaming fails → fall back to gather streaming. **The path does not hold without changing the streaming mechanism.** Side effect: the old VMM path rounded each of the 19 directions up to 2 MB, which is the direct reason it used more VRAM on small cases and one of the reasons `CyclicArray` is the default.

### 8.6 The cheapest VRAM lever is already in hand

FP16 storage has landed (see `docs/FP16_PLAN.md`, 4113 MLUPS): 93 B → about 55 B per cell, **scaling the computable size directly by 1.69x**, with no tile/halo/launch complexity. Use this lever to the full before touching multi-block.

### 8.7 The two measurements that decide this (both quick)

1. **Decide whether 8.1 holds**: with the existing code unchanged, set the flag of a subset of cells in cavity3d to void (keeping N constant) and see whether MLUPS changes. **If MLUPS barely moves → it is confirmed that void cells consume bandwidth in full**, and §8.2 is then the single biggest gain available today.
2. **Get φ and A/V for a real geometry**: for an actual complex geometry, count `getTotalCellNum()` vs the active count by flag. φ > 0.5 and it fits within 400³ ⇒ stop after §8.2, never do multi-block.

Until these two numbers exist, every performance discussion about multi-block is an estimate.

## Appendix: probe

`benchmarks/multiblock_probe/gpu_multiblock_probe.cu` — empty-kernel launch cost, the cost of splitting a fixed workload, and the CUDA Graph comparison, all in one. Every number in §3.1 comes from it and can be reproduced directly.
