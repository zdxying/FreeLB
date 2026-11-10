# GPU POP storage and policy (RegPop/DirectPop, CSE, FP16)

Status: implemented (FP16 off by default; RegPop is the default of `benchmarks/cavity3d_cu`) | Updated: 2026-11-09 | Target: the four combinations FP32/FP16 × RegPop/DirectPop for D3Q19 are available and measured in `benchmarks/cavity3d_cu`; FP16 changes POP storage only, and FP32 stays bit-identical while it is off by default.
Exit condition: when FP16S (half2 packing) lands and is re-measured, FP16 is turned on by default, the POP storage format changes, or `cavity3d_cu` unifies `--reg`/`--base` into a single policy or adds a third POP policy, the storage and policy decision moves to the CUDA mirror in `docs/ARCHITECTURE.md` §6 and the numbers here are superseded by a new table.

## 0. Conclusions

- FP16 storage is usable and suits scenarios insensitive to extreme-value precision; quantitative fine studies still use FP32 (the default), and the FP32 path stays bit-identical while FP16 is off.
  The performance gain depends on RegPop: RegPop + FP16 is the fastest combination, 87% faster than FP32 + RegPop and 2.9x faster than FP16 + DirectPop; under DirectPop FP16 gains almost nothing.
- The policy is the 4th template parameter of `cudev::Cell`, and one CSE artifact covers both policies; the RegPop CSE specialization is selected (FP64 19 → 0).
- Resolving the `p[Q]` addresses once is a performance precondition; the per-direction recompute fallback is measured in §7.

## 1. Storage type and compute type are separate

`PopStorage<T>` (`src/utils/alias.h`) defaults to `= T`; when `FREELB_POP_FP16` is defined and `T = float`, it is `= __half`. **Only the storage element type of the POP field changes**; RHO/VELOCITY/flag and every macroscopic field stay float.

- The FP32 path runs 2211 MLUPS ≈ 336 GB/s effective bandwidth, already at the nominal peak of the RTX 3060 Laptop (336 GB/s); the minimum traffic per cell per step is 152 B (19 reads + 19 writes × 4 B). FP16 halves the traffic at 2 B/pop while the arithmetic stays FP32.
- Precision cost: a 10-bit mantissa, ~5e-4 relative rounding per write-back. FluidX3D has already validated feasibility in its FP16S/FP16C modes; known failure modes: very long runs, high Ma, multiphase with a large density ratio. On the same card FluidX3D FP16S is 4012 MLUPS, this implementation FP16C about 4113 MLUPS — comparable.
- The shared-memory route has been rejected by an ablation experiment (gather at the same index, zero cross-cell reuse); FP16 is the only next step whose gain is of the right order.
- Only POP is touched: RHO/VELOCITY/flag stay float (under 5% of the traffic, no gain). `FREELB_POP_FP16` is not defined by default; `benchmarks/cavity3d_cu/Makefile` has `FP16 ?= 0`, and `make FP16=1` when needed.

### PopRef

`cudev::Cell::operator[]` (`src/data_struct/cell.h`) returns `decltype(auto)`: `T&` when the storage type equals the compute type; otherwise the proxy `PopRef<PopStorage<T>, T>` (it wraps a pointer only, converts implicitly to float and assigns from float). Dynamics code needs no change.

The RegPop path needs almost no change: `PopCache<T, Q, RegPop>::v` is `T` (float), and the storage is read and written through the `__half` implicit conversion in `p[d][0]`. That conversion goes through `__half` rather than `PopRef`; both land on the same SASS instruction, so unifying them is for readability only and not worth the change.

`PopRef` must be transparent to assignment, or the wall populations are silently lost; the failure mode and the evidence are in §7.

## 2. The policy is a template parameter of Cell

### D1 — the policy is a type tag appended to the end of the template parameter list, with the default on the forward declaration

`src/data_struct/pop_policy.h` defines two empty tags, each providing only `block_size` (`DirectPop` 32, `RegPop` 128), and gives the single default:

```cpp
template <typename T, typename LatSet, typename TypePack,
          typename POPPOLICY = DirectPop> class Cell;   // the only place with the default
```

A default argument may appear only on the first declaration. The declaration lives in `pop_policy.h` because code such as `CellType` in `data_struct/cuda_block_lattice.h` spells the 3-argument form directly, and `data_struct` must not depend on `lbm` in reverse; `cell.h` satisfies it through `pop_policy.h`, independent of include order.

The CSE specialization pattern must put `POPPOLICY` **last**, aligned with the alias arity; putting it in the middle misaligns it against the alias, and every specialization silently fails to match (it compiles, the numbers are correct, it just falls back to the loop version).

### CSE exact type matching

The partial specializations emitted by CSE match on the **exact cell type**:

```cpp
template <typename T, typename TypePack, typename POPPOLICY, bool WriteToField>
struct rhoUImpl<CELL<T, D3Q19<T>, TypePack, POPPOLICY>, WriteToField> { ... };
```

Partial-specialization matching performs no base-class conversion: a derived cell class does not match and silently falls back to the primary template with `for (i < LatSet::q)` — it compiles, the numbers are correct, it is just slow. `cudev::RegCell` is only an alias of `Cell<T, LatSet, TypePack, RegPop>` (`cell.h`), not a derived class.

`generated/lbm/*.ur.h` is one artifact shared by all targets; hard-coding the policy into the artifact would deprive a target that uses the other cell of its specializations. The policy is therefore a template parameter of `Cell`, forwarded by csegen: one emission covers every policy, the size does not grow, and a new policy needs no csegen change.

| Path | cell type | CSE specialization |
| --- | --- | --- |
| CPU `benchmarks/cavity3d` | host `Cell` | in effect (`-D_UNROLLFOR` in `Makefile`) |
| `examples/cavity2d_cu` | `cudev::Cell` | in effect (`using CELL` in `cavity2d.cu`) |
| `cavity3d_cu --base` | `cudev::Cell` | in effect (`using CELL` in `cavity3d.cu`) |
| `cavity3d_cu --reg` (default) | `cudev::RegCell` = `Cell<...,RegPop>` | in effect (`g_UseRegCell` / `using RegCELL` in `cavity3d.cu`) |

## 3. Design invariants D2–D5

### D2 — the register cache lives in a separate storage class; private inheritance eats EBO

The primary `PopCache<T, Q, Policy>` template is an empty class; `PopCache<T, Q, RegPop>` stores `PopStorage<T>* p[Q]` (q element addresses) and `T v[Q]` (q values). C++17 has no `[[no_unique_address]]` (`-std=c++17` in `make.mk`), and an empty member still occupies bytes, so a private base class is required; otherwise DirectPop would waste 19 registers and destroy the baseline path (A1 exists for exactly this).

### D3 — `getNeighbor` always returns `DirectPop`

A neighbour is returned by value and the caller usually reads one direction only, so giving the neighbour a 19-register cache is pure waste; a form like `BounceBackMovingWall` that fetches several neighbours inside one kernel multiplies that further. The alias `NeighborCell = Cell<T, LatSet, TypePack, DirectPop>` inside `Cell` states this explicitly and must not be "changed to inherit the policy" in passing.

### D4 — resolve the addresses once at construction, keep the `p[Q]` pointer array

The `RegPop` constructor calls `BlockLatticeBase::getPopArray` once to fill `p`, and `flush()` writes `v` back; under DirectPop both are no-ops, so one kernel serves both policies. Addresses must **not** be recomputed per direction inside the `load`/`store` loops; the measured regression on FP16 and its cause are in §7.

### D5 — `--reg` / `--base` is a runtime switch, the policy is a compile-time parameter

`g_UseRegCell` in `cavity3d.cu` selects between the two instantiated kernels; `block_lattice.h` keeps only two kernel templates (with/without flag), and the launch geometry is determined by `CELLTYPE::block_size`. The argument parsing of the `.cu` file needs no rewrite, and the ability to run the controlled comparison is preserved.

Cost: the cell type cannot be inferred back from `CELLDYNAMICS` — `tmp.h` has two coexisting `TaskSelector` definitions, one with a fixed arity and one variadic, so inferring the cell from the task type is unreliable. Call sites therefore pass `CELLTYPE` explicitly, kernels are backed by `static_assert(std::is_same_v<typename CELLDYNAMICS::CellType, CELLTYPE>)`, and both `TaskSelector` definitions provide `using CellType = CELL`.

## 4. CPU path and fallback

Every entity carrying storage/instructions (`PopCache`, `PopRef`, `cudev::Cell`, `RegCell`) sits inside the `#ifdef __CUDACC__` block of `cell.h`; the host `Cell` and `GenericCell` are unaffected. `PopRef` is instantiated only when the storage type differs from the compute type, and under CPU `PopStorage<T> == T`. `pop_policy.h` holds only empty tags and a forward declaration: visible on the host side but generating no code. The generated host alias `using CELL = Cell<T, LatSet, TypePack>` ignores the 4th parameter, `POPPOLICY` degrades to a free unused template parameter (legal), and it still matches the CSE specializations for the host `Cell`.

Evidence: preprocessing a CPU translation unit of `examples/cavity2d` (g++ + `-D_UNROLLFOR`) shows neither `PopCache`, `PopRef`, `RegCell`, `cudev::Cell` nor `CuDevApplyCellDynamics`; all that enters the CPU compilation is the empty tags of `pop_policy.h`, the `Cell` forward declaration and the two `block_size` constants.

The fallback rule is in `make.mk` and `lbm.h`: `_FLB_UR_CSE` is defined only when `-D_UNROLLFOR` is set and csegen exists; otherwise `lbm.h` compiles the loop version only and includes no `.ur.h`. `src/ur_fallback/lbm/*.ur.h` is a hand-written 3-parameter reference that is legal thanks to the default `DirectPop`, but it is not on the include path and serves only `tools/cse` validation.

## 5. Verification results

D3Q19 / `CyclicArray` / RTX 3060 Laptop (sm_86) / 100³; A1–A5 are FP32, 1000 steps; the FP16 rows reference FP32 with the same configuration, steps in the project name. Summary of measurements from 2026-09-26 onward.

| Item | Result |
| --- | --- |
| Build | all targets build; the FP16 × cyclic compile matrix passes |
| Off by default | FP32 checksum `1061231.05 / 4749162.144 / 33.7224`, bit-identical to the FP32 baseline |
| A1 DirectPop (BGK / RhoU) | 1712 instructions / 216 FLOP / 0 FP64 / REG:44; 680 / 76 / 0 / REG:28, both exactly the same as the baseline |
| A2 RegPop FLOP | 314 → 226 (-28%) |
| A2 RegPop FP64 | 19 → 0, the specialization is selected |
| A3 RegPop registers / code generation | REG:96 / LOCAL:0, no spill; the mangled name confirms POP is `CyclicArray<__half>` under FP16 |
| A4 numerics | `--reg` and `--base` checksums bit-identical: `sum 1061231.284 / sumsq 4749162.209 / maxabs 33.72239685 / NaN 0` |
| A5 `moment.ur.h` size | 192555 B → 195378 B (+1.5%); a twin-artifact approach would add +200KB |
| Container | the `CyclicArray` device mirror is bit-identical to the exact-modulo reference |
| Kernel boundary | the cell-dynamics kernel takes an `n` parameter and returns early on `idx >= n`; that guard does not change the result, FP32/FP16 checksums are bit-identical |
| FP32 GPU vs CPU | sum relative difference 2.3e-7; field max velocity difference 1e-6, RMS 9.4e-8, max u_x 0.2 (lid velocity) |
| FP16 1000 steps | sum relative difference 8.8e-5, sumsq +1.1e-2, maxabs 33.0 vs 33.7224 (boundary-layer extreme, trajectory divergence) |
| FP16 3000 steps | no NaN/blow-up, mass conservation 1.4e-4, sumsq -5.3%, maxabs 95.5 vs 100.4 |
| FP16 10k steps | no NaN/blow-up, mass conservation 1.2e-3, sumsq -17%, maxabs -23%; MLUPS FP32 1998 / FP16 3421 |
| Test gate | `tests/run_tests.sh` does not build FP16; FP16 needs a manual `make FP16=1` run of `tests/cuda/checkpoint` |

The fingerprints come from `src/utils/field_checksum.h` (including the NaN count and the first flat index); the SASS statistics use `tools/sass_census.py`.

**Conclusion: FP16 storage is usable and suits scenarios insensitive to extreme-value precision; quantitative fine studies still use FP32 (the default).**

## 6. Performance matrix

D3Q19 / `CyclicArray` / RTX 3060 (sm_86) / 100³ lid-driven cavity / 1000 steps, 3 runs per table cell. Unit MLUPS.

| | RegPop (register-resident) | DirectPop (streaming) |
| --- | --- | --- |
| **FP32** | 2179 / 2179 / 2175 | 1396 / 1338 / 1364 |
| **FP16** | 4162 / 4113 / 3960 | 1444 / 1413 / 1348 |

- Three-run means: FP32 DirectPop 1366 / RegPop 2178; FP16 DirectPop 1402 / RegPop 4078 (range 3960–4162, about 1.84x FP32 RegPop).
- **RegPop + FP16 is the fastest combination**, **+87%** faster than FP32 + RegPop and **2.9x** faster than FP16 + DirectPop. Memory access is still the main bottleneck; "cache the pop in registers" and "halve the pop" stack and are not mutually exclusive.
- On FP32 the arithmetic gain of CSE does not show: RegPop is 1.6x faster than DirectPop, but the whole gap comes from memory access (CSE removes only 28% of the floating-point work, and under a memory bottleneck that does not change throughput). **Under DirectPop FP16 gains almost nothing** (1366 → 1402, +2.6%): it re-resolves the container address on every element access, so the data traffic is not saved; FP16's gain therefore depends on RegPop, and the default of `benchmarks/cavity3d_cu` (RegPop) is the right default.

## 7. Known limitations and open items

### FP16

1. **Long-run precision**: FP16 rounding amplifies with steps/gradients. At 10k steps there is still no NaN/blow-up, but sumsq -17% / maxabs -23% (§5); quantitative fine studies use FP32 (the default), and beyond 10k steps is not characterised. This is the failure boundary of the optional speedup.
2. **DirectPop + FP16 is CSE-path only**: `latset::c<LatSet>(i) * cell[i]` in the loop version `src/lbm/moment.h` has no `Vector<int,D> * PopRef` overload, so without `-D_UNROLLFOR` FP16 + DirectPop fails to compile; the CSE-generated version (`generated/lbm/moment.ur.h`) first lifts `cell[i]` into `const T` and is unaffected. RegPop compiles on both paths.
3. **PopRef must be transparent to assignment**: `collision::BounceBack` (`src/lbm/collision.h`; the unrolled version `src/lbm/collision_unroll.h`) is written as `cell[i] = cell[iopp]`; if `PopRef` has only `operator=(CT)`, overload resolution picks the implicit copy assignment, the write lands in a discarded temporary, and the wall reflection never reaches device memory. Without a write-through `operator=(const PopRef&)`, FP16 turns 95% of the cells into NaN within 100 steps; after the fix the same case goes from NaN count 19,151,655 → 0, and the checksum agrees with the RegPop path within FP16 quantization noise. FP32 does not instantiate `PopRef` and is unaffected.
4. **MPI × FP16 is unverified**: refblock is compile-checked only; the interaction of multi-block GPU halo communication × FP16 is unknown. Exit condition: a multi-block GPU MPI runtime case builds with FP16 storage and passes, with mass conservation recorded.
5. **Per-direction address recomputation**: dropping `RegPop`'s `p[Q]` and recomputing the address per direction in load/flush is on par for FP32 and 1.75x slower for FP16 (2400 vs 4000), 64-bit address loads 90 → 217; once the payload is halved, address resolution becomes the dominant term. The invariant is §3 D4.

### CSE coverage and code generation

- **Asymmetric coverage**: `src/lbm/collision_unroll.h` and `collisionMRT_unroll.h` are hand-written resident headers, not in `UR_CSE_BASES` of `make.mk`, and still use the 3-parameter `cudev::Cell<...>` → `DirectPop`. `BounceBack<RegCell>` / `BounceBackMovingWall<RegCell>` therefore get no unrolled specialization and fall back to the loop version; this is a pre-existing coverage asymmetry, unrelated to the policy parameter, and its effect is confined to wall cells.
- Bringing collision/MRT into CSE has been evaluated by measurement and is deferred: the engine front end does not support how their primary templates are written (C array declarations, `unsigned int`) and does not unroll stride-2 loops; MRT has no real user, and BounceBack is a pure swap with 0 FLOP. The engine gap list is in `tools/cse/Makefile`.
- `collision::BGK` itself is still the loop version and is not marked `@cse`; `latset::w<LatSet>(k)` is still not folded (several places in equilibrium/force), and the device side reads `__constant__` through `Fraction::operator()`, whose body contains `std::cout` + `exit(-1)` (`utils/util.h`), so it currently compiles only because nvcc constant-folds the literal subscripts.
- `formatConst` of csegen still emits bare double literals with `setprecision(17)`; nvcc folds them in the measurements, and the dynamics kernels have 0 FP64, but there is no guarantee.
- **Comments**: the `PopCache` / `cudev::RegCell` comments in `cell.h` and `getPopArray` in `cuda_block_lattice.h` all describe "resolve once at construction", consistent with the code.

## 8. Phase 2

- FP16S: half2 packing + a 2 cell/thread restructure, closing the ~12% gap from FP16C to FP16S (on the same card FP16C about 4113 MLUPS; FluidX3D FP16S 4012 MLUPS).
- Explicit `__float2half2` vectorisation, `cp.async`/L2 policy evaluation.
