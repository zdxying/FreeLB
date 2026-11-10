# Checkpoint/restart: implementation status and the P4 plan

Status: P1–P3 implemented and merged into `dev` | Updated: 2026-11-09 | Target: P4 (CA / runtime AMR)

> Exit condition: once the P4 acceptance runs for `cazsblock2d` / `cazsAref2d` pass and land in §3, this
> document is compressed into a present-state section of `docs/ARCHITECTURE.md`; until then it is the
> decision record for checkpoint changes.

## 0. Conclusion

- `src/io/checkpoint.{h,hh}` (namespace `flbckpt`) is implemented, works on CPU and GPU, and supports restart
  across a **change of rank count** (preconditions in §1.3). Phases: P1 CPU serialization and restart → P2 GPU
  single-block restart (including FP16) → P3 restart across a rank count → P4 CA / runtime AMR (§5, not done).
- Wired examples: `cavity2d`, `cavity3d` (CPU, writes geometry identity), `cavity2d_cu` (GPU single block); all
  three `.ini` files default to off (`Restart = 0`, `CheckpointStep = 0`).
- Since format v2 a snapshot can carry the **geometry identity** (grid / origin / voxelSize / per-cell entity /
  boundary patch / digest); a restart verifies the digest and raises a hard error on mismatch; callers that do
  not pass a `VoxelGeometry` (2D) still round-trip, but there is nothing to verify.
- Verification in §3. Not supported: CA / runtime AMR restart (P4), cross-architecture / byte order, cross
  FP16/FP32 builds (hard error only, no conversion).

## 1. Implemented capabilities

### 1.1 Serialization scope

| State | Disposition | Notes |
|---|---|---|
| Global step counter and per-block metadata | stored | restored by `Counter::set()`; `blockid / level / overlap / mesh / n` verified per block on load |
| All fields of every `BlockLatticeManager` / `BlockFieldManager` | stored | POP, RHO, VELOCITY, free-surface STATE/MASS/VOLUMEFRAC/MASSEX/PREVIOUS_VELOCITY, CA STATE/FS/EXCESSC, `Data<>` scalars |
| Geometry identity (v2, optional) | stored | when a `VoxelGeometry` is passed, writes grid / origin / voxelSize / per-cell entity / boundary patch / digest; `readWithGeometry` verifies the digest and the grid, mismatch is a hard error |
| BC / task definitions | not stored | replayed by the example's own setup with the same parameters |
| `UOld/RhoOld` (tolerance history) | not stored | the example calls `EnableToleranceU()` after loading to reseed, so the first `res` after a restart is already the correct value |
| CA `Interface` index table, runtime AMR block table | not stored | non-field persistent state; the geometry identity describes one static grid only, P4 handles this |

### 1.2 Wiring and usage

Three optional keys in `[Simulation_Settings]` (absent means off, old `.ini` files are unaffected):

```ini
Restart = 0            ; 0 = new run, 1 = restart from the latest snapshot
CheckpointStep = 0     ; > 0 writes every N steps (0 = write only on exit)
RestartDir = ./ckpt    ; snapshot directory
```

The API is in `src/io/checkpoint.h`: `write` / `read` / `checkCompatible` / `readHeader` /
`latestStep` (`:136-168`) and the geometry-aware `writeWithGeometry` / `readWithGeometry` /
`checkCompatibleWithGeometry` (`:152-162`). `read` with `step < 0` reads `latest` and returns the step
it restored to; the variadic arguments may mix any number of block-lattice managers (including CA-derived
classes) and block field managers. `read` / `checkCompatible` reject files that carry a geometry identity
(there is nothing to verify them against); those must use `*WithGeometry`.

### 1.3 Restart semantics

- **Quiescent-point rule**: `write` may only be called once communication has finished and before the next
  `Stream()`; the examples place it in the output block of the main loop and never inside an AMR
  refine/coarsen step.
- **Logical view**: POP is written through the logical addresses of `getdataPtr()` and the file holds no
  rotation state; on load the data is `normalize()`d first and only then laid down in the logical view,
  otherwise it would be permuted by the stale rotation state.
- **GPU mirror**: `copyToHost()` before writing, `copyToDevice()` after reading, plus
  `CyclicArray::refreshDeviceObject()` to re-upload the device object (`copyToDevice()` does not push
  `start[]`).
- **Bitwise identity**: with the same parameters and the same block set, a restarted trajectory is
  bit-identical to an uninterrupted run, not merely "consistent within error".
- **Changing the rank count**: the read path picks its own blocks out of all rank files by global
  `blockid`, so restarting across a different rank count or a different block-to-process assignment
  works; the write path still produces one file per rank. The block set must be independent of the rank
  count: `AdaptiveOptimization(mpi().getSize())` repartitions the domain by process count and is
  explicitly rejected when the rank count changes; use `GeoHelper.CreateBlocks(bx, by[, bz])` +
  `LoadBalancing(np)` instead.
- **Hard error rather than silence**: on a geometry digest, block table or POP storage precision mismatch
  the first differing block is printed and the run `exit(1)`s; when the block counts differ it reports
  "a block cannot be split" (the decomposition can only be made coarser, never finer).

## 2. File format

Directory: `<RestartDir>/latest` (the step number last written by rank0 as text; its presence means that
step is complete) and `<RestartDir>/step_<s>/rank<r>.ck` (one binary per rank). The version number, record
layout and byte layout are in `src/io/checkpoint.h` as `kFormatVersion` / `BlockRecord` / `FieldRecord` /
`GeometryRecord`; files use the host byte order (no cross-architecture / cross-byte-order support).

- v2 adds an **optional geometry identity block** after the fixed header, written only by `*WithGeometry`;
  when `present = 0` the file still round-trips. `entity` is an `nx*ny*nz` per-cell record index; patches
  use record indices, not coordinates.
- Matching is by `field_index + array_dim + elem_tag + count` (`name` is used for logging only);
  `Data<>` scalars have no count prefix, their size is determined by `elem_tag/elem_dim`
  (`array_dim == 0` in the metadata).
- **Atomicity**: each rank writes `rank<r>.ck.tmp` → barrier → `rename()` → barrier → rank0 writes
  `latest` as `tmp + rename`; a crash never leaves a half snapshot that is taken as complete. Old
  snapshots are not cleaned up automatically.

## 3. Verification status

| Test | Coverage | Result |
|---|---|---|
| `tests/feature/checkpoint` | CPU bitwise restart, CyclicArray rotation-state round trip, header verification, negative block-table checks | passes, `PopChecksumHost` sum/sumsq/maxabs bit-identical |
| `tests/feature/checkpoint_geometry` | 3D geometry identity header round trip, digest rejects a different record, a snapshot without a record is rejected by `checkCompatible()` | checked in, run by `tests/run_tests.sh` (not re-run this round) |
| `tests/cuda/checkpoint` | GPU: FP32 and `make FP16=1`; cross-storage rejection | passes; an FP16 snapshot read by FP32 reports an element-type mismatch and exits |
| `tests/block/checkpoint_decomp` | 2 ranks write → 1 rank reads and restarts (`MPI_Comm_spawn`) | passes, per-block interior-field checksums show 0 difference |
| `tests/block/checkpoint_geometry_mpi` | 2 ranks each write an identity block; on read the identity of each rank file is verified to agree | checked in, run by `tests/run_tests.sh` (not re-run this round) |
| `tests/regression/cavity2d`, `examples/cavity2d` / `cavity3d` / `cavity2d_cu` | numeric path untouched; end-to-end continuous run vs interrupted-and-restarted | passes, res identical |

Build rules before verifying MPI paths are in `docs/ARCHITECTURE.md` §5.

## 4. Known limitations

1. **GPU supports a single block only** (`src/data_struct/block_lattice.hh:991-1018`, `CuDev*` iterates
   `BlockLats[0]` only); P2's scope is exactly "single-block GPU restart", consistent with that limitation.
2. **Changing the rank count requires a block set independent of the rank count** (§1.3); the default
   `AdaptiveOptimization` path is rejected.
3. **partial refreshes only the direction set of `getCommPopDir()`** (pre-existing defect, unrelated to
   checkpoint, not fixed; it is not "only missing corner cells", the direction set itself is not closed).
   Mechanism, unfixed items and the applicability criterion are in `docs/INTP_COMM_OVERFLOW_BUG.md`
   §18 / §18.9 (index out of bounds on AMR in §18.8.3).
4. **CA / runtime AMR examples are not wired**, see P4; **old snapshots are not cleaned up** (§2).

## 5. P4 plan (CA / runtime AMR)

### 5.1 CA `Interface` index table persistence

Facts: each block's `Interface` is seeded once by `Setup()` (`src/ca/zhu_stefanescu2d.hh`) and rebuilt at
runtime by `UpdateInterface()`; the boundary manager holds a **pointer to that vector**
(`ZhuStefanescu2DManager::Interfaces`, `src/ca/zhu_stefanescu2d.h`).

- Add a read/write entry point for "non-field additional state" to `flbckpt` (per-rank
  `extra_<name>.bin`, or an optional callback on `write/read`), with the CA manager providing the
  serialization of `Interface`.
- On restore the data must be **in-place `assign`ed into the existing vector**, never by replacing the
  vector object, otherwise the pointer held by the boundary manager dangles; the restart path **must not
  replay `CA.Setup()`**, otherwise it re-marks the neighbourhood of nuclei as `Interface` and destroys the
  growth state. Verify by blockid + in-block index on load; out of range is a hard error.

Acceptance: `cazsblock2d` interrupted-and-restarted is bit-identical to a continuous run.

### 5.2 Runtime AMR geometry

`cazsAref2d` calls `WillRefineBlockCells()` / `GeoRefine()` at `RefineCheckStep`
(`examples/cazsAref2d/cazsArefblock2d.cpp`), so the block set changes at runtime; `BlockReader` has only
a read side, no writer (`src/io/block_reader.h`). The v2 geometry identity describes one static grid only,
`cazsAref2d` does not call `*WithGeometry` yet, and the runtime block table is not part of it.

- Have the snapshot record the complete block table (blockid / level / mesh / overlap / parent / child
  relations), or add a writer symmetric to the `BlockReader` text format; on restart rebuild the geometry
  directly from the snapshot, skipping `CheckRefine()/CreateBlocks()/AdaptiveOptimization()` and
  `CA.Setup()`.
- Refinement must keep blockids stable, otherwise the cross-rank mapping by blockid does not hold.

Acceptance: take a snapshot before and after a refinement and the restart is bit-identical to a continuous
run; restart across a different rank count is only declared supported when the block set is independent of
the rank count.

## 6. Explicitly out of scope

- Restores that are format/build incompatible: cross-architecture and cross-byte-order (raw bytes, no
  endianness marker), cross FP16/FP32 builds (the element type follows the field record, hard error only,
  no conversion), multi-version format compatibility (a `format_version` bump means a rewrite; there is no
  need to read old formats side by side).
- Persisting the `UOld/RhoOld` tolerance history: the example calls `EnableToleranceU()` after loading to
  reseed, so the first `res` after a restart is already the correct value.
- Any change to the solver's numeric path: restart belongs to the IO layer only; changing the numeric path
  would move the existing regression baselines and is unrelated to the restart goal.
