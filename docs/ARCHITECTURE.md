# FreeLB architecture and module map

Status: current (factual; no implementation phase) | Updated: 2026-10-31 | Target: let a change find its place before any code is written

> This file records the present: where things are, which way dependencies point, how data is
> organised. Judgement rules and gates are in `AGENTS.md`; the trade-offs, evidence and todos
> of a topic are in `docs/<TOPIC>_PLAN.md`.
> Exit condition: a `src/` change that makes §1 or §3 disagree with the code. The code wins;
> fix this file in the same commit.

## 0. How to read

- Which layer a change belongs to: §1.
- Where a file lives: §3.
- Changing data structures or adding a field: §2.
- Examples, build, CUDA: §4–§6.
- Where a defect or a todo goes: §7.

## 1. Layering and dependency direction

The intended order, bottom to top:

```
shared vocabulary (leaves) -> data_struct -> geometry / boundary / lbm
                           -> parallel glue -> io -> examples / tests
```

**Shared vocabulary** means the headers almost anything may include and that depend on no
business module: `src/head.h` (std/OpenMP/thrust), `utils/mpi_macros.h`, `utils/alias.h`,
`utils/tmp.h`, `data_struct/Vector.h`, `data_struct/pop_policy.h`, `parallel/comm_plan.h`
(standard library only). `data_struct/cell.h` sits here too: its file header states that
data_struct must not depend on lbm.

**The code is not a strict DAG.** `src/` contains historical cross-layer includes, typically:

- `data_struct/block_lattice.h` -> `lbm/unit_converter.h`
- `data_struct/field_struct.h` -> `parallel/communicator.h`
- `geometry/block_geometry2d.h` -> `io/block_reader.h`
- `utils/fdm_solver.h` -> `data_struct/field.h`

These are historical shape, not precedent. **New code must not add another reverse
dependency**; when you meet an existing one, check whether it is part of the
`tools/check_headers.sh` baseline, and do not widen it.

The tooling guarantees exactly one thing: **every header is self-contained** (compiles on its
own). `tools/check_headers.sh` sweeps with g++ and `mpic++ -DMPI_ENABLED`; its exit semantics
are "do not add a failure". Why the two baseline headers (`utils/cuda_device.h`,
`parallel/mpi_manager.h`) fail is stated at the top of that script. Fix one and delete its
baseline line; never add a line to silence a failure.

**Header organisation**: `.h` declares, `.hh` implements; `src/freelb.h` + `src/freelb.hh`
are the umbrella (active modules only, not the PopLattice path — see README "Umbrella
headers"). Includes are always spelled from `src/` (`#include "lbm/moment.h"`).

## 2. Data model

Blocks and cells are two layers:

```
BlockGeometry (2D/3D, owns Blocks: AABB + shapes + material/flag)
  └─ BlockLatticeManager<T, LatSet, TypePack>
       ├─ BlockLats[]          one BlockLattice per block
       │    └─ Cell<T, LatSet, FIELDS>   a per-cell view, no storage
       └─ communication plans (below)
```

`Cell` is a **view**: `get<FieldType>()` returns a reference into that block's field slot.
`FIELDS` is a `TypePack` and decides which fields a cell can see; the flag byte decides which
dynamics applies. `BlockLatticeManager::addField()` wires a `BlockFieldManager` into every
block.

Fields are organised as "one type across blocks, one instance per block":

```
BlockFieldManager<FieldType>   owns one field type across all blocks
  └─ BlockField[]              one per block, with independent storage
```

`BlockFieldManager`'s rebuild paths (`Init`, `InitAndComm`, `NonFieldInit`, `InitCopy`, the
`BlockGeometryHelper` overloads) reset the comm plans before clearing or after the swap,
depending on the overload (`field_struct.hh`). Both `BlockLatticeManager::Init()` overloads
clear `_dirPlans`/`_fullPlans`/`_allPlans`. **Any path that rebuilds blocks or fields must
invalidate the plans** for every depth — a missed site reads stale block indices: no crash,
just a wrong number.

The basic-block pointer table is the same kind of derived state: `LoadBalancing` builds it at the
call site, so every block replacement must be followed by it before `BlockGeometry3D` is constructed.

**POP storage and policy**: `POP<T,q>` wraps `FREELB_POP_ARRAY<PopStorage<T>>`, where
`PopStorage` is `__half` under FP16. The policies are `DirectPop` (recompute the address from
the offset) and `RegPop` (resolve the q addresses once in the constructor). The policy is a
template parameter, never a runtime branch.

**Communication plans** belong to the managers and are built per `(CommMode, depth)`:
`Partial` (Dir, same level), `Full`, `All`, each split by refinement depth. The population plans
live in `BlockLatticeManager`; `BlockFieldManager` keeps a second stack for scalar/vector fields
with its own tag range. `flbcomm::CommPlan` owns its `Buffer` so the send buffer outlives
`MPI_Waitall`; do not build a second buffer with a shorter lifetime.

**Unit conversion** is in `lbm/unit_converter.h`: `BaseConverter` (mechanical quantities),
`TempConverter` (thermal), `ConcConverter` (concentration), and the phase-diagram converters
`ZSConverter`/`GandinConverter` built on both; every physical/lattice conversion uses one.

## 3. Module map

Paths are relative to `src/`.

| Concern | Files |
| --- | --- |
| std/OpenMP/thrust umbrella | `head.h` |
| lattice sets (`c`/`w`/`opp`), `LatSet` | `lbm/lattice_set.h` |
| moments, equilibrium, force terms | `lbm/moment.h`, `lbm/equilibrium.h`, `lbm/force.h` |
| collision schemes (BGK/MRT/LES/force) | `lbm/collision.h`, `lbm/collisionMRT.h`, `lbm/collisionLES.h` |
| boundary conditions (bounce-back variants, `anti_pressure`) | `boundary/*.h` |
| free surface, non-Newtonian | `lbm/freeSurface.h`, `lbm/nonNewtonian.h` |
| units <-> lattice | `lbm/unit_converter.h` |
| Cell, fields, POP policy | `data_struct/cell.h`, `data_struct/field.h`, `data_struct/pop_policy.h` |
| block manager, field manager, comm plan | `data_struct/block_lattice.h`, `data_struct/field_struct.h`, `parallel/comm_plan.h` |
| field aliases (`RHO`/`VELOCITY`/`FLAG`/`POP`...) | `utils/alias.h` |
| blocks, AABB, shapes, decomposition | `geometry/block_geometry2d.h`, `geometry/block_geometry3d.h`, `geometry/basic_shape2d.h`, `geometry/basic_shape3d.h` |
| MPI runtime, `mpi()` | `parallel/mpi_manager.h`, `parallel/communicator.h` |
| `.ini`, writer family, checkpoint | `io/ini_reader.h`, `io/vtkWriter.h`, `io/vtm_writer.h`, `io/vtu_writer.h`, `io/plot_writer.h`, `io/checkpoint.h` |
| CA solidification | `ca/zhu_stefanescu2d.h`, `ca/cazs.h` |
| STL -> field | `offLattice/marchingCube.h` |
| metaprogramming (`TypePack`, `TaskSelector`) | `utils/tmp.h` |

Outside `src/`:

| Directory | Contents |
| --- | --- |
| `examples/<name>/` | one case = `main.cpp/.cu` + `Makefile` + `.ini`; reference for assembly and loop |
| `tests/<category>/<name>/` | self-contained test targets, built and run by `tests/run_tests.sh` |
| `benchmarks/<name>/` | performance benchmarks; `hidden_comm`/`plain_comm` dispatch into `2d/`/`3d/` |
| `generated/lbm/` | csegen output, never edited by hand |
| `third_party/cse/` | the CSE code generator (optional) |
| `tools/cse/` | csegen design and bridge, see `tools/cse/DESIGN.md` |

## 4. Assembly shape

Every example has this shape; the reference implementation is `examples/cavity2d`. The listing
is schematic — signatures are elided and it is not a compilable unit:

```cpp
using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
using CELL   = Cell<T, LatSet, FIELDS>;
using BulkTask = tmp::Key_TypePair<AABBFlag, collision::BGK<...>>;   // flag -> dynamics
using NSTask   = tmp::TaskSelector<tmp::TupleWrapper<...>, std::uint8_t, CELL>;
NSLattice.ApplyCellDynamics<NSTask>(FlagFM);                         // per-cell dispatch
```

Three things are separate:

- `FIELDS` decides the fields a cell can see; `CELL` is only its view type.
- `BulkTask` maps a geometry flag to dynamics; `TaskSelector` dispatches per cell on the flag byte.
- **The schedule is not part of the shape**. The *hidden step* (collision split in two around
  a non-blocking exchange, so the MPI latency hides behind the interior collision;
  `BlockLatticeManager::ApplyCellDynamicsWithHiddenComm`) runs after `Stream()` in `cavity2d`;
  `couette2d` runs `ApplyCellDynamics`/`Stream`/`NormalFullCommunicate`. They are not
  interchangeable, for the reason in `docs/PARTIAL_COMM_PLAN.md` §6; per-example call sites
  are in `docs/INTP_COMM_OVERFLOW_BUG.md` §18.

## 5. Build system

There is no top-level library build. Each target's `Makefile` ends with
`include $(ROOT)/make.mk`; the dispatch Makefiles (`benchmarks/hidden_comm`,
`benchmarks/plain_comm`) recurse into `2d/`/`3d/`, and `src/xcore/` builds its own library.

A target needs `ROOT`, and `TARGET = name.exe` on its own literal line; the silent-failure
traps around these (`-D_UNROLLFOR`, `make mpi`, CUDA arch) are `AGENTS.md` §2.

- Any directory a target writes goes in `OUTPUT += <dir>`, so `make clean` removes it.
- `-D_UNROLLFOR` is the cse switch: csegen uses `generated/lbm/*.ur.h`; without it the plain
  loop compiles with a warning (`src/ur_fallback/` is comparison-only, `make.mk:61-93`).
- xcore is optional: no target may require `libxcore.a` to exist; the fallback is README
  "Build".
- CUDA arch: `-arch=native`, or a cross-build that keeps PTX; the `-gencode`/`-rdc` traps are
  `AGENTS.md` §2.
- `make omp` / `make mpi` / `make info` / `make time` exist. Compiler flags are not tracked
  as dependencies, so switching between them needs `make clean` first.
- A new example copies the closest `Makefile`, ships its own `.ini`, and defaults to
  `Restart = 0`, `CheckpointStep = 0`.

## 6. CUDA mirror

Every host class has a `cudev::` twin: `dev_*` raw pointers plus one `dev_<Type>*` handle.
`InitDeviceData()` does malloc + `copyToDevice()` + `constructInDevice()`; kernels use the
mirror returned by `get_devObj()`, never the host object.

The mirror is `__device__`-only except for host-facing constructors, which are `__any__`
(e.g. `FieldPtrCollection`, `BlockLatticeBase`, `BlockLattice` in `cuda_block_lattice.h`,
and the same pattern in `cuda_field.h`/`cuda_field_struct.h`); it holds pointers only.

`copyToDevice()` **does not push container rotation state**: the `normalize()`/
`refreshDeviceObject()` rule is `AGENTS.md` §2 (CUDA); the stale `start[]` cached in
`cudev::CyclicArray` is why.

`check_cuda_status()` covers memcpy/malloc, **not kernel launches**; the launch/status rule
is `AGENTS.md` §2 (CUDA).

FP16 changes only the POP element type. The scope, the `RegPop` requirement and the
`DirectPop`/CSE limit are in `docs/FP16_PLAN.md`; `tests/run_tests.sh` never builds FP16, so
that path always needs a manual `make FP16=1` run of `tests/cuda/checkpoint`.

## 7. Defects and todos

Where a defect or todo goes and what it must carry is `AGENTS.md` §5; there is no global TODO.
`docs/INTP_COMM_OVERFLOW_BUG.md` is the existing defect doc (items marked fixed/unfixed with a
section number); a new defect goes into an existing doc or opens `<TOPIC>_BUG.md`.
