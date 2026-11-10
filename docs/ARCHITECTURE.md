# FreeLB architecture and module map

Status: current (factual; no implementation phase) | Updated: 2026-11-09 | Target: let a change find its place before any code is written

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
are "do not add a failure". Why the single baseline header (`utils/cuda_device.h`) fails is
stated at the top of that script. Fix it and delete its baseline line; never add a line to
silence a failure.

**Header organisation**: `.h` declares, `.hh` implements; `src/freelb.h` + `src/freelb.hh`
are the umbrella (active modules only, not the PopLattice path — see README "Umbrella
headers"). Includes are always spelled from `src/` (`#include "lbm/moment.h"`).

**Geometry record and its readers.** `geometry/voxel_geometry.h` is a domain-sized record
(role, entity, material, patches, band) depending only on `data_struct/Vector.h` and
`geometry/basic_geometry.h`: no IO, no block notion. Rasterising it
(`io/voxel_geometry_builder.h`) and persisting it (`io/geometry_file.h`) belong to IO.
`writeVoxelFlags` sits beside the record but takes the field manager as a template parameter, so
geometry reaches the field layer by instantiation, not by an include. `tools/geometry_build` is a
tool, not a `src/` module.

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

**Voxel geometry** is an input to the flag field, sized by the domain rather than by the block
decomposition:

```
VoxelGeometry --writeVoxelFlags--> BlockFieldManager --InitVoxelMap--> active map
```

The record has no block or rank notion; `writeVoxelFlags` maps it onto blocks and halos, and
`InitVoxelMap` derives the owned/active counts. Nothing writes the record back.

For the decomposition the record can also fill a `GeometryFlagField` directly
(`GeometryFlagField::Init(VoxelGeometry)`), which is the input
`BlockGeometryHelper3D::FillAwareDecomposition` counts; `writeVoxelFlags` remains the field-layer
path. The split is a helper method, not a storage policy: it accepts a cut only when the whole
split subtree lowers the dense footprint (base box plus overlap halo), so the block table is
derived from the geometry and is the same for every rank count. `GetBlockLayout()` records it in
memory and `ValidateActiveCoverage()` enforces one owner per active cell.

**POP storage and policy**: `POP<T,q>` wraps `FREELB_POP_ARRAY<PopStorage<T>>`, where
`PopStorage` is `__half` under FP16. The policies are `DirectPop` (recompute the address from
the offset) and `RegPop` (resolve the q addresses once in the constructor). The policy is a
template parameter, never a runtime branch.

**Communication plans** belong to the managers and are built per `(CommMode, depth)`:
`Partial` (Dir, same level), `Full`, `All`, each split by refinement depth. The population plans
live in `BlockLatticeManager`; `BlockFieldManager` keeps a second stack for scalar/vector fields
with its own tag range. `flbcomm::CommPlan` owns its `Buffer` so the send buffer outlives
`MPI_Waitall`; do not build a second buffer with a shorter lifetime.

**Plan stamps**: the generation a plan was built against lives on `flbcomm::CommPlan`, not on the
manager. One manager owns plans for several depths and modes, and a manager-wide value let one
build vouch for the others. The generation advances on block-set creation/rebuild and on
`InitVoxelMap`; every exchange checks its own plan. Runtime free-surface reclassification is not a
geometry rebuild and does not advance it.

**One message per peer, receives first.** A plan merges every link to one peer rank into one
message, so the message count is `O(peer rank)`, not `O(block x neighbour)`. `runPlan()` and the
hidden step post every receive before packing and posting sends: the handshake overlaps the pack,
and a receive that is already posted cannot depend on a send completing. Packing runs per group,
and group buffers are disjoint, so groups are the unit of parallelism; above a group-count/size
threshold pack and unpack run under OpenMP, with every MPI call kept outside the parallel region.
Below the threshold receives are waited per group, so the first group can be unpacked before the
last arrives. The hidden step schedules adjacent-cell collision, post-receive + pack/send, interior
collision, wait for sends, same-rank halo copy, then wait/unpack for remote halos.

**Tags: fixed slots, receiver-derived order.** Every exchange uses one of seven constant slots in
`parallel/comm_plan.h`: `DirTag/FullTag/AllTag/AverTag/IntpTag/AllAverTag/AllIntpTag`. A
`BlockFieldManager` draws a per-manager base from `allocateFieldTagBase()` above
`FieldTagBase = 1 << 14`, and its normal/all messages use base +0/+1. Tags only have to be unique
among messages in flight on one communicator between one pair of ranks, so a constant per exchange
kind is enough. Order *inside* a coalesced message is a `flbcomm::LinkKey{Receiver, Sender, Kind,
Index}` (`parallel/link_key.h`), a lexicographic total order that both sides derive from geometry
alone. `checkTagBudget()` therefore only checks that `MPI_TAG_UB` covers the reserved field range.

- Rejected: **block-id tags**. They couple tag width to block-id width, and two field managers
  exchanging concurrently can collide; the same argument rejects `MPI_Comm_dup` per subsystem
  (extra communicator lifetime for no matching benefit).
- Rejected: **splitmix64 hash for `Order`**. Its only gain was dropping a width assumption that
  the tag check already enforces, at the cost of a silent collision surface. A 64-bit packed key
  would be smaller but ties the sort key to `FieldTagBase`'s width; four fields need no bit budget.

The premise is that every rank constructs the same field-manager set in the same order. Violating
it hangs or errors; it does not silently mismatch data.

**MPI thread level and the hybrid build.** `MpiManager::init()` calls
`MPI_Init_thread(..., MPI_THREAD_FUNNELED, &provided)` and aborts unless the returned level is at
least `MPI_THREAD_FUNNELED` (`parallel/mpi_manager.hh`). Every MPI call is made by the thread that
enters the step; OpenMP wraps only the pack/unpack bodies, which contain no MPI calls. `make hybrid`
in `make.mk` is the first-class `mpic++ + -DMPI_ENABLED + -fopenmp` target; `make mpi` and
`make omp` remain the single-mode builds.

**Open maintenance items.** The sub-cycling level test is still written inline as
`count % static_cast<int>(std::pow(2, deLevel))` at each call site
(`data_struct/block_lattice.hh`), while plans are selected through `flbcomm::activeDepth()`.
There is no single `isActive(level, count, role)`, so changing the filter means touching every
site. Separately, `BasicCommSet`/`LatticeCommSet` still live in `parallel/communicator.h` beside
the plan stack: they serve the periodic boundary (`boundary/periodic_boundary.h`), a local-only
path with no MPI exchange. They should either fold into the plan structures or be marked
explicitly legacy/local-only; `tools/check_headers.sh` and `block/refblock` are the gate either
way.

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
| blocks, AABB, shapes, decomposition (`RCB`, `FillAwareDecomposition`), active mask (`InitVoxelMap`) | `geometry/block_geometry2d.h`, `geometry/block_geometry3d.h`, `geometry/basic_shape2d.h`, `geometry/basic_shape3d.h` |
| voxel geometry record (role/entity/material, patches, band) -> flag mask | `geometry/voxel_geometry.h`, `geometry/voxel_geometry.hh` |
| STL/primitive -> `VoxelGeometry` | `io/voxel_geometry_builder.h`, `io/voxel_geometry_builder.hh` |
| `.flbgeo` geometry file and its build recipe | `io/geometry_file.h`, `io/geometry_file.hh` |
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
| `tools/geometry_build/` | `.geo` recipe -> `.flbgeo`; its own `Makefile`, no library target |
| `tools/casegen/` | `.case` -> case scaffold, incl. `[geometry]`/`[geometry.source.*]`; see `tools/casegen/README.md` |

The geometry path is pinned by `tests/io/geometry_file`, `feature/voxel_entity`,
`feature/voxel_mask`, `feature/casegen_geo_box3d`, and `block/compact_mpi`; the Python
generator has no C++ target and is run by `make test-casegen`.

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
  interchangeable, for the reason in `docs/INTP_COMM_OVERFLOW_BUG.md` §18.9; per-example call
  sites are in §18.

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

**Not mirrored**: `flbgeo::VoxelGeometry`, `VoxelGeometryBuilder` and `GeometryFile` are
host-only. They have no `cudev::` twin and no device copy, so a CUDA case builds its flag field on
the host and must not assume the record reaches a kernel. The device side arrives with the storage
decision (`feat/active-cuda`, B6c, in `docs/COMPLEX_GEOMETRY_PLAN.md`).

FP16 changes only the POP element type. The scope, the `RegPop` requirement and the
`DirectPop`/CSE limit are in `docs/FP16_PLAN.md`; `tests/run_tests.sh` never builds FP16, so
that path always needs a manual `make FP16=1` run of `tests/cuda/checkpoint`.

## 7. Defects and todos

Where a defect or todo goes and what it must carry is `AGENTS.md` §5; there is no global TODO.
`docs/INTP_COMM_OVERFLOW_BUG.md` is the existing defect doc (items marked fixed/unfixed with a
section number); a new defect goes into an existing doc or opens `<TOPIC>_BUG.md`.
