# FreeLB: Freely Coupled Lattice Boltzmann Code

This code(FreeLB) aims to facilitate the free dissemination and advancement of knowledge regarding the lattice Boltzmann (LB) coupled with other methods/models(cellular automata (CA), free surface, non-newtonian fluid...).

FreeLB should be helpful for beginners who know the basic of LBM but have not started to implement it yet, especially, for those who want to build their own LBM program.

The following pictures show the developing feature of adaptive mesh refinement(AMR) and the capability of FreeLB to couple with CA method to simulate the solidification process of a dendrite with block structure.

![Streamline](https://github.com/zdxying/FreeLB/blob/main/Streamline_of_Lid_Driven_Cavity_with_Refined_Block_Structure.png "Streamline of Lid-Driven Cavity with Refined Block Structure.png")

![Partof_Concentration_Field](https://github.com/zdxying/FreeLB/blob/main/Part_of_Concentration_Field_of_a_Solidifying_Dendrite_with_AMR_Block_Structure.png "Part of Concentration Field of a Solidifying Dendrite with adaptive refined mesh Structure")

The author is inspired by OpenLB project, but provides a more concise and clear code structure, also, a relatively direct data access pattern, which is more suitable for beginners to understand and modify the code. The following files are modified from [OpenLB](https://www.openlb.net/) project: 
  > - 'src/data_struct/octree.h' 
  > - 'src/data_struct/octree.hh' 
  > - 'src/data_struct/Vector.h' 
  > - 'src/io/base64.h' 
  > - 'src/io/stlreader.h' 
  > - 'src/io/stlreader.hh' 
  > - 'src/parallel/mpi_manager.h' 
  > - 'src/parallel/mpi_manager.hh'
  > - 'src/lbm/freeSurface.h'

The author also used files modified from [palabos](https://palabos.unige.ch/) project:
  > - 'src/offLattice/triangleSet.h' 
  > - 'src/offLattice/marchingCube.h' 
  > - 'src/offLattice/marchingCube.hh'

This code is written in C++ for clarity and efficiency, methods of LB and CA are stored in classes in ./src folder, specificlly, the declaration of the methods are in .h files and the implementations are in .hh files. To implement the simulation, one should prepare an .ini file storing simulation parameters and write his/her own .cpp/.cu file which calls the methods of this code, also, defines the initial and boundary conditions. There are some examples in ./examples folder, which can be used as a reference.

### Umbrella headers: `freelb.h` / `freelb.hh`

`freelb.h` (declarations) and `freelb.hh` (implementations) wrap all **active**
FreeLB modules: a simulation `.cpp` normally needs only these two includes (see
`examples/cavity3d`).  Covered: the core single-phase LBM path, free surface
(`lbm/freeSurface.h`), STL geometry to field (`offLattice/marchingCube.h(.hh)`),
the BlockLattice-based CA solidification (`ca/zhu_stefanescu*.h(.hh)`, used by
`examples/cazsblock2d` / `examples/cazsAref2d`), `io/plot_writer.h` and
`utils/field_checksum.h`.

The **legacy PopLattice path** is not usable through the umbrella:
`data_struct/lattice.hh` (the PopLattice/RhoLattice method definitions),
`lbm/buoyancy.h` and `lbm/cell_dynamics.h` (a dead header with zero consumers)
are not wrapped, so legacy simulations include them directly
(`examples/legacy/*`).  What still reaches the umbrella is declaration-only:
`data_struct/lattice.h` and the PopLattice-based `CA::ZhuStefanescu2D/3D` class
definitions arrive transitively through `ca/zhu_stefanescu*.h`, which the active
BlockLattice CA path needs, and the unused `legacy::` schemes come with
`lbm/collision.h`.  All of it is pending removal together with
`examples/legacy/`.

Every header is self-contained (compiles on its own).  Two dead headers with
zero references, `data_struct/tree.h` and `utils/exception.h`, were removed
rather than wrapped.


OpenMP and MPI is enabled through the code. Cuda for multi-block structure is not supported for now

If you have any suggestions, please contact the author by email: ymmanyuan@outlook.com


## Dependencies
- C++ compiler supporting C++17
- GNU Make \
Optional(If parallel computing is needed, install):
- OpenMP
- MPI
- Cuda toolkit


## Build
To build the code, make sure you installed the above dependencies.

If you cloned without submodules, fetch them first (needed for the full build;
Quick Start works without them):
```bash
git submodule update --init third_party/cse
```

### Quick Start (no xcore/csegen required)
```bash
cd ./examples/cavity2d
make
```
The executable file will be generated in the same directory with .exe extension (please note that .exe extension is used to make the executable file recognized by .gitignore, the author did NOT test the code on Windows, but WSL is OK!)

### Full Build (with xcore memory pool + CSE code generation)
For optimal performance, build the xcore library and CSE tools first from the root directory:
```bash
make              # builds xcore library + csegen tool
cd ./examples/cavity2d
make
```
- **xcore** (`src/xcore/`): High-performance pool allocator. If not built, the build system automatically falls back to `std::allocator` (with a warning).
- **csegen** (`tools/cse/`): CSE code generator for `-D_UNROLLFOR` examples. Without `-D_UNROLLFOR`, the plain loop-based schemes compile and no `.ur.h` is touched. With `-D_UNROLLFOR`, the cse-generated specializations in `generated/lbm/` are used; if csegen is missing, the build compiles the plain schemes too (with a warning) -- the hand-written copies in `src/ur_fallback/lbm/` are a comparison reference only and never compiled.

You can also build xcore and csegen separately:
```bash
make lib          # build xcore library only
make tools        # build csegen tool only
```

### Available Make Targets
| Target | Description |
|--------|-------------|
| `make` | Build xcore + csegen (from root) |
| `make lib` | Build xcore library only |
| `make tools` | Build csegen tool only |
| `make libclean` | Clean xcore build artifacts |
| `make tools-clean` | Clean csegen build artifacts |
| `make test` | Build + run the `tests/` suite (CUDA targets auto-skipped if no `nvcc`) |
| `make test-cpu` | Same, CPU + MPI only (no CUDA toolchain / GPU required) |


## Test

`tests/` holds small, self-contained targets grouped by area
(`block/`, `cuda/`, `feature/`, `io/`, `regression/`).  Run them all with:

```bash
make test          # everything available (skips CUDA when nvcc is missing)
make test-cpu      # CPU + MPI only
tests/run_tests.sh --no-cuda --no-mpi   # finer control
```

`tests/run_tests.sh` builds and runs every `tests/<category>/<name>/` target and
checks its exit status, so a test that only prints is not enough -- it must
`return flbtest::summary(...)` (see `tests/test_util.h`) or otherwise exit
non-zero on failure.  `tests/regression/cavity2d` pins the solver's
distribution-function checksum against a baseline; changing the LBM scheme,
equilibrium or boundary treatment will (intentionally) fail it.


## Checkpoint / Restart

A run can be stopped and continued later, or resumed after a crash, without
losing progress.  The solver itself is untouched: a new io module
(`src/io/checkpoint.{h,hh}`, namespace `flbckpt`) serializes **every field of
every lattice/flag manager plus the step counter**; the geometry, boundary
conditions and task definitions are re-created by the example's own setup code,
so a restart assumes **the same parameters and the same MPI size**.

Enabled per example through three optional keys in `[Simulation_Settings]`
(absent = off, so older `.ini` files keep working):

```ini
Restart = 0            ; 0 = fresh run, 1 = resume from the newest checkpoint
CheckpointStep = 0     ; > 0: write a checkpoint every N steps (and at normal exit)
RestartDir = ./ckpt    ; checkpoint directory
```

Wired up so far: `examples/cavity2d`, `examples/cavity3d` (CPU) and
`examples/cavity2d_cu` (GPU, single block).  The CA / runtime-AMR examples
(`cazsAref2d`) are a later phase.

Resuming with a **different number of ranks** works too: the reader matches
blocks by their global `blockid` across every rank's file and restores only the
blocks it owns (`tests/block/checkpoint_decomp`).  This needs a block set that
does not depend on the rank count -- `AdaptiveOptimization(mpi().getSize())`
re-splits the domain per process count and is therefore refused on a rank-count
change; split the domain with `GeoHelper.CreateBlocks(bx, by[, bz])` +
`LoadBalancing(np)` instead.

Layout (one binary per rank, written atomically through a temporary file plus a
rename; `latest` is written last, so a partially written checkpoint is never
picked up):

```
ckpt/
  latest                  # step number of the newest complete checkpoint
  step_<s>/
    rank<r>.ck            # per-rank binary: block table + field metadata + payload
```

Restoring normalizes the population container's rotation state -- on the GPU it
also re-uploads the device mirror object, because `cudev::CyclicArray` caches
`start[]` inside the object and `copyToDevice()` does not push it -- and
re-seeds the convergence history, so a resumed run reproduces an uninterrupted
one step-for-step (`tests/feature/checkpoint` on the CPU, `tests/cuda/checkpoint`
on the GPU, both bit-exact; the latter also runs with `make FP16=1` for `__half`
POP storage).  A checkpoint taken with a different geometry, block set, or POP storage precision
is **rejected with an error** rather than silently mismatched (a different rank
count with a matching block set is accepted).


## Benchmark
The following table shows the performance of FreeLB in simulating the lid-driven cavity flow problem(/benchmarks/cavity3d and /benchmarks/cavity3d_cu) with D3Q19 lattice set on a 100x100x100 lattice grid(in million lattice updates per second, MLUPs).

| Device | Platform | FP32 MLUPs |
| :----- | :---- | :---- |
| i9-14900K       | Ubuntu 22.04      | 651 | 
| Ryzen 9 9950X (Thanks to Wu Kai) | Ubuntu 24.04 WSL2 | 645 |
| E5-2699v4       | Ubuntu 20.04 WSL2 | 553 |
| i7-11800H       | Ubuntu 22.04 WSL2 | 173 |
| i7-8750H        | Ubuntu 22.04 WSL2 | 124 |
| RTX 4080        | Ubuntu 22.04      | 2939 |
| RTX 3060 Laptop | Ubuntu 22.04 WSL2 | 1581 |


## License
FreeLB is licensed under the GNU GPL v3 license.

## Acknowledgement
The author would like to thank OpenLB project for the inspiration of this code and thank Prof. Jing Tao for his support in the development of this code.