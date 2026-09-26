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

### Quick Start (without CSE optimization)
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
- **csegen** (`tools/cse/c`): CSE code generator for `-D_UNROLLFOR` examples. If not built, the build system falls back to hand-written `.ur.h` files in `src/lbm/` (with a warning).

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
| `make install-ur` | Generate and install `.ur.h` specializations into `src/lbm/` |
| `make libclean` | Clean xcore build artifacts |
| `make tools-clean` | Clean csegen build artifacts |


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