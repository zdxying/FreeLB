# FreeLB: Freely Coupled Lattice Boltzmann Code

This code(FreeLB) aims to facilitate the free dissemination and advancement of knowledge regarding the lattice Boltzmann (LB) and cellular automata (CA) methods(may be extended to other methods in the future).

FreeLB should be helpful for beginners who know the basic of LBM but have not started to implement it yet, especially, for those who want to build their own LBM program.

The following picture shows the capability of FreeLB to couple with CA method 

![Concentration_Field](https://github.com/zdxying/FreeLB/blob/main/Concentration_Field_of_a_Solidifying_Dendrite_with_Block_Structure.png "Concentration Field of a Solidifying Dendrite with Block Structure")

The following picture shows the developing feature of adaptive mesh refinement(AMR) in FreeLB

![Streamline](https://github.com/zdxying/FreeLB/blob/main/Streamline_of_Lid_Driven_Cavity_with_Refined_Block_Structure.png "Streamline of Lid-Driven Cavity with Refined Block Structure.png")

The author is inspired by OpenLB project, but provides a more concise and clear code structure, also, a relatively direct data access pattern, which is more suitable for beginners to understand and modify the code. The following files are modified from OpenLB project: 
  > - 'src/data_struct/octree.h' 
  > - 'src/data_struct/octree.hh' 
  > - 'src/data_struct/Vector.h' 
  > - 'src/io/base64.h' 
  > - 'src/io/stlreader.h' 
  > - 'src/io/stlreader.hh' 
  > - 'src/parallel/mpi_manager.h' 
  > - 'src/parallel/mpi_manager.hh' 

This code is written in C++ for clarity and efficiency, methods of LB and CA are stored in classes in ./src folder, specificlly, the declaration of the methods are in .h files and the implementations are in .hh files. To implement the simulation, one should prepare an .ini file storing simulation parameters and write his/her own .cpp file which calls the methods of this code, also, defines the initial and boundary conditions. There are some benchmark(example) cases in ./Benchmarks folder, which can be used as a reference(Some may not compiled correctly due to changes of source files recently, but will be fixed soon).

OpenMP is enabled through the code, MPI is available for a few methods, which is still under development. Implementation of CUDA will be considered in the future.

If you have any suggestions, please contact the author by email: ymmanyuan@outlook.com

## Dependencies
- C++ compiler supporting C++17
- GNU Make \
If parallel computing is needed, install:
- OpenMP
- MPI

## License
FreeLB is licensed under the GNU GPL v3 license.

## Acknowledgement
The author would like to thank OpenLB project for the inspiration of this code and thank Prof. Tao Jing for his help in the development of this code.