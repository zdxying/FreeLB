TARGET = instantiation_smoke.exe
ROOT := ../../..

FLAGS += -DFLOAT_TYPE=double
FLAGS += -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs

# no vtk output needed for this compile-only smoke test
OUTPUT +=

include $(ROOT)/make.mk

# ---------------------------------------------------------------------------
# MPI build: compiles the same TU with -DMPI_ENABLED so the MPI-guarded code in
# data_struct/block_lattice.hh is type-checked (the explicit instantiation of
# DynamicBlockLatticeHelper2D lives behind #ifdef MPI_ENABLED).  Kept separate
# from the default CPU-only target so no MPI toolchain is required normally.
# ---------------------------------------------------------------------------
MPI_TARGET = instantiation_smoke_mpi.exe
MPI_OBJ    = instantiation_mpi.o

mpi: $(MPI_TARGET)

$(MPI_OBJ): instantiation.cpp
	mpic++ $(FLAGS) -DMPI_ENABLED -I$(ROOT)/src/ -MMD -MP -MF $(@:.o=.d) -c $< -o $@

$(MPI_TARGET): $(MPI_OBJ)
	mpic++ $(FLAGS) -DMPI_ENABLED -o $@ $^ $(LINKFLAGS)

clean-mpi:
	$(RM) $(MPI_OBJ) $(MPI_OBJ:.o=.d) $(MPI_TARGET)

-include instantiation_mpi.d

.PHONY: mpi clean-mpi

