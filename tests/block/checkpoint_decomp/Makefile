
TARGET = checkpoint_decomp.exe
ROOT := ../../..

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DThread_Num=1

# MPI: a cross-decomposition restore reads files written by other ranks, and
# the test re-runs itself through mpirun.
# Note: CXXC := mpic++ only swaps the compiler; -DMPI_ENABLED is the switch
# that turns on the MPI code paths
# (without it mpi() is a stub, getSize() is always 1, and both ranks write
# rank0.ck).
FLAGS += -DMPI_ENABLED
CXXC := mpic++

include $(ROOT)/make.mk
