TARGET = normal_comm_correctness.exe
ROOT := ../../..

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DThread_Num=1

# The hidden path only exists under MPI_ENABLED, even for a single process.
# Note: CXXC := mpic++ only swaps the compiler; -DMPI_ENABLED is the switch that
# turns on the MPI code paths.
FLAGS += -DMPI_ENABLED
CXXC := mpic++

include $(ROOT)/make.mk
