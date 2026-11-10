TARGET = plain_comm_mpi.exe
ROOT := ../../..

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DThread_Num=1

# The coalesced plans only exist under MPI_ENABLED.  CXXC := mpic++ only swaps
# the compiler; -DMPI_ENABLED is the switch that turns on the MPI code paths.
FLAGS += -DMPI_ENABLED
CXXC := mpic++

# per-rank checksums the spawned children write; removed on success, but a
# failing run leaves them behind
OUTPUT += child_*.txt ref_*.txt

include $(ROOT)/make.mk
