TARGET = checkpoint_geometry_mpi.exe

ROOT := ../../..

FLAGS += -O2
FLAGS += -Wall -Wextra -Wno-unused-parameter
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DThread_Num=1
FLAGS += -DMPI_ENABLED
CXXC := mpic++
MPI_RANKS = 2

OUTPUT += ckpt

include $(ROOT)/make.mk
