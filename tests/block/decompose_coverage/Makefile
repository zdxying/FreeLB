TARGET = decompose_coverage.exe

ROOT := ../../..

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DMPI_ENABLED
FLAGS += -D_VOX_ENABLED
CXXC := mpic++

# run_tests.sh reads this line and runs the target with mpirun -np 2.
MPI_RANKS = 2

include $(ROOT)/make.mk
