
TARGET = dambreak3d.exe
ROOT := ../..

FLAGS += -O3 -march=native
FLAGS += -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs
# FLAGS += -fopenmp
FLAGS += -DMPI_ENABLED

# FLAGS += -g

FLAGS += -DFLOAT_TYPE=double
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput STLoutput

# mpi enabled
CXXC := mpic++

include $(ROOT)/make.mk
