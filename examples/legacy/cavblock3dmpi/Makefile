
TARGET = cavblock3dmpi.exe
ROOT := ../..

FLAGS += -O3 -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs -march=native
# FLAGS += -O3 -march=native
# FLAGS += -g -DMPI_DEBUG
# FLAGS += -g
# FLAGS += -pg
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DMPI_ENABLED
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput

# use mpic++ to compile
CXXC := mpic++

include $(ROOT)/make.mk
