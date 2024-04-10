
TARGET = cavblock2dmpi.exe
ROOT := ../..

FLAGS += -O3 -march=native
# FLAGS += -O3 -fopenmp
# FLAGS += -O3 -march=native
# FLAGS += -g -DMPI_DEBUG
# FLAGS += -g
# FLAGS += -pg
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DMPI_ENABLED

OUTPUT += vtkoutput

CXXC := mpic++

include $(ROOT)/make.mk
