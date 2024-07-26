
TARGET = genarr.exe
ROOT := ../../..

CXXC := nvcc

# FLAGS += -O3 -march=native
# FLAGS += -O3
# FLAGS += -g -DMPI_DEBUG
FLAGS += -g
# FLAGS += -pg
FLAGS += -DFLOAT_TYPE=float
# FLAGS += -DMPI_ENABLED
# FLAGS += -gencode arch=compute_89,code=sm_89
FLAGS += -gencode arch=compute_86,code=sm_86
FLAGS += --expt-relaxed-constexpr
FLAGS += -rdc=true

OUTPUT += vtkoutput

include $(ROOT)/make.mk



