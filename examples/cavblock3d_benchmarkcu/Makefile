
TARGET = cavblock3d.exe
ROOT := ../..

CXXC := nvcc

FLAGS += -O3 
# FLAGS += -O3
# FLAGS += -g
# FLAGS += -pg
FLAGS += -DFLOAT_TYPE=float
FLAGS += -gencode arch=compute_89,code=sm_89
FLAGS += --expt-relaxed-constexpr
FLAGS += -rdc=true

OUTPUT += vtkoutput

# CXXC := mpic++

include $(ROOT)/make.mk
