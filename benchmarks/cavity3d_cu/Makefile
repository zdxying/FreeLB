
TARGET = cavity3d.exe
ROOT := ../..

CXXC := nvcc

FLAGS += -Xptxas
FLAGS += -O3
# FLAGS += -use_fast_math
# FLAGS += -g
# FLAGS += -pg
FLAGS += -DFLOAT_TYPE=float

FLAGS += -gencode arch=compute_89,code=sm_89
# FLAGS += -gencode arch=compute_86,code=sm_86

FLAGS += -D_UNROLLFOR
FLAGS += --expt-relaxed-constexpr
FLAGS += -rdc=true

OUTPUT += vtkoutput

include $(ROOT)/make.mk
