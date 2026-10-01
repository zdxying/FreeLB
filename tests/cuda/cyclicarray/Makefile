
TARGET = cyclicarray_cuda.exe
ROOT := ../../..

CXXC := nvcc

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=float

FLAGS += -gencode arch=compute_86,code=sm_86

FLAGS += --expt-relaxed-constexpr

include $(ROOT)/make.mk
