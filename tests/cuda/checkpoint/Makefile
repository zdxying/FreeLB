
TARGET = checkpoint.exe
ROOT := ../../..

# GPU architecture: ARCH=native (default) targets the local GPU.
ARCH ?= native

CXXC := nvcc

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=float
FLAGS += -arch=$(ARCH)
FLAGS += --expt-relaxed-constexpr
FLAGS += -DThread_Num=1

# POP storage: default FP32 (CyclicArray<float>).  FP16=1 builds the same test
# against __half storage, which is the other branch of the checkpoint element
# tag (ElemTag::half) -- run both.
FP16 ?= 0
ifeq ($(FP16),1)
FLAGS += -DFREELB_POP_FP16
endif

include $(ROOT)/make.mk
