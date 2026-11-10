TARGET = containers.exe
ROOT := ../../..

CXXC := nvcc

# GPU architecture: ARCH=native (default) asks nvcc to target the local GPU.
# Override for cross-compilation or a GPU-less build host, e.g. make ARCH=sm_89.
ARCH ?= native

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=float
FLAGS += -arch=$(ARCH)
FLAGS += --expt-relaxed-constexpr
FLAGS += -rdc=true

include $(ROOT)/make.mk
