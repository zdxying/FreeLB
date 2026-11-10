TARGET = containers.exe
ROOT := ../../..

CXXC := nvcc

# GPU architecture: ARCH=native (default) targets the local GPU.  Override for
# cross-compilation, e.g. make ARCH=sm_89.
# Do not add -rdc=true; nothing in this target needs relocatable device code.
ARCH ?= native

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=float
FLAGS += -arch=$(ARCH)
FLAGS += --expt-relaxed-constexpr

include $(ROOT)/make.mk
