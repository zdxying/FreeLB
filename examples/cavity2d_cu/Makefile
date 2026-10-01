TARGET = cavity2d.exe
ROOT := ../..

CXXC := nvcc

# ---------------------------------------------------------------------------
# 2026-09-27: same failure mode as benchmarks/cavity3d_cu (fixed there first).
#
# `FLAGS += -gencode arch=compute_89,code=sm_89` produced a binary whose only
# embedded image is sm_89 SASS with no PTX, so on an sm_86 GPU (RTX 3060
# Laptop) no kernel could ever be launched.  `-rdc=true` compounded it: make.mk
# links -lcuda but not -lcudadevrt, so the first D2H in the VTK path failed
# with "named symbol not found" (the launch errors themselves were silent --
# launches were never checked before cavity3d.cu got cudaGetLastError()).
#
# -arch=native lets nvcc pick for whatever device is present; override with
# e.g. `make ARCH=86` to cross-build.  Do NOT re-add -rdc without also adding
# -lcudadevrt.
# ---------------------------------------------------------------------------
ARCH ?= native
ifeq ($(ARCH),native)
FLAGS += -arch=native
else
FLAGS += -gencode arch=compute_$(ARCH),code=sm_$(ARCH)
endif

FLAGS += -O3
FLAGS += -DFLOAT_TYPE=float
FLAGS += -D_UNROLLFOR
FLAGS += --expt-relaxed-constexpr

# ---------------------------------------------------------------------------
# Distribution-function container for the POP field.
#   make POP=cyclic      (default) CyclicArray
#   make POP=streammap             StreamMapArray
# ---------------------------------------------------------------------------
# FP16 POP storage (FP32 compute): traffic halves, ~2x MLUPS.
# NOTE: FP16 is only verified on the register-resident path; the baseline
# PopRef path produces NaN with FP16 storage (see docs/FP16_PLAN.md).
FP16 ?= 0
ifeq ($(FP16),1)
FLAGS += -DFREELB_POP_FP16
endif

POP ?= cyclic
ifeq ($(POP),streammap)
FLAGS += -DFREELB_POP_STREAMMAP
endif
ifeq ($(POP),cyclic)
FLAGS += -DFREELB_POP_CYCLIC
endif

OUTPUT += vtkoutput

include $(ROOT)/make.mk
