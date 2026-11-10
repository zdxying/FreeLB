TARGET = cavity2d.exe
ROOT := ../..

CXXC := nvcc

# ---------------------------------------------------------------------------
# CUDA build arch.
#
# Correct: -arch=native, or a cross-build that keeps PTX in the image:
#     -gencode arch=compute_$(ARCH),code=[sm_$(ARCH),compute_$(ARCH)]
# Wrong:   -gencode arch=compute_$(ARCH),code=sm_$(ARCH)   (SASS only; the
#          binary then fails on any other device with "no kernel image").
# Wrong:   -rdc=true.  Nothing here needs relocatable device code.
# Wrong:   unchecked launches.  Call cudaGetLastError() after the main loop
#          (see the end of main), or a dead kernel run still prints
#          "Calculation Complete!".
# ---------------------------------------------------------------------------
ARCH ?= native
ifeq ($(ARCH),native)
FLAGS += -arch=native
else
FLAGS += -gencode arch=compute_$(ARCH),code=[sm_$(ARCH),compute_$(ARCH)]
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
# PopRef path produces NaN with FP16 storage.
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
