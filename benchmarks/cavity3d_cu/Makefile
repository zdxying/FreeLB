TARGET = cavity3d.exe
ROOT := ../..

CXXC := nvcc

# ---------------------------------------------------------------------------
# 2026-09-26: this target did not run at all on an sm_86 part.
#
# `FLAGS += -gencode arch=compute_89,code=sm_89` produced a binary whose only
# embedded image is sm_89 SASS with no PTX (confirmed with
# `cuobjdump --list-elf` / `--list-ptx`), so on an sm_86 GPU (RTX 3060 Laptop)
# no kernel could ever be launched.  `-rdc=true` compounded it: make.mk links
# -lcuda but not -lcudadevrt, so every launch failed with "named symbol not
# found".
#
# Neither failure was visible, because the launches were never checked: the run
# printed "Calculation Complete!" and a plausible-looking
# "Time Elapsed: 0.001 s" for 1000 steps.  cavity3d.cu now checks
# cudaGetLastError() after the main loop.
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
# Distribution-function container for the POP field.  Same solver either way;
# useful for comparing cache behaviour of the different cudev arrays.
#   make POP=cyclic      (default) CyclicArray
#   make POP=streammap             StreamMapArray
# ---------------------------------------------------------------------------
# FP16 POP storage (FP32 compute): traffic halves, ~2x MLUPS
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
