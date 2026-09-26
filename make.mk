# -------------rm--------------
ifeq ($(OS),Windows_NT)
    RM = del /F /Q
	DM = rmdir /S /Q
else
    RM = rm -rf
endif

OUTPUT ?= output

# use g++ as default c++ compiler
CXXC ?= g++

ifeq ($(CXXC),nvcc)
SRC_EXT = cu
else
SRC_EXT = cpp
endif

SRCS = $(wildcard *.${SRC_EXT})
OBJS = $(SRCS:.$(SRC_EXT)=.o)
DEPS = $(SRCS:.$(SRC_EXT)=.d)

FLAGS += -std=c++17
FLAGS += -flto

ifneq ($(CXXC),nvcc)
FLAGS += -fno-diagnostics-show-template-tree
endif

# linker flags
# LINKFLAGS := -L$(ROOT)/lib -lrt 
# use the following when compiling error: 
# undefined reference to `tbb::detail::r1::execution_slot(tbb::detail::d1::execution_data const*)'
# with -g -fopenmp flag
# LINKFLAGS := -L$(ROOT)/lib -lrt -ltbb

# xcore memory pool library (optional — falls back to std::allocator when absent)
XCORE_LIB := $(ROOT)/src/xcore/build/install/lib/libxcore.a
ifneq (,$(wildcard $(XCORE_LIB)))
    LINKFLAGS := $(XCORE_LIB)
    FLAGS += -DXCORE_ENABLED
else
    $(warning [FreeLB] xcore library not found at $(XCORE_LIB))
    $(warning [FreeLB] Building without xcore — using std::allocator fallback.)
    $(warning [FreeLB] Build xcore with: make -C $(ROOT)/src/xcore)
    LINKFLAGS :=
endif

ifeq ($(CXXC),nvcc)
	LINKFLAGS += -lcuda
endif

all: $(TARGET)

# ------------cse code generation----------------
# Examples built with -D_UNROLLFOR use the .ur.h specializations.  Headers
# listed in UR_CSE_BASES carry `// @cse` markers and are generated into
# $(GEN_DIR) by the tools/cse/csegen source-to-source translator; the
# -I$(GEN_DIR) flag below shadows src/lbm/*.ur.h for those files only
# (unlisted hand-written .ur.h files fall back to src/lbm/ untouched).
ifneq (,$(findstring -D_UNROLLFOR,$(FLAGS)))
CSEGEN := $(ROOT)/tools/cse/csegen
ifneq (,$(wildcard $(CSEGEN)))
# csegen found — enable code generation
UR_CSE_BASES ?= lbm/moment lbm/equilibrium lbm/force
GEN_DIR ?= $(ROOT)/generated
UR_GEN_FILES := $(addprefix $(GEN_DIR)/,$(UR_CSE_BASES:%=%.ur.h))
FLAGS += -I$(GEN_DIR)
else
# csegen missing — fall back to hand-written .ur.h in src/lbm/
$(warning [FreeLB] csegen not found at $(CSEGEN))
$(warning [FreeLB] Falling back to hand-written .ur.h files in src/lbm/.)
$(warning [FreeLB] Build csegen with: make -C $(ROOT)/tools/cse)
endif
endif

ifneq (,$(strip $(UR_GEN_FILES)))
$(GEN_DIR)/%.ur.h: $(ROOT)/src/%.h $(CSEGEN)
	@mkdir -p $(dir $@)
	$(CSEGEN) $< $@
endif

# ------------target----------------
%.o: %.$(SRC_EXT)
	$(CXXC) $(FLAGS) -I$(ROOT)/src/ -c $< -o $@

ifneq (,$(strip $(UR_GEN_FILES)))
# rebuild objects whenever the generated specializations change
$(OBJS): $(UR_GEN_FILES)
endif

# order-only: ensure generation runs first, but do not pass the .ur.h
# headers to the linker
$(TARGET): $(OBJS) | $(UR_GEN_FILES)
	$(CXXC) $(FLAGS) -o $@ $^ $(LINKFLAGS)
# $(CXXC) $(FLAGS) -o $@ $^ $(LDFLAGS) -lname
-include $(DEPS)
#-------------clean----------------
clean:
# rm -f $(OBJS) $(DEPS) $(TARGET)
# output folder is created by the program
	$(RM) $(OBJS) $(DEPS) $(TARGET) $(OUTPUT)
ifeq ($(OS),Windows_NT)
	$(RM) *.exe
	$(DM) $(OUTPUT)
endif
del:
	$(RM) $(OUTPUT)
ifeq ($(OS),Windows_NT)
	$(DM) $(OUTPUT)
endif

#-------------info----------------
info:
	@echo "CXXC" = $(CXXC) 
	@echo "FLAGS" = $(FLAGS) 

#-------------omp----------------
omp: FLAGS += -fopenmp
omp: all

#-------------mpi----------------
mpi: FLAGS += -DMPI_ENABLED
mpi: CXXC := mpic++
mpi: all

#-------------time----------------
time: 
	@echo "Starting timed compilation..."
	@/usr/bin/time -f "User time: %U\nSystem time: %S\nElapsed time: %E\nCPU usage: %P" $(MAKE) all
# $(MAKE) all
# $(MAKECMDGOALS)