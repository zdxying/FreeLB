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

ifneq ($(CXXC),nvcc)
FLAGS += -fno-diagnostics-show-template-tree
endif

# linker flags
LINKFLAGS := -L$(ROOT)/lib -lrt 
# use the following when compiling error: 
# undefined reference to `tbb::detail::r1::execution_slot(tbb::detail::d1::execution_data const*)'
# with -g -fopenmp flag
# LINKFLAGS := -L$(ROOT)/lib -lrt -ltbb

ifeq ($(CXXC),nvcc)
	LINKFLAGS += -lcuda
endif

all: $(TARGET)

# ------------target----------------
%.o: %.$(SRC_EXT)
	$(CXXC) $(FLAGS) -I$(ROOT)/src/ -c $< -o $@

$(TARGET): $(OBJS)
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