# -------------rm--------------
ifeq ($(OS),Windows_NT)
    RM = del /F /Q
	DM = rmdir /S /Q
else
    RM = rm -rf
endif

OUTPUT ?= vtkoutput

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

# linker flags
LINKFLAGS := -L$(ROOT)/lib -lrt
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