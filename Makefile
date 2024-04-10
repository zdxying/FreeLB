# -------------rm--------------
ifeq ($(OS),Windows_NT)
    RM = del /F /Q
	DM = rmdir /S /Q
else
    RM = rm -rf
endif

OUTPUT += output

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.d)

# if not specified, use g++ as default compiler
CXXC ?= g++
# CXXC = mpic++
# FLAGS += -O3 -march=native -fopenmp
FLAGS += -std=c++17

all: $(TARGET)

# ------------target----------------

%.o: %.cpp
	$(CXXC) $(FLAGS) -fPIC -I$(ROOT)/src/ -c $< -o $@    
    # $(CXX) $(FLAGS) -MMD -c $< -o $@

$(TARGET): $(OBJS)
	$(CXXC) $(FLAGS) -o $@ $^
-include $(DEPS)
# ------------libs----------------
# LIBS = libFreeLB.a

# CPPFILES := 
# 	src/parallel/MPImanager.cpp
# OBJ_FILES := $(CPPFILES:.cpp=.o)

# # lib: 
# $(LIBS): $(OBJ_FILES)
# 	ar rc $@ $(OBJ_FILES)

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