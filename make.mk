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

FLAGS += -std=c++17

# linker flags
LINKFLAGS := -L$(ROOT)/lib

all: $(TARGET)
	
# ------------target----------------
%.o: %.cpp
	$(CXXC) $(FLAGS) -fPIC -I$(ROOT)/src/ -c $< -o $@

$(TARGET): $(OBJS)
	$(CXXC) $(FLAGS) -o $@ $^ 
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