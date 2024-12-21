
TARGET = pipeinjection3d.exe
ROOT := ../..

FLAGS += -O3 -march=native -fopenmp
FLAGS += -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs
# FLAGS += -O3 -march=native
# FLAGS += -march=native
# FLAGS += -fopenmp
# FLAGS += -g
FLAGS += -DFLOAT_TYPE=double
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput  STLoutput

include $(ROOT)/make.mk
