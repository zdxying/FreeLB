
TARGET = dambreak2d.exe
ROOT := ../..

FLAGS += -O3 -march=native -fopenmp
FLAGS += -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs
# FLAGS += -O3 -march=native
# FLAGS += -march=native
# FLAGS += -fopenmp
# FLAGS += -fsanitize=address
# FLAGS += -g
FLAGS += -DFLOAT_TYPE=double
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput

include $(ROOT)/make.mk
