
TARGET = cavblock2d.exe
ROOT := ../..

FLAGS += -O3 -Wall -march=native -fopenmp
# FLAGS += -O3 -fopenmp
# FLAGS += -O3 -march=native
# FLAGS += -g
# FLAGS += -pg
FLAGS += -DFLOAT_TYPE=double
FLAGS += -D_UNROLLFOR
# FLAGS += -DUNROLLFOR

OUTPUT += vtkoutput

include $(ROOT)/make.mk
