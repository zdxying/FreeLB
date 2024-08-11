
TARGET = cavblock3d.exe
ROOT := ../..

FLAGS += -O3 -march=native -fopenmp
# FLAGS += -O3 -fopenmp
# FLAGS += -O3 -march=native
# FLAGS += -g
# FLAGS += -pg
FLAGS += -DFLOAT_TYPE=double -DSingleBlock_OMP
FLAGS += -D_UNROLLFOR
# FLAGS += -DUNROLLFOR

# FLAGS += -DMPI_ENABLED

OUTPUT += vtkoutput

# CXXC := mpic++

include $(ROOT)/make.mk
