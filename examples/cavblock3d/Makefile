
TARGET = cavblock3d.exe
ROOT := ../..

FLAGS += -O3 -Wall -march=native -fopenmp
# FLAGS += -O3 -fopenmp
# FLAGS += -O3 -march=native
# FLAGS += -g
# FLAGS += -pg
FLAGS += -DFLOAT_TYPE=double

# FLAGS += -DMPI_ENABLED
FLAGS += -D_UNROLLFOR
# FLAGS += -DUNROLLFOR

OUTPUT += vtkoutput

# CXXC := mpic++

include $(ROOT)/make.mk
