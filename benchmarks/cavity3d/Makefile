# -------------executable name-------------
TARGET = cavity3d.exe

# -------------/src/-------------
ROOT := ../..

# -------------optimization flags-------------
FLAGS += -O3 -march=native -mtune=native
FLAGS += -D_UNROLLFOR
FLAGS += -DSingleBlock_OMP

# -------------warning flags-------------
FLAGS += -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs

# -------------debug flags-------------
# FLAGS += -g
# FLAGS += -pg

# -------------floating point type-------------
FLAGS += -DFLOAT_TYPE=float

# -------------output folder-------------
OUTPUT += vtkoutput

include $(ROOT)/make.mk