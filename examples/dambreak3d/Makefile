# -------------executable name-------------
TARGET = dambreak3d.exe

# -------------/src/-------------
ROOT := ../..

# -------------optimization flags-------------
FLAGS += -O3 -march=native -mtune=native
FLAGS += -D_UNROLLFOR

# -------------warning flags-------------
FLAGS += -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs

# -------------debug flags-------------
# FLAGS += -g
# FLAGS += -DMPI_DEBUG
# FLAGS += -pg

# -------------floating point type-------------
FLAGS += -DFLOAT_TYPE=double

# -------------output folder-------------
OUTPUT += vtkoutput STLoutput

include $(ROOT)/make.mk
