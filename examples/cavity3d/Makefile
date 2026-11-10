# -------------executable name-------------
TARGET = cavity3d.exe

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
OUTPUT += vtkoutput
OUTPUT += cavity3d.flbgeo

include $(ROOT)/make.mk

# the example reads the geometry file at start-up
cavity3d.flbgeo: cavity3d.geo
	$(MAKE) -C $(ROOT)/tools/geometry_build
	$(ROOT)/tools/geometry_build/geometry_build.exe cavity3d.geo -o $@
$(TARGET): | cavity3d.flbgeo
