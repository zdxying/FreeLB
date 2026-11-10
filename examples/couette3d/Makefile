TARGET = couette3d.exe

ROOT := ../..

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=double

OUTPUT += vtkoutput
OUTPUT += couette3d.flbgeo

include $(ROOT)/make.mk

couette3d.flbgeo: couette3d.geo
	$(MAKE) -C $(ROOT)/tools/geometry_build
	$(ROOT)/tools/geometry_build/geometry_build.exe couette3d.geo -o $@
$(TARGET): | couette3d.flbgeo
