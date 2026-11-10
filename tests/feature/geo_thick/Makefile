TARGET = geo_thick.exe

ROOT := ../../..

FLAGS += -O2 -g
FLAGS += -DFLOAT_TYPE=double
FLAGS += -D_VOX_ENABLED

OUTPUT += geo_thick.flbgeo

.DEFAULT_GOAL := all
all: geo_thick.flbgeo

geo_thick.flbgeo: geo_thick.geo
	$(MAKE) -C $(ROOT)/tools/geometry_build
	$(ROOT)/tools/geometry_build/geometry_build.exe geo_thick.geo -o $@

include $(ROOT)/make.mk
