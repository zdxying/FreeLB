TARGET = geo_thin.exe

ROOT := ../../..

FLAGS += -O2 -g
FLAGS += -DFLOAT_TYPE=double
FLAGS += -D_VOX_ENABLED

OUTPUT += geo_thin.flbgeo

.DEFAULT_GOAL := all
all: geo_thin.flbgeo

geo_thin.flbgeo: geo_thin.geo
	$(MAKE) -C $(ROOT)/tools/geometry_build
	$(ROOT)/tools/geometry_build/geometry_build.exe geo_thin.geo -o $@

include $(ROOT)/make.mk
