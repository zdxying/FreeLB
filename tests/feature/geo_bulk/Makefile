TARGET = geo_bulk.exe

ROOT := ../../..

FLAGS += -O2 -g
FLAGS += -DFLOAT_TYPE=double
FLAGS += -D_VOX_ENABLED

OUTPUT += geo_bulk.flbgeo

.DEFAULT_GOAL := all
all: geo_bulk.flbgeo

geo_bulk.flbgeo: geo_bulk.geo
	$(MAKE) -C $(ROOT)/tools/geometry_build
	$(ROOT)/tools/geometry_build/geometry_build.exe geo_bulk.geo -o $@

include $(ROOT)/make.mk
