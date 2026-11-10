# hidden communication benchmarks (2D and 3D)
#
#   make          build both
#   make clean    remove both
#
# Each subdirectory has its own Makefile because both benchmarks define main().

.PHONY: all clean 2d 3d

all: 2d 3d

2d:
	$(MAKE) -C 2d

3d:
	$(MAKE) -C 3d

clean:
	$(MAKE) -C 2d clean
	$(MAKE) -C 3d clean
