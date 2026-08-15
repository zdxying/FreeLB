all:
	cd src/xcore/ && make

libclean:
	cd src/xcore/ && rm -rf build

.PHONY: all libclean