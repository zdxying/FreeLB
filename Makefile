all: lib tools

lib:
	cd src/xcore/ && make

# build the csegen source-to-source generator (needed by examples built
# with -D_UNROLLFOR; see make.mk)
tools:
	cd tools/cse/ && make

libclean:
	cd src/xcore/ && rm -rf build

tools-clean:
	cd tools/cse/ && make clean

# generate the .ur.h specializations and copy them into src/lbm/, replacing
# the hand-written versions (run tools/cse verify first)
install-ur:
	cd tools/cse/ && make install

.PHONY: all lib tools libclean tools-clean install-ur
