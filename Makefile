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

# build + run every tests/<category>/<name> target, checking exit codes.
# CUDA targets are skipped automatically when nvcc is unavailable.
test:
	bash tests/run_tests.sh

# CPU + MPI only (safe on machines without a CUDA toolchain / GPU).
test-cpu:
	bash tests/run_tests.sh --no-cuda

.PHONY: all lib tools libclean tools-clean test test-cpu
