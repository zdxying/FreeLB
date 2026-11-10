
TARGET = checkpoint.exe
ROOT := ../../..

# portable flags: no -march=native, no _UNROLLFOR, no FP16 -> pin the plain
# loop-based schemes and the default POP container so the bit-level comparison
# in the test is reproducible
FLAGS += -O2
FLAGS += -Wall -Wextra -Wno-unused-parameter

# -------------floating point type-------------
FLAGS += -DFLOAT_TYPE=double

# single thread: the checksum comparison is bit-exact, so the run must be
# deterministic (Thread_Num is read from the ini in the examples, but the test
# does not use omp-parallel dynamics)
FLAGS += -DThread_Num=1

include $(ROOT)/make.mk
