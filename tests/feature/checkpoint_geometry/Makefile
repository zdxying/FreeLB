TARGET = checkpoint_geometry.exe
ROOT := ../../..

# portable flags: no -march=native, no _UNROLLFOR, no FP16 -> pin the plain
# loop-based schemes and the default POP container so the bit-level comparison
# is reproducible
FLAGS += -O2
FLAGS += -Wall -Wextra -Wno-unused-parameter
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DThread_Num=1

# checkpoint files the test writes to its own directory
OUTPUT += ckpt

include $(ROOT)/make.mk
