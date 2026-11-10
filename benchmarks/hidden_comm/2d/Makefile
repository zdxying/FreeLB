TARGET = hidden_comm_bench2d.exe
ROOT := ../../..

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DThread_Num=1
FLAGS += -DMPI_ENABLED
FLAGS += -D_HIDDEN_PROFILE
CXXC := mpic++

include $(ROOT)/make.mk
