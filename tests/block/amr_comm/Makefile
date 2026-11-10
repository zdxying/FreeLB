TARGET = amr_comm.exe
ROOT := ../../..

FLAGS += -O2
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DThread_Num=1
FLAGS += -DMPI_ENABLED
CXXC := mpic++

# per-rank checksums the spawned children write; removed on success, but a
# failing run leaves them behind
OUTPUT += child_*.txt

include $(ROOT)/make.mk
