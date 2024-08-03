
TARGET = cazsArefblock2d.exe
ROOT := ../..

# FLAGS += -march=native -ffast-math
# FLAGS += -O3 -march=native #-fopenmp
FLAGS += -O3 -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs -march=native -fopenmp
# FLAGS += -march=native
# FLAGS += -DThread_Num 16 -fopenmp

# unused #pragma warning
# the order of member initialization not matching the order of member declaration in the class definition
# type comparison 

# FLAGS        += -Wall -Wno-unknown-pragmas -Wno-reorder -Wno-sign-compare -Wno-unused-variable

# Uncomment to add debug symbols
# FLAGS        += -g 
# FLAGS        += -D_FLB_DEBUG
FLAGS        += -DFLOAT_TYPE=double

OUTPUT += vtkoutput

include $(ROOT)/make.mk