#
# Include file for Intel Fortran compiler on Linux
#


F90C ?= ifort
F90FLAGS ?= -implicitnone
LD ?= $(F90C)
LDFLAGS ?= 

# avoid stack size problem
# http://software.intel.com/en-us/articles/intel-fortran-compiler-increased-stack-usage-of-80-or-higher-compilers-causes-segmentation-fault/

F90FLAGS += -heap-arrays 

DEBUG_F90FLAGS = -g -check all -traceback

OPTIM_F90FLAGS = -vec-report0 -O3

OPENMP_F90FLAGS = -openmp
OPENMP_LDFLAGS = -openmp

PROFILING_F90FLAGS ?= -p
PROFILING_LDFLAGS ?= -p

#PROFILING_F90FLAGS ?= -prof-gen
#PROFILING_LDFLAGS ?= -prof-gen

PIC_F90FLAGS=-fPIC
PIC_CFLAGS=-fPIC

ifeq ($(PRECISION),double)
  F90FLAGS += -r8
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -convert big_endian
endif  

