#
# Include file for Intel Fortran compiler on Linux
#


F90C ?= ifort
F90FLAGS ?= -implicitnone
LD ?= $(F90C)
LDFLAGS ?= 

# avoid stack size problem, but can be performance hit (7.8 times slower in the case of test/test_nondiag revision 7003)
# http://software.intel.com/en-us/articles/intel-fortran-compiler-increased-stack-usage-of-80-or-higher-compilers-causes-segmentation-fault/

#F90FLAGS += -heap-arrays 

DEBUG_F90FLAGS = -g -fbounds-check

OPTIM_F90FLAGS = -vec-report0 -O3

OPENMP_F90FLAGS = -openmp
OPENMP_LDFLAGS = -openmp

PROFILING_F90FLAGS ?= -profile-functions -profile-loops=all -profile-loops-report=2
PROFILING_LDFLAGS ?=  -profile-functions -profile-loops=all -profile-loops-report=2

PROFILING_F90FLAGS ?= -p
PROFILING_LDFLAGS ?= -p

PIC_F90FLAGS=-fPIC
PIC_CFLAGS=-fPIC

ifeq ($(PRECISION),double)
  F90FLAGS += -r8
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -convert big_endian
endif  

