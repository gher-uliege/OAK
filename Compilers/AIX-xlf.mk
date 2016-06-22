#
# Include file for XLF Fortran compiler on AIX
#


F90C ?= xlf95_r7
F90FLAGS ?= -q64 -qsuffix=f=f90:cpp=F90 -WF,-DNETCDF
LD ?= $(F90C)
LDFLAGS ?= $(F90FLAGS)
CFLAGS=-q64
MPIF90=mpxlf90
EXTRA_LDFLAGS = -brename:.flush,.flush_  -brename:.match,.match_

DEBUG_F90FLAGS += -g -C
OPTIM_F90FLAGS += -O3

OPENMP_F90FLAGS += -qsmp=omp
OPENMP_LDFLAGS += -qsmp=omp

ifeq ($(FORMAT),big_endian)
  F90FLAGS +=
else
  echo "Error: machine format little_endian not available."; exit 1
endif  
