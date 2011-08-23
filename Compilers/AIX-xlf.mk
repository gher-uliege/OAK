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


ifdef OPENMP
  F90FLAGS += -qsmp=omp
  LDFLAGS += -qsmp=omp
endif

ifdef DEBUG
  F90FLAGS += -g -C
else
  F90FLAGS += -O3
endif


ifeq ($(FORMAT),big_endian)
  F90FLAGS +=
else
  echo "Error: machine format little_endian not available."; exit 1
endif  
