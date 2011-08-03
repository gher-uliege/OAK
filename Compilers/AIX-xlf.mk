#
# Include file for XLF Fortran compiler on AIX
#


F90C ?= xlf95_r7
F90FLAGS ?= -q64 -qsuffix=f=f90:cpp=F90 -WF,-DNETCDF -I$(HOME)/netcdf-3.6.0_64/include/ 
#F90FLAGS ?=  -qsuffix=f=f90:cpp=F90 -WF,-DNETCDF -I$(HOME)/netcdf-3.6.0_32/include/ 
LD ?= $(F90C)
LDFLAGS ?= $(F90FLAGS)
CFLAGS=-q64
MPIF90=mpxlf90
EXTRA_LDFLAGS = -brename:.flush,.flush_  -brename:.match,.match_

MPI ?= 
OpenMP ?= 
DEBUG ?= 

ifdef OpenMP
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

include Compilers/libs.mk

