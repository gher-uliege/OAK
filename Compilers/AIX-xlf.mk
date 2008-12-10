#
# Include file for XLF Fortran compiler on AIX
#


F90C := xlf95_r7
F90FLAGS := -q64 -qsuffix=f=f90:cpp=F90 -WF,-DNETCDF -I$(HOME)/netcdf-3.6.0_64/include/ 
#F90FLAGS :=  -qsuffix=f=f90:cpp=F90 -WF,-DNETCDF -I$(HOME)/netcdf-3.6.0_32/include/ 
LD := $(F90C)
LDFLAGS := $(F90FLAGS)
CFLAGS=-q64


OpenMP := 
ifdef OpenMP
  F90FLAGS += -qsmp=omp
  LDFLAGS += -qsmp=omp
endif

MPI := on

ifdef MPI
  F90C=mpxlf90
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS +=
endif  

DEBUG := on
ifdef DEBUG
  F90FLAGS += -g -C
else
  F90FLAGS += -O3
endif

#
# Library locations
#

LAPACK_LIB=-L$(HOME)/lib -llapack_64
BLAS_LIB=-L/usr/lib -lblas

LIBS = -brename:.flush,.flush_  -brename:.match,.match_   -L$(HOME)/netcdf-3.6.0_64/lib -lnetcdf match.o -L$(HOME)/lib -llapack_64 -L/usr/lib -lblas

