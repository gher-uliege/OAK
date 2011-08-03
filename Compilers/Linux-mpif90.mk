#
# Include file for PGI Fortran compiler on Linux
#


F90C := mpif90
F90FLAGS := 
LD := $(F90C)
LDFLAGS := $(F90FLAGS)


ifdef OpenMP
  F90FLAGS += -mp
  LDFLAGS += -mp
endif


ifeq ($(FORMAT),big_endian)
  F90FLAGS += -byteswapio
endif  

DEBUG := 
ifdef DEBUG
  F90FLAGS += -g -C
else
  F90FLAGS += -O3
endif

#
# Library locations
#


F90FLAGS += -I/home/abarth/local64/netcdf-pgi/include

LIBS = -L/home/abarth/local64/netcdf-pgi/lib -lnetcdf -L/home/abarth/lib64 -llapack -lblas

