#
# Include file for PGI Fortran compiler on Linux
#


F90C := pgf90
F90FLAGS :=  -I/home/abarth/local64/netcdf-pgi/include/
LD := $(F90C)
LDFLAGS := 

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -byteswapio
endif  

ifdef DEBUG
  F90FLAGS += -g
else
  F90FLAGS += -u -Bstatic -O1  -r8
endif

#
# Library locations
#

INCLUDES = -I/home/abarth/local64/netcdf-pgi/include/

#DINEOF_LIBRARIES =  -L$(HOME)/lib -larpack_x86_64  -llapack  -lblas  -lnetcdf
DINEOF_LIBRARIES =  -L$(HOME)/lib -larpack_x86_64  -llapack  -lblas -L/home/abarth/local64/netcdf-pgi/lib -lnetcdf

CROSSVAL_LIBRARIES =  -L/home/abarth/local64/netcdf-pgi/lib -lnetcdf
