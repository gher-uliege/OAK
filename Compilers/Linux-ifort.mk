#
# Include file for Intel Fortran compiler on Linux
#


F90C ?= ifort
F90FLAGS := 
LD := $(F90C)
LDFLAGS := 

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -convert big_endian
endif  

ifdef DEBUG
  F90FLAGS += -g
else
  F90FLAGS +=  
endif

# If all libraries are in one folder

INCDIR ?= /usr/local/include
LIBDIR ?= /usr/local/lib

NETCDF_INCDIR ?= $(INCDIR)
NETCDF_LIBDIR ?= $(LIBDIR)

MPI_INCDIR ?= $(INCDIR)
MPI_LIBDIR ?= $(LIBDIR)

LAPACK_LIBDIR ?= $(LIBDIR)

EXTRA_F90FLAGS ?=
EXTRA_LDFLAGS ?=


#
# Library locations
#

F90FLAGS += -I$(NETCDF_INCDIR)
LIBS += -L$(NETCDF_LIBDIR) -lnetcdf

LIBS += -L$(LAPACK_LIBDIR) -llapack -lblas
