#
# Include file for gfortran Fortran compiler on Linux
#
# Warning avoid versions 4.0 and 4.1

F90C := gfortran
F90FLAGS :=
LD := $(F90C)
LDFLAGS := 
MPI := on


#
# Library locations
#

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


F90FLAGS += -I$(NETCDF_INCDIR)
LIBS += -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf

LIBS += -L$(LAPACK_LIBDIR) -llapack -lblas


ifdef MPI
#  F90FLAGS += -I$(MPI_INCDIR)
#  LIBS += -L$(MPI_LIBDIR) -llamf77mpi -lmpi -llam -lutil -lpthread -ldl
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -fconvert=big-endian -frecord-marker=4
endif  


ifdef DEBUG
  F90FLAGS += -g -fbounds-check
else
  F90FLAGS += -O3 -ffast-math
endif
