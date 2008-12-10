#
# Include file for g95 Fortran compiler on Linux
#


F90C := g95
F90FLAGS :=
LD := $(F90C)
LDFLAGS := 
MPI := on


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


ifdef MPI
  F90FLAGS += -I$(MPI_INCDIR)
  LIBS += -L$(MPI_LIBDIR) -llamf77mpi -lmpi -llam -lutil -lpthread -ldl
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -fendian=big
endif  


ifdef DEBUG
  F90FLAGS += -g -fbounds-check -ftrace=full
else
  F90FLAGS += -O3
endif

F90FLAGS += $(EXTRA_F90FLAGS)