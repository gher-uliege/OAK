#
# Include file for PGI Fortran compiler on Linux
#


F90C := pgf90
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


OpenMP := 
ifdef OpenMP
  F90FLAGS += -mp
  LDFLAGS += -mp
endif

ifdef MPI
  F90FLAGS += -I$(MPI_INCDIR)
  LIBS += -L$(MPI_LIBDIR) -llamf77mpi -lmpi -llam -lutil -lpthread -ldl
endif


DEBUG := 
ifdef DEBUG
  F90FLAGS += -g -C
else
#  F90FLAGS += -u -Bstatic -fastsse -Mipa=fast
  F90FLAGS += -O3 -Bstatic -Mflushz
endif

MPI := on

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -byteswapio
endif  


# Append library locations

F90FLAGS += -I$(NETCDF_INCDIR)
LIBS += -L$(NETCDF_LIBDIR) -lnetcdf

LIBS += -L$(LAPACK_LIBDIR) -llapack -lblas
