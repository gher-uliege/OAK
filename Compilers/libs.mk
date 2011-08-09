# name of executable


ASSIM_PROG ?= assim-$(FORT)-$(PRECISION)

ifdef MPI
  ASSIM_PROG := $(ASSIM_PROG)-mpi
endif

ifdef OPENMP
  ASSIM_PROG := $(ASSIM_PROG)-openmp
endif

ifdef DEBUG
  ASSIM_PROG := $(ASSIM_PROG)-debug
endif




#
# Library locations
#

# If all libraries are in one folder

INCDIR ?= /usr/local/include
LIBDIR ?= /usr/local/lib

# netCDF configuration

NETCDF_CONFIG ?= nc-config
NETCDF_INCDIR ?= $(INCDIR)
NETCDF_LIBDIR ?= $(LIBDIR)
NETCDF_LIB ?= -lnetcdf

# MPI configuration

MPIF90 ?= mpif90
MPI_INCDIR ?= $(INCDIR)
MPI_LIBDIR ?= $(LIBDIR)
MPI_LIB ?= -llamf77mpi -lmpi -llam -lutil -lpthread -ldl

# LAPACK configuration

LAPACK_LIBDIR ?= $(LIBDIR)
LAPACK_LIB ?= -llapack

# BLAS configuration

BLAS_LIBDIR ?= $(LIBDIR)
BLAS_LIB ?= -lblas

# Extra parameters

EXTRA_F90FLAGS ?=
EXTRA_LDFLAGS ?=



# MPI

ifdef MPI
ifdef USE_MPIF90
  F90C := $(MPIF90)
  F90FLAGS += -DASSIM_PARALLEL
else
  F90FLAGS += -I$(MPI_INCDIR)
  LIBS += -L$(MPI_LIBDIR) $(MPI_LIB)
endif

ifeq ($(PRECISION),double)
  F90FLAGS += -DDEFAULT_REAL=MPI_DOUBLE_PRECISION 
else
  F90FLAGS += -DDEFAULT_REAL=MPI_REAL
endif

endif

# LAPACK

LIBS += -L$(LAPACK_LIBDIR) $(LAPACK_LIB)

# BLAS

LIBS += -L$(BLAS_LIBDIR) $(BLAS_LIB)



# netCDF library
# * use nc-config script if present (full path can be specified with 
# NETCDF_CONFIG environement variable)
# * if not use variables NETCDF_LIBDIR, NETCDF_INCDIR and NETCDF_LIB

NETCDF_VERSION := $(shell $(NETCDF_CONFIG) --version)

### check presense of nc-config script
ifeq ($(NETCDF_VERSION),)
  F90FLAGS += -I$(NETCDF_INCDIR)
  LIBS += -L$(NETCDF_LIBDIR) $(NETCDF_LIB)
else
  F90FLAGS += -I$(shell $(NETCDF_CONFIG) --includedir)
  LIBS += $(shell $(NETCDF_CONFIG) --flibs)
endif
