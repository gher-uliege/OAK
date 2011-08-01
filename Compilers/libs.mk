

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




LIBS += -L$(LAPACK_LIBDIR) -llapack -lblas


# netCDF library
# * use nc-config script if present (full path can be specified with 
# NC_CONFIG environement variable)
# * if not use variables NETCDF_LIBDIR, NETCDF_INCDIR and NETCDF_LIB

NC_CONFIG ?= nc-config
NETCDF_VERSION := $(shell $(NC_CONFIG) --version)

### check presense of nc-config script
ifeq ($(NETCDF_VERSION),)
  NETCDF_INCDIR ?= $(INCDIR)
  NETCDF_LIBDIR ?= $(LIBDIR)
  NETCDF_LIB ?= -lnetcdf

  F90FLAGS += -I$(NETCDF_INCDIR)
  LIBS += -L$(NETCDF_LIBDIR) $(NETCDF_LIB)
else
  F90FLAGS += $(shell $(NC_CONFIG) --fflags)
  LIBS += $(shell $(NC_CONFIG) --flibs)
endif
