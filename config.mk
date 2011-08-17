# Operating System (Linux (default),AIX or IRIX)

OS ?= Linux

# Compiler (gfortran (default), ifort, pgi, g95)

FORT ?= ifort

# Binary format (big_endian (default) or little_endian)

FORMAT ?= big_endian

# Precision (single or double (default))

PRECISION ?= single

# Suffix for executables

EXEC_SUFFIX ?= 

# uncomment/comment this line to activate/deactivate MPI support
#MPI ?= on

# uncomment/comment this line to activate/deactivate compilation with mpif90 
# Otherwise, the Fortran compiler is directly used the variables 
# MPI_INCDIR, MPI_LIBDIR and MPI_LIB.
# The flag is only used if MPI is activated.

USE_MPIF90 ?= on

# uncomment/comment this line to activate/deactivate OpenMP support

#OPENMP ?= on


# uncomment/comment this line to activate/deactivate debugger information

#DEBUG ?= on

# Library names and locations
# Variables name use the following naming scheme:
# 
# <library>_INCDIR: directory containing the include files and compiled modules (e.g. /home/user/mpi/include)
# <library>_LIBDIR: directory containing the library  (e.g. /home/user/mpi/lib)
# <library>_LIB: name of the library (including dependencies) for the linker (e.g. -lmpi -lpthread)
#

# netCDF library
# Note if the nc-config tool is present and NETCDF_CONFIG is defined, then the values of 
# NETCDF_INCDIR, NETCDF_LIBDIR and NETCDF_LIB are not used


#NETCDF_CONFIG ?= nc-config
#NETCDF_INCDIR ?= 
#NETCDF_LIBDIR ?= 
#NETCDF_LIB ?= -lnetcdf

# MPI library
# Only used if USE_MPIF90 is deactivated
#MPI_INCDIR ?= 
#MPI_LIBDIR ?= 
#MPI_LIB ?= 

# Lapack library

#LAPACK_LIBDIR ?= 
#LAPACK_LIB ?= -llapack

# BLAS library

#BLAS_LIBDIR ?=
#BLAS_LIB ?= -lblas


ifeq ($(FORT),ifort)
MPIF90 ?= /cvos/shared/apps/openmpi/intel/64/1.3.0/bin/mpif90
LDFLAGS ?= -Wl,-R/u/mc/opt/intel-11.1/netcdf-4.1.2/lib -Wl,-R/u/mc/opt/intel-11.1/lapack-3.3.1/lib/ -Wl,-R/cvos/shared/apps/gotoblas/penryn/64/1.26
NETCDF_CONFIG=/u/mc/opt/intel-11.1/netcdf-4.1.2/bin/nc-config
LAPACK_LIBDIR ?= /u/mc/opt/intel-11.1/lapack-3.3.1/lib/ 
endif

ifeq ($(FORT),pgi)
MPIF90 ?= /cvos/shared/apps/openmpi/pgi/64/1.3.0/bin/mpif90
LDFLAGS ?= -Wl,-R/u/mc/opt/pgi-9.0/netcdf-4.1.2/lib -Wl,-R/u/mc/opt/pgi-9.0/lapack-3.3.1/lib/ -Wl,-R/cvos/shared/apps/gotoblas/penryn/64/1.26
NETCDF_CONFIG=/u/mc/opt/pgi-9.0/netcdf-4.1.2/bin/nc-config
LAPACK_LIBDIR ?= /u/mc/opt/pgi-9.0/lapack-3.3.1/lib/ 
endif

ifeq ($(FORT),gfortran)
MPIF90 ?= /cvos/shared/apps/openmpi/gcc/64/1.3.0/bin/mpif90
LDFLAGS ?= -Wl,-R/u/mc/opt/netcdf-4.1.2/lib -Wl,-R/u/mc/opt/lapack-3.3.1/lib/ -Wl,-R/cvos/shared/apps/gotoblas/penryn/64/1.26
NETCDF_CONFIG=/u/mc/opt/netcdf-4.1.2/bin/nc-config
LAPACK_LIBDIR ?= /u/mc/opt/lapack-3.3.1-goto/lib/ 
endif

BLAS_LIBDIR ?= /cvos/shared/apps/gotoblas/penryn/64/1.26 
BLAS_LIB ?= -lgoto 




