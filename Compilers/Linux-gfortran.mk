#
# Include file for gfortran Fortran compiler on Linux
#
# Warning avoid versions 4.0 and 4.1

F90C ?= gfortran
F90FLAGS ?= -fimplicit-none
LD ?= $(F90C)
LDFLAGS ?= 

PROFILING_F90FLAGS ?= -pg
PROFILING_LDFLAGS ?= -pg

ifdef OPENMP
  F90FLAGS += -fopenmp
  LDFLAGS += -fopenmp
endif

ifdef DEBUG
  F90FLAGS += -g -fbounds-check
else
  F90FLAGS += -O3 -ffast-math
endif

ifeq ($(PRECISION),double)
  F90FLAGS += -fdefault-real-8
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -fconvert=big-endian -frecord-marker=4
endif  


