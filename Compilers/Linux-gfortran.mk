#
# Include file for gfortran Fortran compiler on Linux
#


F90C ?= gfortran
F90FLAGS ?= -fimplicit-none
LD ?= $(F90C)
LDFLAGS ?= 

DEBUG_F90FLAGS = -g -fbounds-check

OPTIM_F90FLAGS =  -O3 -ffast-math

OPENMP_F90FLAGS = -fopenmp
OPENMP_LDFLAGS = -fopenmp

PROFILING_F90FLAGS ?= -pg
PROFILING_LDFLAGS ?= -pg

PIC_F90FLAGS=-fPIC
PIC_CFLAGS=-fPIC

# Fortran Run-Time Library
FRTLIB=-lgfortran

ifeq ($(PRECISION),double)
  F90FLAGS += -fdefault-real-8
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -fconvert=big-endian -frecord-marker=4
endif  
