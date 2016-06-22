#
# Include file for g95 Fortran compiler on Linux
#


F90C ?= g95
F90FLAGS ?=
LD ?= $(F90C)
LDFLAGS ?= 

# OpenMP is not availble

DEBUG_F90FLAGS += -g -fbounds-check -ftrace=full

OPTIM_F90FLAGS += -O3

ifeq ($(PRECISION),double)
  F90FLAGS += -fdefault-real-8
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -fendian=big
endif  

