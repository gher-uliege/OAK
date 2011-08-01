#
# Include file for Intel Fortran compiler on Linux
#


F90C ?= ifort
F90FLAGS ?=
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


include Compilers/libs.mk
