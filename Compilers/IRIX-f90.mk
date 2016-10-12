#
# Include file for MIPSpro Fortran compiler on IRIX
#


F90C := f90
F90FLAGS := 
LD := $(F90C)
LDFLAGS := 

ifeq ($(FORMAT),little_endian)
error:
	echo "Error: machine format little_endian not available."; exit 1
endif  

DEBUG_F90FLAGS += -g
OPTIM_F90FLAGS += -O3
