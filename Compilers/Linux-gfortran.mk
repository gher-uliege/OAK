#
# Include file for gfortran Fortran compiler on Linux
#
# Warning avoid versions 4.0 and 4.1

F90C := gfortran
F90FLAGS :=
LD := $(F90C)
LDFLAGS := 
MPI := on


ifdef MPI
#  F90FLAGS += -I$(MPI_INCDIR)
#  LIBS += -L$(MPI_LIBDIR) -llamf77mpi -lmpi -llam -lutil -lpthread -ldl
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -fconvert=big-endian -frecord-marker=4
endif  


ifdef DEBUG
  F90FLAGS += -g -fbounds-check
else
  F90FLAGS += -O3 -ffast-math
endif


include Compilers/libs.mk
