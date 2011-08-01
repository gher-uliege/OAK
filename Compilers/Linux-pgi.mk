#
# Include file for PGI Fortran compiler on Linux
#


F90C := pgf90
F90FLAGS :=
LD := $(F90C)
LDFLAGS := 
MPI := on

OpenMP := 
ifdef OpenMP
  F90FLAGS += -mp
  LDFLAGS += -mp
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

include Compilers/libs.mk
