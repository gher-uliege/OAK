# This file needs GNU make
#generated by: makedep -v OS=Linux;;FORT=pgi;;FORMAT=big_endian -o Makefile assim.F90

#---------#
#  Rules  #
#---------#

.f90.o:
	$(F90C) $(F90FLAGS) -c $<

.F90.o:
	$(F90C) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $<

.SUFFIXES: $(SUFFIXES) .f90 .F90

#-------------------------------#
#  Platform specific variables  #
#-------------------------------#

include config.mk

# default settings
# If you need to adapt, OS and FORT, then make the corresponding changes in config.mk

OS ?= Linux

FORT ?= gfortran

FORMAT ?= big_endian

PRECISION ?= double
USE_MPIF90 ?= on
MPI ?= 
OPENMP ?= 
DEBUG ?= 
JOBS ?= 1

include Compilers/$(OS)-$(strip $(FORT)).mk
include Compilers/libs.mk

#---------#
#  assim  #
#---------#

ASSIM_PROG ?= assim

ASSIM_SRCS = anamorphosis.F90 assim.F90 assimilation.F90 date.F90 grids.F90 \
	initfile.F90 matoper.F90 ndgrid.F90 parall.F90 rrsqrt.F90 \
	ufileformat.F90

ASSIM_OBJS = anamorphosis.o assim.o assimilation.o date.o grids.o initfile.o \
	matoper.o ndgrid.o parall.o rrsqrt.o ufileformat.o match.o

#-----------------#
#  Common macros  #
#-----------------#

PROG = $(ASSIM_PROG)

SRCS = $(ASSIM_SRCS)

OBJS = $(ASSIM_OBJS)

#-------------------#
#  Standard tagets  #
#-------------------#



all: $(PROG)

clean:
	rm -f $(PROG) $(OBJS) *.mod

lib: $(OBJS)
	ar rs liboak.a $(OBJS)

dynlib: $(OBJS)
	$(CC) -shared -Wl,-soname,liboak.so.1 -o liboak.so.1 $(OBJS) $(LIBS) $(FRTLIB)

allbin:
	make -j $(JOBS) DEBUG=on clean
	make -j $(JOBS) DEBUG=on all
	make -j $(JOBS) DEBUG=on OPENMP=on clean
	make -j $(JOBS) DEBUG=on OPENMP=on all
	make -j $(JOBS) DEBUG=on MPI=on clean
	make -j $(JOBS) DEBUG=on MPI=on all


print:
	echo $(LIBS) $(F90FLAGS) $(PRECISION)

#---------------#
#  Executables  #
#---------------#

$(ASSIM_PROG): $(ASSIM_OBJS) 
	$(F90C) $(LDFLAGS) -o $@ $(ASSIM_OBJS) $(LIBS)

#----------------#
#  Dependencies  #
#----------------#

assim.o: assim.F90 assimilation.o initfile.o matoper.o rrsqrt.o ufileformat.o \
	ppdef.h

ndgrid.o: ndgrid.F90 matoper.o ufileformat.o ppdef.h

parall.o: parall.F90 ppdef.h

matoper.o: matoper.F90 ppdef.h

date.o: date.F90 ppdef.h

anamorphosis.o: anamorphosis.F90 initfile.o ufileformat.o ppdef.h

grids.o: grids.F90 ppdef.h

initfile.o: initfile.F90 ppdef.h

rrsqrt.o: rrsqrt.F90 matoper.o parall.o ufileformat.o ppdef.h

assimilation.o: assimilation.F90 anamorphosis.o date.o grids.o initfile.o \
	matoper.o ndgrid.o parall.o rrsqrt.o ufileformat.o ppdef.h

ufileformat.o: ufileformat.F90 ppdef.h

match.o: match.c
