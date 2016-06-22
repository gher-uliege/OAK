# This file needs GNU make
#generated by: makedep -v OS=Linux;;FORT=pgi;;FORMAT=big_endian -o Makefile assim.F90

#---------#
#  Rules  #
#---------#

# .f90.o:
# 	$(F90C) $(F90FLAGS) -c $<

# .F90.o:
# 	$(F90C) $(F90FLAGS) -c $<

# use -o to ensure that the object file is in the same directory that source file
%.o: %.F90
	$(F90C) $(F90FLAGS) -o $@ -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -o $@ -c $<

.SUFFIXES: $(SUFFIXES) .f90 .F90
.PHONY: test

#-------------------------------#
#  Platform specific variables  #
#-------------------------------#

include config.mk

# default settings
# If you need to adapt any of these variabels OS, FORT, ... JOBS, then make the corresponding changes in config.mk

OS ?= Linux
FORT ?= gfortran
FORMAT ?= big_endian
PRECISION ?= double
USE_MPIF90 ?= on
MPI ?= 
OPENMP ?= 
DEBUG ?= 
JOBS ?= 1

VERSION=1.6
OAK_SONAME ?= liboak.so.1
OAK_LIBNAME ?= liboak.a

include Compilers/$(OS)-$(strip $(FORT)).mk
include Compilers/libs.mk

#---------#
#  assim  #
#---------#

ASSIM_PROG ?= assim

ASSIM_SRCS = sangoma_base.f90 \
	anamorphosis.F90 assim.F90 assimilation.F90 date.F90 \
	initfile.F90 matoper.F90 covariance.F90 ndgrid.F90 parall.F90 rrsqrt.F90 \
	ufileformat.F90 random_d.f90 sangoma_ewpf.F90 user_base.f90 oak.F90

ASSIM_OBJS = anamorphosis.o assim.o assimilation.o date.o initfile.o \
	matoper.o covariance.o ndgrid.o parall.o rrsqrt.o ufileformat.o match.o sangoma_ewpf.o \
	random_d.o user_base.o oak.o

MODULES = anamorphosis.mod  assimilation.mod  date.mod initfile.mod  \
        matoper.mod covariance.mod  ndgrid.mod  parall.mod  rrsqrt.mod  ufileformat.mod oak.mod

#-----------------#
#  Common macros  #
#-----------------#

PROG = $(ASSIM_PROG)

SRCS = $(ASSIM_SRCS)

OBJS = $(ASSIM_OBJS)

#-------------------#
#  Standard tagets  #
#-------------------#



all: $(PROG) $(OAK_LIBNAME)

clean: test-clean
	rm -f $(PROG) $(OBJS) $(MODULES)

$(OAK_LIBNAME): $(OBJS)
	ar rs $(OAK_LIBNAME) $(OBJS)

$(OAK_SONAME): $(OBJS)
	@if test $$(getconf LONG_BIT) = 64 -a ! "$(PIC)"; then \
	  echo "Warning: Your system seems to be 64-bit and PIC is not activated. Creating dynamic libraries will probaly fail."; \
	  echo "Use 'make PIC=on ...' or set PIC=on in config.mk"; \
        fi 
	$(CC) -shared -Wl,-soname,$(OAK_SONAME) -o $(OAK_SONAME) $(OBJS) $(LIBS) $(FRTLIB)

allbin:
	make -j $(JOBS) FORT=$(FORT) DEBUG=on clean
	make -j $(JOBS) FORT=$(FORT) DEBUG=on all
	make -j $(JOBS) FORT=$(FORT) DEBUG=on OPENMP=on clean
	make -j $(JOBS) FORT=$(FORT) DEBUG=on OPENMP=on all
	make -j $(JOBS) FORT=$(FORT) DEBUG=on MPI=on clean
	make -j $(JOBS) FORT=$(FORT) DEBUG=on MPI=on all


# create install directories (if needed)

$(OAK_LIBDIR):
	$(MKDIR_P) $(OAK_LIBDIR)

$(OAK_INCDIR):
	$(MKDIR_P) $(OAK_INCDIR)

$(OAK_BINDIR):
	$(MKDIR_P) $(OAK_BINDIR)

install: all $(OAK_LIBDIR) $(OAK_INCDIR) $(OAK_BINDIR)
	cp $(PROG) $(OAK_BINDIR)
	if [ -e $(OAK_LIBNAME) ]; then cp $(OAK_LIBNAME) $(OAK_LIBDIR); fi
	if [ -e $(OAK_SONAME) ]; then cp $(OAK_SONAME) $(OAK_LIBDIR); fi
	cp $(MODULES) $(OAK_INCDIR)

print:
	echo $(LIBS) $(F90FLAGS) $(PRECISION) $(OAK_DIR)

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

ndgrid.o: ndgrid.F90 ndgrid_inc.F90 matoper.o ufileformat.o ppdef.h 

parall.o: parall.F90 ppdef.h

matoper.o: matoper.F90 ppdef.h matoper_inc.F90

covariance.o: matoper.o covariance.F90

date.o: date.F90 ppdef.h

anamorphosis.o: anamorphosis.F90 initfile.o ufileformat.o ppdef.h

initfile.o: initfile.F90 ppdef.h

rrsqrt.o: rrsqrt.F90 matoper.o parall.o ufileformat.o ppdef.h

sangoma_ewpf.o: random_d.o equal_weights_step.f90 quicksort.f90 gen_random.f90 subroutines_for_EWPF.f90 proposal_step.f90

user_base.o: sangoma_base.o

assimilation.o: assimilation.F90 user_base.o sangoma_base.o sangoma_ewpf.o anamorphosis.o date.o initfile.o \
	matoper.o ndgrid.o parall.o rrsqrt.o ufileformat.o ppdef.h covariance.o

ufileformat.o: ufileformat.F90 ppdef.h

match.o: match.c

oak.o: oak.F90 assimilation.o ndgrid.o

covariance.o: covariance.F90 matoper.o

# We use a single Makefile
# http://stackoverflow.com/a/1139423/3801401

# Tests
test-clean:
	rm -f test/*.mod test/*.o test/test_ndgrid test/toymodel test/test_covariance test/test_matoper
	rm -Rf test/Analysis001 test/Analysis002 test/Ens001 test/Ens002


test/toymodel.o: test/toymodel.F90 oak.o assimilation.o
test/toymodel: test/toymodel.o matoper.o covariance.o ndgrid.o assimilation.o rrsqrt.o anamorphosis.o \
               date.o parall.o initfile.o user_base.o  oak.o ufileformat.o random_d.f90 sangoma_base.f90 \
               sangoma_ewpf.o match.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $+ $(LIBS) $(EXTRA_LDFLAGS)

test/test_covariance.o: test/test_covariance.F90 matoper.o covariance.o
test/test_covariance: test/test_covariance.o matoper.o covariance.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $+ $(LIBS) $(EXTRA_LDFLAGS)

test/test_ndgrid.o: test/test_ndgrid.F90 matoper.o ndgrid.o ufileformat.o
test/test_ndgrid: test/test_ndgrid.o matoper.o ndgrid.o ufileformat.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $+ $(LIBS) $(EXTRA_LDFLAGS)

test/test_cellgrid.o: test/test_cellgrid.F90 matoper.o ndgrid.o ufileformat.o
test/test_cellgrid: test/test_cellgrid.o matoper.o ndgrid.o ufileformat.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $+ $(LIBS) $(EXTRA_LDFLAGS)

test/assimtest2.o: test/assimtest2.F90 matoper.o rrsqrt.o
test/assimtest2: test/assimtest2.o matoper.o rrsqrt.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $+ $(LIBS) $(EXTRA_LDFLAGS)

test/test_matoper.o: test/test_matoper.F90 matoper.o 
test/test_matoper: test/test_matoper.o matoper.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $+ $(LIBS) $(EXTRA_LDFLAGS)

test/test_rrsqrt.o: test/test_rrsqrt.F90  matoper.o rrsqrt.o
test/test_rrsqrt: test/test_rrsqrt.o  matoper.o rrsqrt.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $+ $(LIBS) $(EXTRA_LDFLAGS)

test: test/test_covariance test/test_ndgrid test/test_cellgrid test/assimtest2 test/test_matoper test/test_rrsqrt test/toymodel
	./test/test_ndgrid
	./test/test_toymodel
	./test/test_covariance
	./test/test_cellgrid
	./test/test_matoper

release:
	TMPOAK=$$(mktemp -d -t --suffix -OAK); \
	svn export . $$TMPOAK/OAK-$(VERSION); \
	mv $$TMPOAK/OAK-$(VERSION)/config.mk.template $$TMPOAK/OAK-$(VERSION)/config.mk; \
	tar -cvzf OAK-$(VERSION).tar.gz -C $$TMPOAK --exclude-vcs  OAK-$(VERSION)

upload:
	scp OAK-$(VERSION).tar.gz modb:/var/lib/mediawiki/upload/Alex/OAK/release/
