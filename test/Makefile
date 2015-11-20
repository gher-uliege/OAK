include ../config.mk
include ../Compilers/$(OS)-$(strip $(FORT)).mk
include ../Compilers/libs.mk

# add module files from parent directory
F90FLAGS := $(F90FLAGS) -I..

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $<

.SUFFIXES: $(SUFFIXES) .f90 .F90


all: toymodel test_ndgrid

TOYMODEL_OBJ=toymodel.o  ../ndgrid.o ../assimilation.o ../rrsqrt.o ../anamorphosis.o ../date.o ../parall.o ../initfile.o ../user_base.o  ../matoper.o ../oak.o ../ufileformat.o ../random_d.f90 ../sangoma_base.f90 ../sangoma_ewpf.o ../match.o 

toymodel.o: toymodel.F90 ../oak.o ../assimilation.o

toymodel: $(TOYMODEL_OBJ) 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $+ $(LIBS)


test_ndgrid: test_ndgrid.F90 ../ndgrid.F90 ../ndgrid_inc.F90 ../matoper.F90 ../ufileformat.F90
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ test_ndgrid.F90 ../ndgrid.F90 ../matoper.F90 ../ufileformat.F90 $(LIBS)

print:
	echo $(LIBS) $(F90FLAGS) $(PRECISION)

clean:
	rm -f *.mod *.o
