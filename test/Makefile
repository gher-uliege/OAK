include ../config.mk
include ../Compilers/$(OS)-$(strip $(FORT)).mk
include ../Compilers/libs.mk

# add module files from parent directory
F90FLAGS += -Wunused-variable -Wall -I..
LIBS += -L.. -loak

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $<

.SUFFIXES: $(SUFFIXES) .f90 .F90


#all: test_covariance test_ndgrid test_cellgrid locassim assimtest2 test_matoper test_rrsqrt toymodel
all: test_covariance test_ndgrid test_cellgrid assimtest2 test_matoper test_rrsqrt toymodel

TOYMODEL_OBJ=toymodel.o  ../ndgrid.o ../assimilation.o ../rrsqrt.o ../anamorphosis.o ../date.o ../parall.o ../initfile.o ../user_base.o  ../matoper.o ../oak.o ../ufileformat.o ../random_d.f90 ../sangoma_base.f90 ../sangoma_ewpf.o ../match.o 

toymodel.o: toymodel.F90 ../oak.o ../assimilation.o

toymodel: $(TOYMODEL_OBJ) 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $+ $(LIBS)

test_covariance: test_covariance.o ../matoper.o ../covariance.o
	$(F90C) $(F90FLAGS)  -o $@ test_covariance.o ../matoper.o ../covariance.o  $(LIBS) $(EXTRA_LDFLAGS)

locassim: ../anamorphosis.o ../assimilation.o ../date.o ../grids.o ../initfile.o ../matoper.o ../covariance.o ../ndgrid.o ../parall.o ../rrsqrt.o ../ufileformat.o ../match.o test_covariance.F90
	$(F90C) $(F90FLAGS)  -o $@  locassim.F90  $(LIBS) $(EXTRA_LDFLAGS)

test_ndgrid: test_ndgrid.F90 ../ndgrid.F90 ../ndgrid_inc.F90 ../matoper.F90 ../ufileformat.F90
	$(F90C) $(F90FLAGS) $(LDFLAGS)  -o $@ test_ndgrid.F90 ../ndgrid.o ../matoper.o ../ufileformat.o $(LIBS)

test_cellgrid: test_cellgrid.F90
	$(F90C) $(F90FLAGS) $(LDFLAGS)  -o $@ test_cellgrid.F90 ../ndgrid.o ../matoper.o ../ufileformat.o $(LIBS)

assimtest2: assimtest2.F90
	$(F90C) $(F90FLAGS) $(LDFLAGS)  -o $@ assimtest2.F90 ../matoper.o ../rrsqrt.o $(LIBS)

test_matoper: test_matoper.F90
	$(F90C) $(F90FLAGS) $(LDFLAGS)  -o $@ test_matoper.F90 ../matoper.o $(LIBS)

test_rrsqrt: test_rrsqrt.F90
	$(F90C) $(F90FLAGS) $(LDFLAGS)  -o $@ test_rrsqrt.F90 ../matoper.o $(LIBS)

print:
	echo $(LIBS) $(F90FLAGS) $(PRECISION)

clean:
	rm -f *.mod *.o test_ndgrid toymodel test_covariance test_matoper
	rm -Rf Analysis001 Analysis002 Ens001 Ens002

test: test_ndgrid test_covariance test_cellgrid test_matoper
	./test_ndgrid
	./test_toymodel
	./test_covariance
	./test_cellgrid
	./test_matoper
