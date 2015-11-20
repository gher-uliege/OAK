include ../config.mk
include ../Compilers/$(OS)-$(strip $(FORT)).mk
include ../Compilers/libs.mk


F90FLAGS += -Wunused-variable -Wall -I..
all: test_locassim test_ndgrid test_cellgrid locassim assimtest2 test_matoper test_rrsqrt

test_locassim: test_locassim.F90 ../matoper.o ../matoper2.o
	$(F90C) $(F90FLAGS)  -o $@ ../matoper.o ../matoper2.o test_locassim.F90  $(LIBS) $(EXTRA_LDFLAGS)

locassim: ../anamorphosis.o ../assimilation.o ../date.o ../grids.o ../initfile.o ../matoper.o ../matoper2.o ../ndgrid.o ../parall.o ../rrsqrt.o ../ufileformat.o ../match.o test_locassim.F90
	$(F90C) $(F90FLAGS)  -o $@  ../anamorphosis.o ../assimilation.o ../date.o ../grids.o ../initfile.o ../matoper.o ../matoper2.o ../ndgrid.o ../parall.o ../rrsqrt.o ../ufileformat.o ../match.o test_locassim.F90  $(LIBS) $(EXTRA_LDFLAGS)


test_ndgrid: test_ndgrid.F90
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
	rm -f test_locassim test_ndgrid *.o

test: test_locassim
	./test_locassim
	./test_cellgrid
