include ../config.mk
include ../Compilers/$(OS)-$(strip $(FORT)).mk
include ../Compilers/libs.mk


F90FLAGS += -Wunused-variable -Wall
all: test_locassim test_ndgrid locassim

test_locassim: test_locassim.F90 ../matoper.o ../matoper2.o
	$(F90C) $(F90FLAGS) -I.. -o $@ ../matoper.o ../matoper2.o test_locassim.F90  $(LIBS) $(EXTRA_LDFLAGS)

locassim: ../anamorphosis.o ../assimilation.o ../date.o ../grids.o ../initfile.o ../matoper.o ../matoper2.o ../ndgrid.o ../parall.o ../rrsqrt.o ../ufileformat.o ../match.o test_locassim.F90
	$(F90C) $(F90FLAGS) -I.. -o $@  ../anamorphosis.o ../assimilation.o ../date.o ../grids.o ../initfile.o ../matoper.o ../matoper2.o ../ndgrid.o ../parall.o ../rrsqrt.o ../ufileformat.o ../match.o test_locassim.F90  $(LIBS) $(EXTRA_LDFLAGS)


test_ndgrid: test_ndgrid.F90
	$(F90C) $(F90FLAGS) $(LDFLAGS) -I.. -o $@ test_ndgrid.F90 ../ndgrid.o ../matoper.o ../ufileformat.o $(LIBS)

print:
	echo $(LIBS) $(F90FLAGS) $(PRECISION)

clean:
	rm -f test_locassim test_ndgrid *.o

check: test_locassim
	./test_locassim