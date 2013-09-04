include ../config.mk
include ../Compilers/$(OS)-$(strip $(FORT)).mk
include ../Compilers/libs.mk

all: test_locassim test_ndgrid

test_locassim: test_locassim.F90 ../matoper.o ../matoper2.o
	$(F90C) $(F90FLAGS) -I.. -o $@ ../matoper.o ../matoper2.o test_locassim.F90  $(LIBS) $(EXTRA_LDFLAGS)


test_ndgrid: test_ndgrid.F90
	$(F90C) $(LDFLAGS) -I.. -o $@ test_ndgrid.F90 ../ndgrid.o ../matoper.o ../ufileformat.o $(LIBS)

print:
	echo $(LIBS) $(F90FLAGS) $(PRECISION)

clean:
	rm -f test_locassim test_ndgrid

check: all
	./test_locassim