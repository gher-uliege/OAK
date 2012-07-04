include ../config.mk
include ../Compilers/$(OS)-$(strip $(FORT)).mk
include ../Compilers/libs.mk


test_ndgrid: test_ndgrid.F90
	$(F90C) $(LDFLAGS) -I.. -o $@ test_ndgrid.F90 ../ndgrid.o ../matoper.o ../ufileformat.o $(LIBS)

print:
	echo $(LIBS) $(F90FLAGS) $(PRECISION)
