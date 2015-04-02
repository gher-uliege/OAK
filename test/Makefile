include ../config.mk
include ../Compilers/$(OS)-$(strip $(FORT)).mk
include ../Compilers/libs.mk


test_ndgrid: test_ndgrid.F90 ../ndgrid.F90 ../ndgrid_inc.F90 ../matoper.F90 ../ufileformat.F90
	$(F90C) $(F90FLAGS) $(LDFLAGS) -I.. -o $@ test_ndgrid.F90 ../ndgrid.F90 ../matoper.F90 ../ufileformat.F90 $(LIBS)

print:
	echo $(LIBS) $(F90FLAGS) $(PRECISION)

clean:
	rm -f *.mod *.o
