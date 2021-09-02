Installation
------------

To install OAK you need:

* A Fortran 90 compiler
* NetCDF (with Fortran 90 interface)
* LAPACK library
* an implementation of the BLAS library (such as ATLAS, OpenBLAS or MKL) 
* cholmod

Optionally, for parallelization, either one of:

* MPI (with Fortran 90 interface)
* OpenMP (included in Fortran 90 compilers) 


# Compilation

## gfortran and default BLAS

To compile with open-source gfortran, parallelized with MPI, at least the following packages are needed:

```
sudo apt-get install gfortran libatlas-base-dev liblapack-dev openmpi-bin libopenmpi-dev libnetcdf-dev netcdf-bin
cp config.mk.template config.mk
make
```

# ifort and MKL

## Without CLODMOD:

Set `$MKLROOT`:

```bash
$ make FORT=ifort NETCDF_CONFIG=/path/to/bin/nf-config CHOLMOD_LIB= BLAS_LIBDIR=$MKLROOT/lib/intel64/ BLAS_LIB="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core" LAPACK_LIB=
```
