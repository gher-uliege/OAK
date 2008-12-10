#!/bin/bash -e

export FC=gfortran
export CC=gcc
export CXX=g++
export INSTALLROOT=/opt/gcc-4.1.2


export FC=/home/abarth/Download/gcc-4.2/bin/gfortran
export CC=gcc
export CXX=g++
export INSTALLROOT=/opt/gcc-4.2.1-soft

#export FC=/opt/g95-0.9/bin/x86_64-unknown-linux-gnu-g95
#export CC=gcc
#export CXX=g++
#export INSTALLROOT=/opt/g95-0.9-soft/

#export FC=pgf90
#export CC=gcc
#export CXX=g++
#export INSTALLROOT=/opt/pgi-6.0


BUILDDIR=$(pwd)
JOBS=6

# NETCDF

rm -Rf netcdf
cd $BUILDDIR
wget -O - ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf.tar.gz  | tar xvzf -
DIR=$(ls -d netcdf-* | tail -n1)
cd $DIR
./configure --prefix=$INSTALLROOT/$DIR
perl -pi -e 's/f2cFortran/gFortran/' config.h
perl -pi -e 's/undef NF_INT_IS_C_INT/define NF_INT_IS_C_INT 1/' nfconfig.inc
make -j $JOBS
make check
make install

# LAM MPI

rm -Rf lam-7.1.3
cd $BUILDDIR
wget -O - http://www.lam-mpi.org/download/files/lam-7.1.3.tar.bz2 | tar jxf -
DIR=$(ls -d lam-* | tail -n1)
cd $DIR
./configure --prefix=$INSTALLROOT/$DIR
make -j $JOBS
make install

# OPeNDAP

URL=http://www.opendap.org/pub/source/libdap-3.7.8.tar.gz
NAME=$(basename $URL | sed 's/^\(.*\)-[^-]*$/\1/')
cd $BUILDDIR
wget -O - $URL | tar zxf -
DIR=$(ls -d $NAME-* | tail -n1)
cd $DIR
./configure --prefix=$INSTALLROOT/$DIR
make -j $JOBS
make install
export PATH=$INSTALLROOT/libdap-*/bin:$PATH

# NCDAP

URL=ftp://ftp.unidata.ucar.edu/pub/opendap/source/libnc-dap-3.7.0.tar.gz
NAME=$(basename $URL | sed 's/^\(.*\)-[^-]*$/\1/')
cd $BUILDDIR
wget -O - $URL | tar zxf -
DIR=$(ls -d $NAME-* | tail -n1)
cd $DIR
./configure F77=$FC --prefix=$INSTALLROOT/$DIR
make -j $JOBS
make install
mv $INSTALLROOT/$DIR/bin/ncdump $INSTALLROOT/$DIR/bin/dncdump

# LAPACK and BLAS

cd $BUILDDIR
wget -O - http://www.netlib.org/lapack/lapack-3.1.1.tgz | tar zvxf -
DIR=$(ls -d lapack-* | tail -n1)
cd $DIR
cp make.inc.example  make.inc
mkdir -p $INSTALLROOT/$DIR/lib
make FORTRAN=$FC LOADER=$FC PLAT= BLASLIB=libblas.a LAPACKLIB=liblapack.a blaslib lapacklib 
cp BLAS/SRC/libblas.a liblapack.a $INSTALLROOT/$DIR/lib/

# ARPACK

cd $BUILDDIR
wget -O - http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz | tar xzvf -
wget -O - http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz | tar xzvf -
cd ARPACK
mkdir -p $INSTALLROOT/arpack96/lib
make FC=$FC FFLAGS= MAKE=/usr/bin/make ARPACKLIB=$INSTALLROOT/arpack96/lib/libarpack.a home=$BUILDDIR/ARPACK lib

# make symlinks

rm -Rf $INSTALLROOT/bin
rm -Rf $INSTALLROOT/lib
rm -Rf $INSTALLROOT/include

mkdir $INSTALLROOT/bin
mkdir $INSTALLROOT/lib
mkdir $INSTALLROOT/include

for d in bin lib include; do
    for i in $INSTALLROOT/*/$d/*; do
	t=$INSTALLROOT/$d/$(basename $i)
	if [ -e $t ]; then
	    echo "File $t already exists; skipped;"
	else
	    ln -s $i $INSTALLROOT/$d
	fi
    done
done


#ln -s $INSTALLROOT/*/bin/* $INSTALLROOT/bin
#ln -s $INSTALLROOT/*/lib/* $INSTALLROOT/lib
#ln -s $INSTALLROOT/*/include/* $INSTALLROOT/include