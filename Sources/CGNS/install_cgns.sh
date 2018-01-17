#!/bin/bash

# current dir which will contain CNGS/ HDF5/ and optionally MPICH/ folders
CGNS_DIR=$PWD

# decide if you wish to build mpich by yourself
BUILD_MPI=true

# put your compilers here (gcc, gfortran are allowed if $MPI is built anyway)
export CC="gcc";
export FC="gfortran"; # or ifort, mpif90, mpifort, mpiifiort

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


#------ 1 : MPICH 3.2.1
if [ $BUILD_MPI == true ]; then

	# download
	wget -N http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz; tar -zxvf mpich-3.2.1.tar.gz; mv mpich-3.2.1/ MPICH/;
		
	# configure
	cd MPICH/;
	./configure --prefix=$PWD/install_dir/ --enable-fast=all,O3 --enable-fortran=all --disable-shared --disable-dependency-tracking

	# build
	make

	# install
	make install

	# and now your mpicc and mpif90 compilers are:
	export CC=$CGNS_DIR/MPICH/install_dir/bin/mpicc
	export FC=$CGNS_DIR/MPICH/install_dir/bin/mpif90

	# return
	cd $CGNS_DIR
	rm mpich-3.2.1.tar.gz
fi

#------ 2 : HDF5 5.1.8 (Paraview 5.4.1 & Visit 2.12.3 work with HDF5 5.1.8 and not with 5.1.10)

# download
git clone --depth=1  https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git --branch hdf5_1_8 --single-branch ./HDF5

# configure
cd HDF5/; 
FCFLAGS=-O3 ./configure --prefix=$PWD/install_dir/ --enable-fortran  --enable-parallel --disable-shared

# build
make

# install
make install

cd $CGNS_DIR

#------ 3 : CGNS (latest version)

# download
git clone --depth=1  https://github.com/CGNS/CGNS.git CGNS/; cd CGNS/; rm -rf .git; cd src/

# configure
FLIBS=-Wl,--no-as-needed\ -ldl\ -lz  LIBS=-Wl,--no-as-needed\ -ldl\ -lz ./configure --prefix=$CGNS_DIR/CGNS/install_dir --with-hdf5=$CGNS_DIR/HDF5/install_dir --with-fortran --enable-lfs --enable-64bit --disable-shared --disable-debug --with-zlib --disable-cgnstools --enable-64bit --enable-parallel

# build
make
# install
make install

# run tests
cd tests                                   # for self-confidence
make                                       # for self-confidence
make test                                  # for self-confidence
cd ../examples/fortran                     # for self-confidence
make                                       # for self-confidence
make test                                  # for self-confidence
cd ../../Test_UserGuideCode/Fortran_code   # for self-confidence
make                                       # for self-confidence
make test                                  # for self-confidence
cd ../C_code                               # for self-confidence
make                                       # for self-confidence
make test                                  # for self-confidence


cd $CGNS_DIR

#---------

echo CNGS is now built in "$CGNS_DIR/CGNS/install_dir".