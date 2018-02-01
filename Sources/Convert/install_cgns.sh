#!/bin/bash

CGNS_DIR=$PWD # this is top dir for this lib
INSTALL_DIR=$CGNS_DIR/install_dir # this dir contain CNGS/ HDF5/ and optional folders
SRC_DIR=$CGNS_DIR/src_dir # dir with sources

# decide if you wish to build mpich by yourself
BUILD_MPI=true
# decide if you wish to build gui cgns tools (cgnsview)
CGNS_TOOLS=true

# put your compilers here (gcc, gfortran are allowed if $MPI is built anyway)
export CC="gcc";
export FC="gfortran"; # or ifort, mpif90, mpifort, mpiifiort

mkdir -p $SRC_DIR/
mkdir -p $INSTALL_DIR/

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


#--------- functions definitions

function build_mpi_lib {
#------ MPICH 3.2.1


	# download sources
	cd $SRC_DIR/
	wget -N http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz
	tar -zxvf mpich-3.2.1.tar.gz
	mkdir -p $SRC_DIR/MPICH/
	rsync -azvh mpich-3.2.1/* $SRC_DIR/MPICH/
	rm -r mpich-3.2.1/
		
	# configure
	cd $SRC_DIR/MPICH/
	./configure \
	--prefix=$INSTALL_DIR/MPICH \
	--enable-fast=all,O3 \
	--enable-fortran=all \
	--disable-shared \
	--disable-dependency-tracking

	# build
	make
	# install
	make install

	# return
	cd $CGNS_DIR
}

#------ HDF5 5.1.8 (Paraview 5.4.1 & Visit 2.12.3 work with HDF5 5.1.8 and not with 5.1.10)
function build_hdf5_lib {

	# download sources
	cd $SRC_DIR/
	mkdir -p $SRC_DIR/HDF5/
	git clone --depth=1 https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git --branch hdf5_1_8 ./hdf5
	rsync -azvh hdf5/* $SRC_DIR/HDF5/
	rm -rf hdf5

	# configure
	cd HDF5/; rm -rf .git
	FCFLAGS=-O3 \
	./configure \
	--prefix=$INSTALL_DIR/HDF5 \
	--enable-fortran \
	--enable-parallel \
	--disable-shared \
	--enable-production

	# build
	make
	# install
	make install

	cd $CGNS_DIR
}

#------ TCL
function build_tcl_lib {

	cd $SRC_DIR/

	# download sources
	cd $SRC_DIR/
	wget -N https://prdownloads.sourceforge.net/tcl/tcl8.6.8-src.tar.gz
	tar -zxvf tcl8.6.8-src.tar.gz;
	mkdir -p $SRC_DIR/TCL/
	rsync -azvh tcl8.6.8/* $SRC_DIR/TCL/
	rm -r tcl8.6.8/

	# configure
	cd TCL/unix/
	./configure \
	--prefix=$INSTALL_DIR/TCL

	# build
	make
	# install
	make install

	cd $CGNS_DIR
}

#------ TK (requires libx11-dev from repo)
function build_tk_lib {

	cd $SRC_DIR/

	# download sources
	cd $SRC_DIR/
	wget -N https://prdownloads.sourceforge.net/tcl/tk8.6.8-src.tar.gz
	tar -zxvf tk8.6.8-src.tar.gz;
	mkdir -p $SRC_DIR/TCL/
	rsync -azvh tk8.6.8/* $SRC_DIR/TK/
	rm -r tk8.6.8/

	# configure
	cd TK/unix/
	./configure \
	--prefix=$INSTALL_DIR/TK \
	--with-tcl=$INSTALL_DIR/TCL/lib/

	# build
	make
	# install
	make install

	cd $CGNS_DIR

}

function build_cgns_lib {
#------ parallel CGNS with HDF5 (latest version; requires zlib)

	# download sources
	cd $SRC_DIR/
	mkdir -p $SRC_DIR/CGNS/
	git clone --depth=1 https://github.com/CGNS/CGNS.git ./cgns
	rsync -azvh cgns/* $SRC_DIR/CGNS/
	rm -rf cgns

if [ $CGNS_TOOLS == true ]; then

	# make links to TCL and TK where cgns searches for them
	cd $INSTALL_DIR
	rm -rf TCL/unix TCL/generic TK/unix TK/generic

	ln -s -r -f TCL/lib     TCL/tmp; mv TCL/tmp TCL/unix
	ln -s -r -f TCL/include TCL/tmp; mv TCL/tmp TCL/generic

	ln -s -r -f TK/lib      TK/tmp;  mv TK/tmp  TK/unix
	ln -s -r -f TK/include  TK/tmp;  mv TK/tmp  TK/generic

	# configure
	cd $SRC_DIR/CGNS/; rm -rf .git; cd src/

	FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	./configure \
	--prefix=$INSTALL_DIR/CGNS_HDF5 \
	--with-hdf5=$INSTALL_DIR/HDF5 \
	--with-fortran \
	--enable-lfs \
	--enable-64bit \
	--disable-shared \
	--disable-debug \
	--with-zlib \
	--enable-parallel \
	--enable-cgnstools \
	--with-tcl=$INSTALL_DIR/TCL \
	--with-tk=$INSTALL_DIR/TK \
	--datarootdir=$INSTALL_DIR/CGNS/tcl_scripts

else # no cgns GUI tools
	
	# configure
	cd CGNS/; rm -rf .git; cd src/

	FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	./configure \
	--prefix=$CGNS_DIR/CGNS/ \
	--with-hdf5=$CGNS_DIR/HDF5/ \
	--with-fortran \
	--enable-lfs \
	--enable-64bit \
	--disable-shared \
	--disable-debug \
	--with-zlib \
	--enable-parallel \
	--disable-cgnstools

fi

# build
make
# install
make install

# build and run tests for self-confidence

# cd tests
# make
# make test
# cd ../examples/fortran
# make
# make test
# cd ../../Test_UserGuideCode/Fortran_code
# make
# make test
# cd ../C_code
# make
# make test

#------ sequential CGNS without HDF5 (->ADF5)

FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
./configure \
--prefix=$INSTALL_DIR/CGNS_ADF5 \
--with-fortran \
--enable-lfs \
--enable-64bit \
--disable-shared \
--disable-debug \
--with-zlib \

# build
make
# install
make install

cd $CGNS_DIR

}

#--------- script with functions defined above

if [ $BUILD_MPI == true ]; then

	build_mpi_lib

	# mpicc and mpif90 compilers from now are:
	export CC=$INSTALL_DIR/MPICH/bin/mpicc
	export FC=$INSTALL_DIR/MPICH/bin/mpif90
fi

build_hdf5_lib

if [ $CGNS_TOOLS == true ]; then
	build_tcl_lib
	build_tk_lib
fi

build_cgns_lib

echo \n
echo parallel   CNGS with HDF5         is now installed in $INSTALL_DIR/CGNS_HDF5
echo sequential CNGS with without HDF5 is now installed in $INSTALL_DIR/CGNS_ADF5
echo CNGS tools are installed in $INSTALL_DIR/CGNS_HDF5/bin/
echo You can make relative links to them for convinience
echo \n
echo \n
echo you can safely remove $SRC_DIR/ folder with its content