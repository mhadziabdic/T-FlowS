#!/bin/bash

# it is an useful script to automatically build and run most cases in T-Flows.

# put your compilers here
FCOMP="gnu"; # or ifort/gfortran/mpif90/mpifort/mpiifort
DEBUG="yes"; # run tests in debug mode

# folder structure
TEST_DIR=$PWD                      # dir with tests
GENE_DIR=$PWD/../Sources/Generate  # generator src folder
CONV_DIR=$PWD/../Sources/Convert   # generator src folder
DIVI_DIR=$PWD/../Sources/Divide    # divisor src folder
PROC_DIR=$PWD/../Sources/Processor # processor src folder
BINA_DIR=$PWD/../Binaries/         # binaries folder

# executables
GENE_EXE=$BINA_DIR/Generate
CONV_EXE=$BINA_DIR/Convert
DIVI_EXE=$BINA_DIR/Divisor
PROC_EXE=$BINA_DIR/Processor
#CGNS_DIR=/home/l_palkin_e/Development/T-Flows-modern/Sources/Libraries/install_dir/cgnslib_latest_linux_64_hdf5_seq/bin/

#-------------------------------------------------#
#---------   READ ABOVE UP TO THIS ROW   ---------#
#-------------------------------------------------#

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


#------------------------------------------------------------------------------#
# generator tests
#------------------------------------------------------------------------------#
function generator_tests {
  #-- seq, no cgns
  cd $GENE_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;             $GENE_EXE < generate.scr
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;          $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;              $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_28000;                  $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;   $GENE_EXE < generate.scr

  #-- seq, cgns(adf5)
  cd $GENE_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_ADF5=yes
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;             $GENE_EXE < generate.scr
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;          $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;              $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_28000;                  $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;   $GENE_EXE < generate.scr

  #-- seq, cgns(hdf5)
  cd $GENE_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=yes
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;             $GENE_EXE < generate.scr
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;          $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;              $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_28000;                  $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;   $GENE_EXE < generate.scr
}
#------------------------------------------------------------------------------#
# convert tests
#------------------------------------------------------------------------------#
#-- seq, no cgns
cd $CONV_DIR; make clean
make FORTRAN=$FCOMP DEBUG=$DEBUG
cd $TEST_DIR/Rans/Channel_Re_Tau_590;                    $CONV_EXE < convert.scr

#-- seq, cgns (adf5)
cd $CONV_DIR; make clean
make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_ADF5=yes
cd $TEST_DIR/Rans/Channel_Re_Tau_590;                    $CONV_EXE < convert.scr

#-- seq, cgns (hdf5)
cd $CONV_DIR; make clean
make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=yes
cd $TEST_DIR/Rans/Channel_Re_Tau_590;                    $CONV_EXE < convert.scr

#------------------------------------------------------------------------------#
# divisor tests
#------------------------------------------------------------------------------#
cd $DIVI_DIR; make clean
make FORTRAN=$FCOMP DEBUG=$DEBUG
