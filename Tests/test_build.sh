#!/bin/bash

# it is an useful script to automatically build and run most cases in T-Flows.

# put your compilers here
FCOMP="gnu"; # or ifort/gfortran/mpif90/mpifort/mpiifort
DEBUG="yes"; # run tests in debug mode

# folder structure
TEST_DIR=$PWD                      # dir with tests
GENE_DIR=$PWD/../Sources/Generate  # Generate src folder
CONV_DIR=$PWD/../Sources/Convert   # Convert  src folder
DIVI_DIR=$PWD/../Sources/Divide    # Divide   src folder
PROC_DIR=$PWD/../Sources/Process   # Process  src folder
BINA_DIR=$PWD/../Binaries/         # binaries folder

# executables
GENE_EXE=$BINA_DIR/Generate        # Generate ex
CONV_EXE=$BINA_DIR/Convert         # Convert  ex
DIVI_EXE=$BINA_DIR/Divide          # Divide   ex
PROC_EXE=$BINA_DIR/Process         # Process  ex

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
# make links
#------------------------------------------------------------------------------#
function make_links {
  ln -rsf $BINA_DIR/* $TEST_DIR/Laminar/Backstep_Orthogonal/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Laminar/Backstep_Nonorthogonal/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Backstep_Re_26000_Rsm/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Backstep_Re_28000/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Channel_Re_Tau_590/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Fuel_Bundle/;
}
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
function convert_tests {
  # unpacking geometry
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                    tar -zxvf chan.tar.gz
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;      tar -zxvf jet.tar.gz
  cd $TEST_DIR/Rans/Fuel_Bundle;           tar -zxvf subflow_LowRe_coarse.tar.gz

  #-- seq, no cgns
  cd $CONV_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                  $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                         $CONV_EXE < convert.scr

  #-- seq, cgns (adf5)
  cd $CONV_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_ADF5=yes
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                  $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                         $CONV_EXE < convert.scr

  #-- seq, cgns (hdf5)
  cd $CONV_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=yes
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                  $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                         $CONV_EXE < convert.scr
}
#------------------------------------------------------------------------------#
# Divide tests
#------------------------------------------------------------------------------#
function divide_tests {
  #-- seq
  cd $DIVI_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;               $DIVI_EXE < divide.scr
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;            $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;                $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Backstep_Re_28000;                    $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;     $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                   $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;    $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                          $DIVI_EXE < divide.scr
}
#------------------------------------------------------------------------------#
# processor tests
#------------------------------------------------------------------------------#
function processor_tests {
  #-- seq, no cgns
  cd $PROC_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG MPI=no
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;
  cd $TEST_DIR/Rans/Backstep_Re_28000;
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;

  #-- seq, cgns(adf5)
  cd $PROC_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_ADF5=yes MPI=no
  cd $TEST_DIR/Laminar/Backstep_Orthogonal
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;
  cd $TEST_DIR/Rans/Backstep_Re_28000;
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;

  #-- seq, cgns(hdf5)
  cd $PROC_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=yes MPI=no
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;
  cd $TEST_DIR/Rans/Backstep_Re_28000;
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;

  #-- par, no cgns
  cd $PROC_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=no MPI=yes
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;
  cd $TEST_DIR/Rans/Backstep_Re_28000;
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;

  #-- par, cgns(hdf5)
  cd $PROC_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=yes MPI=yes
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;
  cd $TEST_DIR/Rans/Backstep_Re_28000;
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;
}
#------------------------------------------------------------------------------#
# actual script
#------------------------------------------------------------------------------#
#generator_tests
make_links
convert_tests
divide_tests
processor_tests
