#===============================================================================
# 
#   Process Makefile   
#   
#-------------------------------------------------------------------------------

#--------------------------
#   Variable definitions
#--------------------------

# Fortran compiler ("gnu", "intel" or "portland")
FORTRAN = gnu

# Debugging ("yes" or "no")
DEBUG = no

# Compile for MPI ("yes" or "no")
MPI = no

# Link with CGNS ("yes" or "no")
CGNS_ADF5 = no
CGNS_HDF5 = no

# Directories for objects and modules. (No need to change.)
DIR_BINARY = ../../Binaries
DIR_SHARED = ../Shared
DIR_MODULE = ../.Modules
DIR_OBJECT = .Objects
DIR_MODELS = Turbulence
DIR_COMMUN = Communication/Sequential
ifeq ($(MPI), yes)
  DIR_COMMUN = Communication/Parallel/Double
endif
VPATH = $(DIR_SHARED):$(DIR_MODELS):$(DIR_COMMUN):.

# Program name (This should hardly change)
PROGRAM_NAME = Process
PROGRAM_FILE = $(DIR_BINARY)/$(PROGRAM_NAME)

#-------------------------------------------------------------------------
#   Compiler and linker options
#-------------------------------------------------------------------------
#   Note: Changes only when support to a new Fortran compiler is added.
#-------------------------------------------------------------------------

# Fortran == gnu
ifeq ($(FORTRAN), gnu)
  FC = gfortran
  ifeq ($(DEBUG),yes)
    OPT_COMP = -J $(DIR_MODULE) -fdefault-real-8 -fdefault-integer-8 -O0 -g -Wall -fcheck=all
  else
    OPT_COMP = -J $(DIR_MODULE) -fdefault-real-8 -fdefault-integer-8 -O3
  endif 
  OPT_LINK = $(OPT_COMP)
endif 

# Fortran == intel
ifeq ($(FORTRAN), intel)
  FC = ifort
  ifeq ($(DEBUG),yes)
    OPT_COMP = -module $(DIR_MODULE) -r8 -O0 -warn all -check all -debug all
  else
    OPT_COMP = -module $(DIR_MODULE) -r8 -O3
  endif
  OPT_LINK = $(OPT_COMP)
endif 

# Fortran == portland
ifeq ($(FORTRAN), portland)
  FC = pgfortran
  ifeq ($(DEBUG),yes)
    OPT_COMP = -module $(DIR_MODULE) -r8 -O0
  else
    OPT_COMP = -module $(DIR_MODULE) -r8 -O3
  endif
  OPT_LINK = $(OPT_COMP)
endif 

# It seems that for all supported compilers ("gnu", "intel" 
# or "portland") command to compile in parallel is "mpif90"
ifeq ($(MPI), yes)
  FC = mpif90
endif

# If compiling with CGNS, add path to include directory
# (One day, when you want to work with compressed 
#  files, the following line might still be needed:
#  CGNS_LIB   += -Wl,--no-as-needed -ldl -lz
ifeq ($(CGNS_ADF5), yes)
  DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Adf5_Seq/Include
  CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Adf5_Seq/Lib/libcgns.a
  OPT_COMP   += $(addprefix -I, $(DIR_CGNS_H))
endif

ifeq ($(CGNS_HDF5), yes)
  DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par/Include
  CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par/Lib/libcgns.a \
                ../Libraries/Hdf5_Par/Lib/libhdf5_fortran.a                        \
                -Wl,--no-as-needed -ldl -lz

  #DIR_CGNS_H  = /home/l_palkin_e/eclipse/T-FlowS-Egor/Library/CGNS/install_dir/cgnslib_3.2.1_linux_64_hdf5_par/include
  #CGNS_LIB    = /home/l_palkin_e/eclipse/T-FlowS-Egor/Library/CGNS/install_dir/cgnslib_3.2.1_linux_64_hdf5_par/libcgns.a \
  #              /home/l_palkin_e/eclipse/T-FlowS-Egor/Library/CGNS/install_dir/HDF5_par/lib/libhdf5.a                    \
  #              -Wl,--no-as-needed -ldl -lz
  OPT_COMP   += $(addprefix -I, $(DIR_CGNS_H))
endif

$(info #=================================================================)
$(info # Compiling $(PROGRAM_NAME) with compiler $(FORTRAN)              )
$(info #-----------------------------------------------------------------)
$(info # Usage:                                                          )
$(info #   make <FORTRAN=gnu/intel/portland> <DEBUG=yes/no> <MPI=yes/no> )
$(info #        <CGNS_ADF=yes/no> <CGNS_HDF5=yes/no>                     )
$(info #-----------------------------------------------------------------)

#------------------------------------------------------
#   List of sources for modules and functions
#------------------------------------------------------
#   Modules' order must obey their dependency 
#   This list should therefore be written "by hand".
#   Note: Modules written in lower case 
#         letters are candidates for deletion.
#------------------------------------------------------

# Sources for modules in shared directory
SRC_MOD =  allp_mod.f90			\
           all_mod.f90			\
           par_mod.f90			\
           Grid_Mod.f90			\
           Work_Mod.f90			\
           Tokenizer_Mod.f90    \
           Material_Mod.f90

# Sources for modules in the local directory
SRC_MOD += Constants_Pro_Mod.f90\
           Info_Mod.f90 		\
           Bulk_Mod.f90 		\
           Var_Mod.f90			\
           Matrix_Mod.f90		\
           Solvers_Mod.f90		\
           pro_mod.f90			\
           User_Mod.f90

# Sources for modules in sub-directories
SRC_MOD += Turbulence_Mod.f90	\
           les_mod.f90			\
           rans_mod.f90			\

# Sources for modules in the local directory
SRC_MOD += Cgns_Mod.f90

# If CGNS is supported, add the right module
ifeq ($(CGNS_ADF5), $(CGNS_HDF5))  # meaning both are "no"
  SRC_MOD += Cgns_Mod_Void.f90
endif

# Sources for all functions are obtained by a shell command
SRC_RAW =  $(shell ls -1 *.f90               | grep -v -i _Mod)
SRC_RAW += $(shell ls -1 $(DIR_MODELS)/*.f90 | grep -v -i _Mod)
SRC_RAW += $(shell ls -1 $(DIR_COMMUN)/*.f90 | grep -v -i _Mod)

# Remove directory names
STEP_A  = $(subst $(DIR_MODELS),,$(SRC_RAW)) # remove directory name
STEP_B  = $(subst $(DIR_COMMUN),,$(STEP_A))  # remove directory name
SRC_FUN = $(subst /,,$(STEP_B))              # remove last slash 

# Add files from "Shared" directory, but one by one, to avoid clutter 
SRC_FUN += Distance.f90			 \
           Load_Geo.f90			 \
           Name_File.f90		 \
           Sort_Int_Carry_Int.f90\
           Sort_Real_Carry_Int.f90\
           Approx.f90\
           Swap_Int.f90\
           To_Upper_Case.f90

# ... to allow conditional inclusion of CGNS support
#ifeq ($(CGNS_ADF5), $(CGNS_HDF5))
#  SRC_FUN += Save_Cgns_Results_Void.f90
#else
#  SRC_FUN += Save_Cgns_Results.f90
#endif

#----------------------------------------------------------------------
#   List of objects generated from the list of modules and functions  
#----------------------------------------------------------------------
#   Note: This doesn't need editing.
#----------------------------------------------------------------------
OBJ_MOD = $(SRC_MOD:%.f90=$(DIR_OBJECT)/%.o)
OBJ_FUN = $(SRC_FUN:%.f90=$(DIR_OBJECT)/%.o)
OBJ = $(OBJ_MOD) $(OBJ_FUN)

#-------------------------------------------------------
#   List of modules currently used for target "clean" 
#-------------------------------------------------------
#   Note: This doesn't need editing.
#-------------------------------------------------------
SRC_MOD_LOW = $(shell echo $(SRC_MOD) | tr A-Z a-z)
MOD = $(SRC_MOD_LOW:%.f90=$(DIR_MODULE)/%.mod)

#---------------------------------------------------------
#   Default rule to build Fortran modules and functions
#---------------------------------------------------------
#   Note: This doesn't need editing.
#---------------------------------------------------------

# Modules
$(DIR_OBJECT)/%.o: %.f90 %/*.f90
	$(FC) $(OPT_COMP) -c -o $@ $<

# Functions
$(DIR_OBJECT)/%.o: %.f90
	$(FC) $(OPT_COMP) -c -o $@ $<

#-----------------------------------
#   Rule to build main program
#-----------------------------------
#   Note: Should not be modified.
#-----------------------------------
$(PROGRAM_FILE): $(OBJ) 
	$(FC) $(OPT_LINK) -o $(PROGRAM_FILE) $(OBJ) $(CGNS_LIB)

#---------------------------------------
#   Explicit dependencies for modules
#---------------------------------------
#   These should be inserted by 
#   hand for tuning of dependencies.
#---------------------------------------
$(DIR_OBJECT)/all_mod.o:	$(DIR_OBJECT)/allp_mod.o
$(DIR_OBJECT)/Grid_Mod.o:	$(DIR_OBJECT)/Bnd_Cond_Mod.o\
				$(DIR_OBJECT)/Material_Mod.o
$(DIR_OBJECT)/Domain_Mod.o:	$(DIR_OBJECT)/Point_Mod.o	\
  				$(DIR_OBJECT)/Line_Mod.o	            \
  				$(DIR_OBJECT)/Block_Mod.o	            \
  				$(DIR_OBJECT)/Region_Mod.o
$(DIR_OBJECT)/Cgns_Mod.o:	$(DIR_OBJECT)/Grid_Mod.o

#---------------------
#   Explicit target.
#---------------------
clean:
	rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod $(PROGRAM_FILE)
