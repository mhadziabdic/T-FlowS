#===============================================================================================================
#
#   Process Makefile
#
#---------------------------------------------------------------------------------------------------------------

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
DIR_MODULE = .Modules
DIR_OBJECT = .Objects
DIR_MODELS = Turbulence
DIR_USER   = User_Mod
DIR_CASE   = no
VPATH = $(DIR_SHARED):$(DIR_MODELS):.

# Program name (This should hardly change)
PROGRAM_NAME = Process
PROGRAM_FILE = $(DIR_BINARY)/$(PROGRAM_NAME)

$(info #=======================================================================================================)
$(info # Compiling $(PROGRAM_NAME) with compiler $(FORTRAN)                                                    )
$(info #-------------------------------------------------------------------------------------------------------)
$(info # Usage:                                                                                                )
$(info #   make <FORTRAN=gnu/intel/portland> <DEBUG=yes/no> <MPI=yes/no>  <DIR_CASE=relative_path_to_case> ... )
$(info #        ... <CGNS_ADF5=yes/no> <CGNS_HDF5=yes/no>                                                      )   
$(info #                                                                                                       )
$(info # Examples:                                                                                             )
$(info #   make                            for sequential version with gnu compiler                            )
$(info #   make MPI=yes                    for parallel version with gnu compiler                              )
$(info #   make FORTRAN=intel MPI=yes      for parallel version with intel compiler                            )
$(info #   make MPI=yes CGNS_ADF5=yes      for parallel version with CGNS/ADF5 support                         )
$(info #-------------------------------------------------------------------------------------------------------)

#-------------------------------------------------------------------------
#   Compiler and linker options
#-------------------------------------------------------------------------
#   Note: Changes only when support to a new Fortran compiler is added.
#-------------------------------------------------------------------------

# Fortran == gnu
ifeq ($(FORTRAN), gnu)
  FC = gfortran
  ifeq ($(DEBUG),yes)
    OPT_COMP = -J $(DIR_MODULE) -fdefault-real-8 -fdefault-integer-8 -O0 -g -Wall -fcheck=all -ffixed-line-length-80
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
ifeq ($(CGNS_ADF5), yes)
  DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Adf5_Seq/Include
  CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Adf5_Seq/Lib/libcgns.a
endif

ifeq ($(CGNS_HDF5), yes)
  ifeq ($(MPI), yes)
    DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par/Include
    CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Par/Lib/libcgns.a	\
                  ../Libraries/Hdf5_Par/Lib/libhdf5.a				\
                  -Wl,--no-as-needed -ldl -lz
  else
    DIR_CGNS_H  = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Include
    CGNS_LIB    = ../Libraries/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Lib/libcgns.a	\
                  ../Libraries/Hdf5_Seq/Lib/libhdf5.a				\
                  -Wl,--no-as-needed -ldl -lz
  endif
endif
OPT_COMP += $(addprefix -I, $(DIR_CGNS_H))

#------------------------------------------------------
#   List of sources for modules and functions
#------------------------------------------------------
#   Modules' order must obey their dependency
#   This list should therefore be written "by hand".
#   Note: Modules written in lower case
#         letters are candidates for deletion.
#------------------------------------------------------

#--------------------
#   User functions
#--------------------
ifneq ($(DIR_CASE), no)
  ifeq ('$(wildcard $(DIR_CASE))','')
    $(info Case directory $(DIR_CASE) does not exist!)
  else
    $(info Compling for the case in $(DIR_CASE))
    ifeq ('$(wildcard $(DIR_CASE)/$(DIR_USER))','')
      $(info Directory with user functions $(DIR_CASE)/$(DIR_USER) does not exist!)
    else
      $(info Linking user functions from directory $(DIR_CASE)/$(DIR_USER))

      # Fetch all user functions from the case directory
      STEP_A  = $(shell ls -1 $(DIR_CASE)/$(DIR_USER)/*.f90)

      # Append "../" to each of the file names
      SRC_USR = $(patsubst %.f90, ../%.f90, $(STEP_A))

      # Make soft links
      $(shell ln -f -s $(SRC_USR) $(DIR_USER))
    endif
  endif
endif

#-------------
#   Modules
#-------------

# Modules in shared directory
SRC_MOD = Const_Mod.f90		\
          Grid_Mod.f90		\
          Work_Mod.f90		\
          Name_Mod.f90		\
          Tokenizer_Mod.f90	\
          Material_Mod.f90

# Conditional support for MPI communication
ifeq ($(MPI), yes) 
  SRC_MOD += Comm_Mod_Par.f90
else
  SRC_MOD += Comm_Mod_Seq.f90
endif

# Sources for modules in the local directory
SRC_MOD += Info_Mod.f90		\
           Bulk_Mod.f90		\
           Var_Mod.f90		\
           Matrix_Mod.f90	\
           Numerics_Mod.f90	\
           Solvers_Mod.f90	\
           Flow_Mod.f90		\
           Grad_Mod.f90

# Sources for modules in sub-directories
SRC_MOD += Turbulence_Mod.f90	\
           les_mod.f90		\
           rans_mod.f90		\

# User_Mod should be here
SRC_MOD += User_Mod.f90

# Control_Mod should be last since it reads for other modules
SRC_MOD += Control_Mod.f90

#---------------
#   Functions   
#---------------

# Sources for all functions are obtained by a shell command
SRC_RAW =  $(shell ls -1 *.f90 | grep -v -i _Mod | grep -v -i Save_Cgns)
SRC_RAW += $(shell ls -1 $(DIR_MODELS)/*.f90 | grep -v -i _Mod)

# Remove directory names
STEP_A  = $(subst $(DIR_MODELS),,$(SRC_RAW)) # remove directory name
SRC_FUN = $(subst /,,$(STEP_A))              # remove last slash

# Add files from "Shared" directory, but one by one, to avoid clutter
SRC_FUN += Distance.f90 		\
           Load_Geo.f90 		\
           Name_File.f90		\
           Sort_Int_Carry_Int.f90	\
           Sort_Real_Carry_Int.f90	\
           Approx.f90	 		\
           Swap_Int.f90			\
           To_Upper_Case.f90

#-------------------------------------------
#   Conditional inclusion of CGNS support 
#-------------------------------------------
ifeq ($(CGNS_ADF5), $(CGNS_HDF5))  # meaning both are "no"
  SRC_MOD += Cgns_Mod_Void.f90
  SRC_FUN += Save_Cgns_Results_Void.f90	\
             Save_Cgns_Cells_Void.f90
else
  ifeq ($(MPI), yes)
    SRC_MOD += Cgns_Mod_Par.f90
  else
    SRC_MOD += Cgns_Mod_Seq.f90
  endif
  SRC_FUN += Save_Cgns_Results.f90	\
             Save_Cgns_Cells.f90
endif

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
$(DIR_OBJECT)/all_mod.o:    $(DIR_OBJECT)/allp_mod.o
$(DIR_OBJECT)/Grid_Mod.o:   $(DIR_OBJECT)/Bnd_Cond_Mod.o	\
                            $(DIR_OBJECT)/Material_Mod.o
$(DIR_OBJECT)/Domain_Mod.o: $(DIR_OBJECT)/Point_Mod.o		\
                            $(DIR_OBJECT)/Line_Mod.o		\
                            $(DIR_OBJECT)/Block_Mod.o		\
                            $(DIR_OBJECT)/Region_Mod.o
$(DIR_OBJECT)/Cgns_Mod.o:   $(DIR_OBJECT)/Grid_Mod.o

#---------------------
#   Explicit target.
#---------------------
clean:
	rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod $(PROGRAM_FILE)
	rm -f $(DIR_USER)/*.f90
	git checkout $(DIR_USER)
