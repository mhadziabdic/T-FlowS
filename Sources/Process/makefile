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

# Directories for objects and modules. (No need to change.)
DIR_BINARY = ../../Binaries
DIR_MODULE = ../.Modules
DIR_OBJECT = ../.Objects
DIR_SHARED = ../Shared
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
    OPT_COMP = -J $(DIR_MODULE) -fdefault-real-8 -O0 -g -Wall -fcheck=all
  else
    OPT_COMP = -J $(DIR_MODULE) -fdefault-real-8 -O3
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

$(info #---------------------------------------------------)
$(info # Compiling $(PROGRAM_NAME) with compiler $(FORTRAN))
$(info #---------------------------------------------------)

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
           Tokenizer_Mod.f90

# Sources for modules in the local directory
SRC_MOD += Parameters_Mod.f90		\
           Info_Mod.f90 		\
           Bulk_Mod.f90 		\
           Var_Mod.f90			\
           Matrix_Mod.f90		\
           Solvers_Mod.f90		\
           pro_mod.f90			\
           User_Mod.f90

# Sources for modules in sub-directories
SRC_MOD += les_mod.f90			\
           rans_mod.f90			\

# Sources for all functions are obtained by a shell command
SRC_RAW =  $(shell ls -1 *.f90               | grep -v -i _Mod)
SRC_RAW += $(shell ls -1 $(DIR_MODELS)/*.f90 | grep -v -i _Mod)
SRC_RAW += $(shell ls -1 $(DIR_COMMUN)/*.f90 | grep -v -i _Mod)

# Remove directory names
STEP_A  = $(subst $(DIR_MODELS),,$(SRC_RAW)) # remove directory name
STEP_B  = $(subst $(DIR_COMMUN),,$(STEP_A))  # remove directory name
SRC_FUN = $(subst /,,$(STEP_B))              # remove last slash 

# Add files from "Shared" directory, but one by one, to avoid clutter 
SRC_FUN += Distance.f90			\
           Load_Geo.f90			\
           Name_File.f90		\
           Sort_Int_Carry_Int.f90	\
           To_Upper_Case.f90

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
	$(FC) $(OPT_LINK) -o $(PROGRAM_FILE) $(OBJ)

#---------------------------------------
#   Explicit dependencies for modules
#---------------------------------------
#   These should be inserted by 
#   hand for tuning of dependencies.
#---------------------------------------
$(DIR_OBJECT)/all_mod.o:	$(DIR_OBJECT)/allp_mod.o
$(DIR_OBJECT)/Grid_Mod.o:	$(DIR_OBJECT)/Boundary_Condition_Mod.o	\
				$(DIR_OBJECT)/Material_Mod.o
$(DIR_OBJECT)/Domain_Mod.o:	$(DIR_OBJECT)/Point_Mod.o		\
  				$(DIR_OBJECT)/Line_Mod.o		\
  				$(DIR_OBJECT)/Block_Mod.o		\
  				$(DIR_OBJECT)/Region_Mod.o

#---------------------
#   Explicit target.
#---------------------
clean:
	rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod $(PROGRAM_FILE)
