#===============================================================================
# 
#   Neu2TflowS Makefile   
#   
#-------------------------------------------------------------------------------

#--------------------------
#   Variable definitions
#--------------------------

# Fortran compiler ("gnu", "intel" or "portland")
FORTRAN = gnu

# Debugging ("yes" or "no")
DEBUG = no

# Directories for objects and modules. (No need to change.)
DIR_BINARY = ../../Binaries
DIR_SHARED = ../Shared
DIR_MODULE = ../.Modules
DIR_OBJECT = .Objects
VPATH = $(DIR_SHARED):.

# Program name (This should hardly change)
PROGRAM_NAME = Neu2TflowS
PROGRAM_FILE = $(DIR_BINARY)/$(PROGRAM_NAME)

#-------------------------------------------------------------------------
#   Compiler and linker options
#-------------------------------------------------------------------------
#   Note: Changes only when support to a new Fortran compiler is added.
#-------------------------------------------------------------------------

# Fortran = gnu
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

$(info #====================================================)
$(info # Compiling $(PROGRAM_NAME) with compiler $(FORTRAN) )
$(info #----------------------------------------------------)
$(info # Usage:                                             )
$(info #   make <FORTRAN=gnu/intel/portland> <DEBUG=yes/no> )
$(info #----------------------------------------------------)

#------------------------------------------------------
#   List of sources for modules and functions
#------------------------------------------------------
#   Modules' order must obey their dependency 
#   This list should therefore be written "by hand".
#   Note: Modules written in lower case 
#         letters are candidates for deletion.
#------------------------------------------------------

# Sources for modules in the local directory
SRC_MOD =  neu_mod.f90

# Sources for modules in other directories
SRC_MOD += allp_mod.f90			\
           all_mod.f90			\
           gen_mod.f90			\
           div_mod.f90			\
           Tokenizer_Mod.f90		\
           Material_Mod.f90		\
           Boundary_Condition_Mod.f90	\
           Grid_Mod.f90

# Sources for all functions are obtained by a shell command
SRC_FUN_WITH_DIR = $(shell ls -1 $(DIR_SHARED)/*.f90 *.f90 | grep -v -i _Mod)
SRC_FUN = $(subst $(DIR_SHARED)/, ,$(SRC_FUN_WITH_DIR))

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

#---------------------
#   Explicit target.
#---------------------
clean:
	rm -f $(DIR_OBJECT)/*.o $(DIR_MODULE)/*.mod $(PROGRAM_FILE)
