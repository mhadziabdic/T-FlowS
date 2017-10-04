#//////////////////////////////////////////////////////////////////////#
#              ___________ __________ _______________  ___ 
#              \__    ___/ \______   \\_   _____/\   \/  /  
#                |    |     |       _/ |    __)_  \     /   
#                |    |     |    |   \ |        \ /     \   
#                |____|     |____|_  //_______  //___/\  \  
#                                  \/         \/       \_/  
#                             Divisor Makefile
#----------------------------------------------------------------------#
#  You must have enviroment variable TARGET set to specific platform.
#  In csh you can do do it by inserting the line:
#    setenv TARGET [machine]
#  where [machine is: linux, hp, cray-t3e.
#----------------------------------------------------------------------#
  INC_SRCDIRSH = ../Library
  INCDIR = $(INC_SRCDIRSH)
  vpath %.f90 $(INC_SRCDIRSH)
  vpath %.c $(INC_SRCDIRSH)
  vpath %.f90 $(INC_SRCDIRSH)

#---------------------------
# Default target is cluster
#---------------------------
TARGET = cluster

#-------------------------
PROGRAM = RUN
#-------------------------

#=============================
#....Fortran compiler & linker
#-----------------------------

#----
#- hp
#----
  ifeq ($(TARGET),hp)
    FC     = f90  
    FL     = $(FC)
    ifeq ($(PRECISION),single)
      FCOPTS = 
    else
      FCOPTS = -R8 
    endif
    FLOPTS =   
  endif

  ifeq ($(TARGET),cluster)
    FC     = gfortran
    FL     = $(FC)
    ifeq ($(PRECISION),single)
      FCOPTS =  -fdefault-real-8 -fdefault-double-8 -O3
    else
      FCOPTS = -fdefault-real-8 -fdefault-double-8 -O3
    endif
    FLOPTS =
  endif


#======================= 
#....C compiler & linker
#-----------------------
  CC     = cc  
  ifeq ($(TARGET),hp)
    CCOPTS = -Aa   
  else
    CCOPTS =     
  endif

#//////////////////////////////////////////////////////////////////////#
  INC_SRCDIRS= $(subst :, ,$(INC_SRCDIRSH))

#....Fortran sources
#-------------------
  FSRCS0 := Divisor.f90                                             \
            DivLoa.f90   Factor.f90   Number.f90   CasSav.f90       \
            IniDiv.f90   Split.f90    RISort.f90   Distance.f90     \
            ReadC.f90    EpsPar.f90   GenSav.f90   ComSav.f90       \
            GeoLoa.f90   GeoSav.f90   TestLn.f90   NamFil.f90       \
            INSort.f90   BCelLoa.f90  ToUppr.f90   AA.f90           \
            Inertia.f90  Jacobi.f90   CouMat.f90   RNSort.f90
  FSRCS  := $(sort $(FSRCS0))

#....C sources
#-------------
  CSRCS0 := logo_div.c 
  CSRCS  := $(sort $(notdir $(CSRCS0)))

#....header files
#----------------
  INC0 := allp_mod.f90 all_mod.f90 gen_mod.f90 par_mod.f90 div_mod.f90
  INC  := $(sort $(notdir $(INC0)))

#-----
  FOBJS := $(FSRCS:.f90=.o)
  COBJS := $(CSRCS:.c=.o)

#//////////////////////////////////////////////////////////////////////#
all      : $(PROGRAM)
#----------------------------------------------------------------------#
$(PROGRAM): $(FOBJS) $(COBJS)
	@echo "   Linking $(PROGRAM) ..." 
	$(FL) -o $(PROGRAM) $(FLOPTS) $(sort $(FOBJS)  $(COBJS) )

$(FOBJS): %.o: %.f90 $(INC)
	$(FC) -c $(FCOPTS) -I$(INCDIR) $<

$(COBJS): %.o: %.c $(INC)
	$(CC) -c $(CCOPTS) -I$(INCDIR) $<

#======================================================================#
#...explicit targets
#-----------------------------
clean:
	rm -f *.o *.vo *.inc *.mod RUN Divisor.real
