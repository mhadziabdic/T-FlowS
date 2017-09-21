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

#-------------------------
  PROGRAM = RUN 
#-------------------------

#=============================
#....Fortran compiler & linker
#-----------------------------

  ifeq ($(TARGET),cluster)
    FC     = gfortran
    FL     = $(FC)
    ifeq ($(PRECISION),single)
      FCOPTS =  -fdefault-real-8 -fdefault-double-8 -g -fcheck=all
    else
      FCOPTS = -fdefault-real-8 -fdefault-double-8 -g -fcheck=all
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
  FSRCS0 := Neu2TFlowS.f90                                              \
            Calc4.f90 Dist.f90 FindSides.f90 Read9Ln.f90 DISort.f90   \
            ReadFluentNeu.f90  AA.f90 TetVol.f90 TopolM.f90 ReadC.f90 \
            GeoSav.f90 GenSav.f90 TestLn.f90 NamFil.f90 EpsSav.f90    \
            CasSav.f90 EpsWho.f90 ToUppr.f90 RISort.f90 Approx.f90    \
            CouMat.f90 Dist2.f90 Connect3.f90 Probe1D.f90 Probe2D.f90 SSort.f90\
            Probe1D_nodes.f90 Sort2.f90
  FSRCS  := $(sort $(FSRCS0))

#....C sources
#-------------
  CSRCS0 := logo_neu.c 
  CSRCS  := $(sort $(notdir $(CSRCS0)))

#....header files
#----------------
  INC0 := allp_mod.f90 all_mod.f90 gen_mod.f90 par_mod.f90 neu_mod.f90
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
	rm -f *.o *.vo *.inc *.mod RUN 
