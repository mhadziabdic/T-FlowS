#//////////////////////////////////////////////////////////////////////#
#              ___________ __________ _______________  ___ 
#              \__    ___/ \______   \\_   _____/\   \/  /  
#                |    |     |       _/ |    __)_  \     /   
#                |    |     |    |   \ |        \ /     \   
#                |____|     |____|_  //_______  //___/\  \  
#                                  \/         \/       \_/  
#                            Generator Makefile
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

#-------
#- linux
#-------
#  ifeq ($(TARGET),linux)
#    FC     = mpixlf90 
#    FC     = f90 
#    FL     = $(FC)
#    ifeq ($(PRECISION),single)
#      FCOPTS =   -qrealsize=8 -dbl8 -O3 #-qautodbl=dbl4  #-qrealsize=8 -qintsize=8 
#    else
#      FCOPTS =   -qrealsize=8 -dbl8 -O3 #-qautodbl=dbl4  #-qrealsize=8 -qintsize=8 
#    endif
#    FLOPTS =   
#  endif

  ifeq ($(TARGET),linux)
    FC     = gfortran
    FL     = $(FC)
    ifeq ($(PRECISION),single)
      FCOPTS =  -fdefault-real-8 -fdefault-double-8 -O3
    else
      FCOPTS = -fdefault-real-8 -fdefault-double-8 -O3
    endif
    FLOPTS =
  endif


  ifeq ($(TARGET),linuxx)
    FC     = ifort
    FL     = $(FC)
    ifeq ($(PRECISION),single)
      FCOPTS =  -r8 -O4
    else
      FCOPTS = -r8 -O4
    endif
    FLOPTS =
  endif


#-------
#- cluster
#-------
  ifeq ($(TARGET),cluster)
    FC     = gfortran  
    FL     = $(FC)
    ifeq ($(PRECISION),single)
      FCOPTS =  -fdefault-real-8 -fdefault-double-8
    else
      FCOPTS = -fdefault-real-8 -fdefault-double-8
    endif
    FLOPTS =   
  endif

#-----
#- sun
#-----
  ifeq ($(TARGET),sun)
    FC     = f90
    FL     = $(FC)
    ifeq ($(PRECISION),single)
      FCOPTS =
    else
      FCOPTS = -r8
    endif
    FLOPTS =
  endif

#----------
#- cray-t3e
#----------
  ifeq ($(TARGET),cray-t3e)
    FC     = f90  
    FL     = $(FC)
    FCOPTS =     
    FLOPTS =   
  endif

#-------
#- teras
#-------
  ifeq ($(TARGET),teras)
    FC     = f90
    FL     = $(FC)
    ifeq ($(PRECISION),single)
      FCOPTS =
    else
      FCOPTS = -r8
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
  FSRCS0 := Generator.f90                                           \
            Atanh.f90     Fuzion.f90    Linija.f90    PrintG.f90    \
            Calc1.f90     GenLoa.f90    Mark.f90      SwapI.f90     \
            Calc2.f90     CasSav.f90    PeriBC.f90    SwapR.f90     \
            EpsSav.f90    IsLine.f90    Refine.f90    TetVol.f90    \
            FinLin.f90    IsTwin.f90    IniGen.f90    TopSys.f90    \
            FinSur.f90    Laplac.f90    Smooth.f90    WchNod.f90    \
            CopyBC.f90    Probe1D.f90   CouMat.f90    ToUppr.f90    \
            Dist.f90      TestLn.f90    ReadC.f90     SSort.f90     \
            GenSav.f90    GeoSav.f90    Dist2.f90     ISort.f90     \
            ToUppr.f90    NamFil.f90    Approx.f90    EpsWho.f90    \
            RISort.f90    AA.f90        Probe2D.f90   Probe1D_nodes_gen.f90 Sort2.f90
  FSRCS  := $(sort $(FSRCS0))

#....C sources
#-------------
  CSRCS0 := logo_gen.c 
  CSRCS  := $(sort $(notdir $(CSRCS0)))

#....header files
#----------------
  INC0 := allp_mod.f90 all_mod.f90 gen_mod.f90 
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
	rm -f *.o *.vo *.inc *.mod RUN Generator.real
