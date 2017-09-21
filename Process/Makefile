#//////////////////////////////////////////////////////////////////////#
#              ___________ __________ _______________  ___ 
#              \__    ___/ \______   \\_   _____/\   \/  /  
#                |    |     |       _/ |    __)_  \     /   
#                |    |     |    |   \ |        \ /     \   
#                |____|     |____|_  //_______  //___/\  \  
#                                  \/         \/       \_/  
#                            Processor Makefile
#----------------------------------------------------------------------#
#  You must have enviroment variable TARGET set to specific platform.
#  In csh you can do do it by inserting the line:
#    setenv TARGET [machine]
#  where [machine is: gnu, cluster, intel, ....
#----------------------------------------------------------------------#

TARGET = cluster

STATIC_LIBS =

# dirs where to search mod files while compiling
HEADER_DIRS = $(DMOD) #../Library/VTKFortran-1.1.1/mod/

# where to store compiled .o and .mod files
DOBJ  = ./obj/
DMOD  = ./mod/

ifeq ($(PRECISION),single)
	INCLUDE_DIRS = ../Parallel/Single
else
	INCLUDE_DIRS = ../Parallel/Double
endif

ifeq ($(FORCE),seq)
	INCLUDE_DIRS = ../Sequential
endif

vpath %.f90 ../Library/ ../User ./  $(INCLUDE_DIRS) # source .f90 files search pattern
vpath %.c   ../Library/ ../User ./  $(INCLUDE_DIRS) # source .c   files search pattern
#-------------------------
PROGRAM = RUN
#-------------------------

#=============================
#....Fortran compiler & linker
#-----------------------------

ifeq ($(TARGET),gnu)
	FC     = mpif90
	FL     = $(FC)
	ifeq ($(DEBUG),yes)
		FCOPTS = -g -O0 -J$(DMOD) -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace -fmax-errors=0 -Wall -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized -fcheck=all -fmodule-private -fimplicit-none -ffpe-trap=invalid,overflow -fbacktrace -fdump-core -finit-real=nan #-fno-range-check  ,precision,denormal,underflow -std=f2003 -fall-intrinsics
	else
		FCOPTS =    -O3 -J$(DMOD) -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none
	endif 
	FLOPTS =    
endif

ifeq ($(TARGET),intel)
	FC     = mpiifort
	FL     = $(FC)
	ifeq ($(DEBUG),yes)
		FCOPTS = -g -O0 -module $(DMOD) -r8 -traceback -warn all -check all -debug all -extend-source 132 -fpe-all=0 -fp-stack-check -fstack-protector-all -ftrapuv -no-ftz -traceback -gen-interfaces -std03 -ipo -inline all -ipo-jobs4 -vec-report1
	else
		FCOPTS =    -O3 -module $(DMOD) -r8
	endif 
	FLOPTS =    
endif

ifeq ($(TARGET),cluster)
  FC     = mpif90
  FL     = $(FC)
  ifeq ($(PRECISION),single)
    FCOPTS =  -J$(DMOD) -fdefault-real-8 -fdefault-double-8 -O3  -I /usr/lib64/mpi/gcc/openmpi/include/
  else
    FCOPTS = -J$(DMOD) -fdefault-real-8 -fdefault-double-8  -O3  -I /usr/lib64/mpi/gcc/openmpi/include/
  endif
  FLOPTS = #-L/home/muhamed/Downloads/openmpi-1.8.3/ompi/include  
endif

#....Fortran sources
#-------------------
  FSRCS := allp_mod.f90\
			all_mod.f90\
			pro_mod.f90\
			les_mod.f90\
			par_mod.f90\
			sol_mod.f90\
			rans_mod.f90\
			Processor.f90                                                   \
            BiCG.f90        CorUVW.f90      CnsLoa.f90      CGS.f90         \
            BouLoa.f90      ProSav.f90      Resid.f90       \
            Calc3.f90       ForAPF.f90      CalcSc.f90      \
            CalBou.f90      GraPhi.f90      ReaCom.f90      TopolM.f90      \
            CalcMn.f90      IniPar.f90      CG.f90          GradP.f90       \
            CalcPF.f90      IniVar.f90      SavRes.f90      SourceKinKEps.f90 \
            CalcPS.f90      LoaRes.f90      CalcG.f90       SourceEpsKEps.f90 \
            ModOut.f90      NewUVW.f90      CalcShear.f90   DatSav.f90      \
            Dist.f90        CalcTurb.f90    ReadC.f90       GeoAloc.f90     \
            ToUppr.f90      Rnd.f90         GeoLoa.f90      ISort.f90       \
            Prec1.f90	    NamFil.f90      GradP3.f90      UnkAloc.f90     \
            StaPar.f90      EndPar.f90      BufLoa.f90      Exchng.f90      \
            GloMax.f90      GloMin.f90      GloSum.f90      IGlSum.f90      \
            FindBad.f90     CorBad.f90      CalcFlux.f90    CalcConvect.f90 \
            Wait.f90        Timex.f90       Approx.f90               \
            UserUPlus.f90   UserURms.f90    UserVRms.f90    UserWRms.f90   \
            UserCutLines_budgets_cylind.f90 UserCutLines_budgets_cylind_HOT.f90\
            CalMinMax.f90   ConvScheme.f90  CalcVISt_KEps.f90              \
            CalcSGS.f90     UserPerturb.f90 UserPerturb2.f90 \
            UserForce.f90   CalcBudgets_cylind.f90  \
            GraCorNew.f90   UserSource.f90 \
            SourceF22KEPSV2F.f90 SourceF22_EBM.f90 UserOut.f90 \
            SourceEpsKEPSV2F.f90 SourceKinKEPSV2F.f90 Scale.f90 UserStream.f90 \
            SourceV2KEPSV2F.f90 CalcF22.f90 DatSavSc.f90 \
            CalcVISt_KEPSV2F.f90 SourceVisSpalart.f90 \
            CalcVISt_SPA_ALL.f90 CalcVort.f90 RNSort.f90 RISort.f90\
            CalcSGS_Dynamic.f90  NearWallCell.f90 LoaIni.f90 SavIni.f90\
            CalcMn_Cylind.f90 UserCutLines_channel.f90 SSort.f90\
            UserCutLines_pipe.f90 UserCutLines_Nu.f90 DistR.f90 CalcWALE.f90\
            UserDiffuser.f90 UserProbe1D_Nusselt_jet.f90 UserProbe1D_jet_rans.f90\
            SavParView.f90 CalcStresses.f90\
	    CalcVISt_EBM.f90 SourcesEBM.f90 CalcSGS_hybrid.f90 UserCutLines_annulus.f90\
            UserCutLines_Horiz_RANS.f90 SourcesHJ.f90 UserPlainGen_Tecplot.f90 UserBackstep.f90\
            UserCutLines_Y_dir.f90 UserBackstep_Y.f90 UserCutLines_jet.f90 UserCutLines_Point_Cart.f90\
            UserCutLines_RB_conv.f90\
            LoaRes_Ini.f90 

#....C sources
#-------------
  CSRCS := logo_pro.c 

#-----
FOBJS := $(FSRCS:%.f90=$(DOBJ)%.o)
COBJS := $(CSRCS:%.c=$(DOBJ)%.o)
  
  

#//////////////////////////////////////////////////////////////////////#
all      : $(PROGRAM)
#----------------------------------------------------------------------#
$(PROGRAM): $(FOBJS) $(COBJS)
	@echo "   Linking $(PROGRAM) ..."
	$(FL) -o $(PROGRAM) $(FOBJS) $(STATIC_LIBS) $(COBJS) $(FLOPTS)

$(FOBJS): $(DOBJ)%.o: %.f90
	$(FC) -c $(FCOPTS) $(addprefix -I, $(HEADER_DIRS)) -o $@ $<

$(COBJS): $(DOBJ)%.o: %.c
	$(CC) -c $(CCOPTS) $(addprefix -I, $(HEADER_DIRS)) -o $@ $<

#======================================================================#
#...explicit targets
#-----------------------------
clean:
	rm -f $(DOBJ)*.o *.vo *.inc $(DMOD)*.mod *__genmod.f90 $(PROGRAM) *.real *.optrpt
