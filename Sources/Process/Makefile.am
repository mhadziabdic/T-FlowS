bin_PROGRAMS = Process

Process_SOURCES =	 			\
  Allocate_Memory.f90				\
  Allocate_Variables.f90			\
  Balance_Mass.f90				\
  CalBou.f90					\
  CalcConvect.f90				\
  CalcG.f90					\
  CalcMn_Cylind.f90				\
  CalcMn.f90					\
  CalcPF.f90					\
  CalcPS.f90					\
  CalcShear.f90					\
  CalcVort.f90					\
  CalMinMax.f90					\
  Compute_Fluxes.f90				\
  Compute_Geometry.f90				\
  Compute_Scalar.f90				\
  ConvScheme.f90				\
  CorBad.f90					\
  Correct_Velocity.f90				\
  FindBad.f90					\
  ForAPF.f90					\
  GraCorNew.f90					\
  GradP3.f90					\
  GradP.f90					\
  GraPhiCor.f90					\
  GraPhi.f90					\
  IniPar.f90					\
  Initialize_Variables.f90			\
  Load_Boundary_Conditions.f90			\
  Load_Cns.f90					\
  Load_Ini.f90					\
  Load_Restart.f90				\
  Load_Restart_Ini.f90				\
  Logo_Process.f90				\
  Main.f90   					\
  NearWallCell.f90				\
  NewUVW.f90					\
  Save_Dat_Results.f90				\
  Save_Dat_Scalar.f90				\
  Save_Ini.f90					\
  Save_Gmv_Results.f90				\
  Save_Restart.f90				\
  ReaCom.f90					\
  SavParView.f90				\
  Scale.f90

SUBDIRS = ../Shared/Modules			\
          Communication				\
          Solvers				\
          Turbulence				\
          User					\
          ../Shared

LDADD = ../Shared/Modules/libSharedModules.a	\
        Communication/libCommunication.a	\
        Solvers/libSolvers.a			\
        Turbulence/libTurbulence.a		\
        User/libUser.a				\
        ../Shared/libShared.a

INCLUDE = -I@top_srcdir@/
