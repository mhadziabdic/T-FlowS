bin_PROGRAMS = Process

Process_SOURCES =	 			\
  CalBou.f90					\
  Calc3.f90					\
  CalcConvect.f90				\
  CalcF22.f90					\
  CalcFlux.f90					\
  CalcG.f90					\
  CalcMn_Cylind.f90				\
  CalcMn.f90					\
  CalcPF.f90					\
  CalcPS.f90					\
  CalcSc.f90					\
  CalcSGS_Dynamic.f90				\
  CalcSGS.f90					\
  CalcSGS_hybrid.f90				\
  CalcShear.f90					\
  CalcStresses.f90				\
  CalcTurb.f90					\
  CalcVISt_EBM.f90				\
  CalcVISt_KEps.f90				\
  CalcVISt_KEPSV2F.f90				\
  CalcVISt_SPA_ALL.f90				\
  CalcVort.f90					\
  CalcWALE.f90					\
  CalMinMax.f90					\
  ConvScheme.f90				\
  CorBad.f90					\
  CorUVW.f90					\
  FindBad.f90					\
  ForAPF.f90					\
  GeoAloc.f90					\
  GraCorNew.f90					\
  GradP3.f90					\
  GradP.f90					\
  GraPhiCor.f90					\
  GraPhi.f90					\
  IniPar.f90					\
  IniVar.f90					\
  Load_Boundary_Conditions.f90			\
  Load_Cns.f90					\
  Load_Ini.f90					\
  Load_Restart.f90				\
  Load_Restart_Ini.f90				\
  Logo_Process.f90				\
  Main.f90   					\
  ModOut.f90					\
  NearWallCell.f90				\
  NewUVW.f90					\
  Save_Dat_Results.f90				\
  Save_Dat_Scalar.f90				\
  Save_Ini.f90					\
  Save_Gmv_Results.f90				\
  Save_Restart.f90				\
  ReaCom.f90					\
  SavParView.f90				\
  Scale.f90					\
  SourceEpsKEps.f90				\
  SourceEpsKEPSV2F.f90				\
  SourceF22_EBM.f90				\
  SourceF22KEPSV2F.f90				\
  SourceKinKEps.f90				\
  SourceKinKEPSV2F.f90				\
  SourcesEBM.f90				\
  SourcesHJ.f90					\
  SourceV2KEPSV2F.f90				\
  SourceVisSpalart.f90				\
  UnkAloc.f90

SUBDIRS = ../Shared/Modules			\
          Communication				\
          Linear_Solvers			\
          User					\
          ../Shared

LDADD = ../Shared/Modules/libSharedModules.a	\
        Communication/libCommunication.a	\
        Linear_Solvers/libLinear_Solvers.a	\
        User/libUser.a				\
        ../Shared/libShared.a

INCLUDE = -I@top_srcdir@/
