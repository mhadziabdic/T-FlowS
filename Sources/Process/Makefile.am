bin_PROGRAMS = Process

Process_SOURCES =	 			\
  BouLoa.f90					\
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
  CnsLoa.f90					\
  ConvScheme.f90				\
  CorBad.f90					\
  CorUVW.f90					\
  DatSav.f90					\
  DatSavSc.f90					\
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
  LoaIni.f90					\
  LoaRes.f90					\
  LoaRes_Ini.f90				\
  Logo_Process.f90				\
  ModOut.f90					\
  NearWallCell.f90				\
  NewUVW.f90					\
  Process.f90					\
  ProSav.f90					\
  ReaCom.f90					\
  SavIni.f90					\
  SavParView.f90				\
  SavRes.f90					\
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

SUBDIRS = ../Modules				\
          Communication				\
          Linear_Solvers			\
          User					\
          ../Shared

LDADD = ../Modules/libModules.a			\
        Communication/libCommunication.a	\
        Linear_Solvers/libLinear_Solvers.a	\
        User/libUser.a				\
        ../Shared/libShared.a

INCLUDE = -I@top_srcdir@/
