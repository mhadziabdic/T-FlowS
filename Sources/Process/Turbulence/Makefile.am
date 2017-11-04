noinst_LIBRARIES = libTurbulence.a

libTurbulence_a_SOURCES =	\
  les_mod.f90			\
  rans_mod.f90			\
  CalcF22.f90			\
  CalcSGS_Dynamic.f90		\
  CalcSGS.f90			\
  CalcSGS_hybrid.f90		\
  CalcVISt_EBM.f90		\
  CalcVISt_KEps.f90		\
  CalcVISt_KEPSV2F.f90		\
  CalcVISt_SPA_ALL.f90		\
  CalcWALE.f90			\
  Compute_Stresses.f90		\
  Compute_Turbulent.f90		\
  SourceEpsKEps.f90		\
  SourceEpsKEPSV2F.f90		\
  SourceF22_EBM.f90		\
  SourceF22KEPSV2F.f90		\
  SourceKinKEps.f90		\
  SourceKinKEPSV2F.f90		\
  SourcesEBM.f90		\
  SourcesHJ.f90			\
  SourceV2KEPSV2F.f90		\
  SourceVisSpalart.f90

SUBDIRS = ../../Shared

LDADD = ../../Shared/libShared.a

INCLUDE = -I@top_srcdir@/
