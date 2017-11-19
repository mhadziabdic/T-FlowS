noinst_LIBRARIES = libTurbulence.a

libTurbulence_a_SOURCES =	\
  les_mod.f90			\
  rans_mod.f90			\
  Compute_Sgs_Dynamic.f90	\
  Compute_Sgs_Hybrid.f90	\
  CalcVISt_KEps.f90		\
  CalcVISt_KEPSV2F.f90		\
  CalcVISt_RSM.f90		\
  CalcVISt_SPA_ALL.f90		\
  CalcWALE.f90			\
  Compute_F22.f90		\
  Compute_Sgs.f90		\
  Compute_Stresses.f90		\
  Compute_Turbulent.f90		\
  Source_Kin_K_Eps.f90		\
  Source_Eps_K_Eps.f90		\
  Source_Eps_K_Eps_V2_F.f90	\
  Source_F22_K_Eps_V2_F.f90	\
  Source_Kin_K_Eps_V2_F.f90	\
  Source_V2_K_Eps_V2_F.f90	\
  Source_Ebm.f90		\
  Source_F22_Ebm.f90		\
  Source_Hanjalic_Jakirlic.f90	\
  Source_Vis_Spalart_Almaras.f90

SUBDIRS = ../../Shared

LDADD = ../../Shared/libShared.a

INCLUDE = -I@top_srcdir@/
