bin_PROGRAMS = Divide

Divide_SOURCES =		\
  Load_Cns.f90			\
  Load_Gmv.f90			\
  Load_Gmv_Faces.f90		\
  Factorize_Number.f90		\
  Inertia.f90			\
  Compute_Eigenvalues.f90	\
  Logo.f90			\
  Main.f90			\
  Number.f90			\
  Save_Com.f90			\
  Save_Eps_Decomposed.f90	\
  Save_Scripts.f90		\
  Split_Subdomain.f90

SUBDIRS = ../Modules		\
          ../Shared

LDADD = ../Modules/libModules.a	\
        ../Shared/libShared.a

INCLUDE = -I@top_srcdir@/
