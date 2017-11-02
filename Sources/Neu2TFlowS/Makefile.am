bin_PROGRAMS = Neu2TFlowS

Neu2TFlowS_SOURCES =		\
  Compute_Geometry.f90		\
  Connect_Domains.f90		\
  Find_Faces.f90		\
  Logo.f90			\
  Main.f90			\
  Load_Neu.f90			\
  Grid_Topology.f90		\
  Read_Line.f90

SUBDIRS = ../Shared/Modules	\
          ../Shared

LDADD = ../Shared/Modules/libSharedModules.a	\
        ../Shared/libShared.a

INCLUDE = -I@top_srcdir@/
