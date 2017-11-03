noinst_LIBRARIES = libSolvers.a

libSolvers_a_SOURCES =		\
  BiCG.f90			\
  CG.f90			\
  CGS.f90			\
  Matrix_Topology.f90           \
  Prec_Form.f90			\
  Prec_Solve.f90		\
  Resid.f90    

SUBDIRS = ../../Shared/Modules

LDADD = ../../SharedModules/libSharedModules.a

INCLUDE = -I@top_srcdir@/
