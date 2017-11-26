noinst_LIBRARIES = libSolvers.a

libSolvers_a_SOURCES =		\
  Bicg.f90			\
  Cg.f90			\
  Cgs.f90			\
  Prec_Form.f90			\
  Prec_Solve.f90		\
  Residual.f90

SUBDIRS = ../../Shared

LDADD = ../../Shared/libShared.a

INCLUDE = -I@top_srcdir@/
