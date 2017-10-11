noinst_LIBRARIES = libUser.a

libUser_a_SOURCES =	\
  UserForce.f90		\
  UserSource.f90           \  
  UserCutLines_channel.f90 \
  UserCutLines_pipe.f90    \
  UserCutLines_annulus.f90 \
  UserCutLines_Nu.f90      \
  UserCutLines_jet.f90     \
  UserProbe1D_Nusselt_jet.f90  \
  UserCutLines_Horiz_RANS.f90  \
  UserBackstep.f90             \
  UserBackstep_Y.f90           \    
  UserCutLines_RB_conv.f90     \
  UserCutLines_Point_Cart.f90

SUBDIRS = ../../Modules

LDADD = ../../Modules/libModules.a

INCLUDE = -I@top_srcdir@/
