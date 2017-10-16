bin_PROGRAMS = Generate

Generate_SOURCES =	 	\
  Are_Nodes_Twins.f90		\
  Atanh.f90			\
  Compute_Node_Coordinates.f90	\
  Compute_Grid_Geometry.f90	\
  CopyBC.f90			\
  Distribute_Nodes.f90		\
  Find_Line.f90			\
  Find_Surface.f90		\
  Connect_Blocks.f90		\
  Load_Domain.f90		\
  Is_Line_In_Block.f90		\
  Laplac.f90			\
  Logo.f90			\
  Main.f90			\
  PeriBC.f90			\
  Print_Grid_Statistics.f90	\
  Refine_Grid.f90		\
  Refine_Marked_Cells.f90	\
  Smooth.f90			\
  Swap_Integers.f90		\
  Swap_Reals.f90		\
  TopSys.f90			\
  Which_Node.f90

SUBDIRS = ../Modules		\
          ../Shared

LDADD = ../Modules/libModules.a	\
        ../Shared/libShared.a

INCLUDE = -I@top_srcdir@/
