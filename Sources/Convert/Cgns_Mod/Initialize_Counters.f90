!==============================================================================!
  subroutine Initialize_Counters
!------------------------------------------------------------------------------!
!   Opens name_in file and return file index                                   !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  cnt_nodes = 0
  cnt_cells = 0
  cnt_hex   = 0
  cnt_pyr   = 0
  cnt_wed   = 0
  cnt_tet   = 0
  cnt_tri   = 0
  cnt_qua   = 0

  cnt_x = 0
  cnt_y = 0
  cnt_z = 0

  last_hex = 0
  last_pyr = 0
  last_wed = 0
  last_tet = 0
  last_tri = 0
  last_qua = 0

  end subroutine
