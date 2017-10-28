!==============================================================================!
  module Cell_Mod
!------------------------------------------------------------------------------!
!   Cells defining blocks used in "Generator"                                  !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Cell type   !
  !---------------!
  type Cell_Type

    integer :: n_nodes  ! number of cell's nodes
    integer :: c(24)    ! cell's neighboring cells (in case of refined cells)
    integer :: n(8)     ! cell's nodes 

    ! These are cell center coordinates
    real :: x
    real :: y
    real :: z

  end type Cell_Type

  end module Cell_Mod
