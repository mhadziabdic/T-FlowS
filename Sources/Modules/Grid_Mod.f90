!==============================================================================!
  module Grid_Mod
!------------------------------------------------------------------------------!
!   Grids defining blocks used in "Generator"                                  !
!------------------------------------------------------------------------------!
  use Point_Mod
  use Cell_Mod
  use Material_Mod
  use Boundary_Condition_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Grid type   !
  !---------------!
  type Grid_Type

    type(Point_Type),              allocatable :: nodes(:)
    type(Cell_Type),               allocatable :: cells(:)  

    type(Material_Type),           allocatable :: materials(:)
    type(Boundary_Condition_Type), allocatable :: boundary_conditions(:)

  end type Grid_Type

  ! If defined like this, one can easily think of multiple grids
  type(Grid_Type) :: grid

  contains

!==============================================================================!
  subroutine Sort_Cells_By_Index(cells, indx, n)
!------------------------------------------------------------------------------!
!   Sorts array of cells according to indx.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Cell_Type) :: cells(n)
  integer         :: n, indx(n)
!-----------------------------------[Locals]-----------------------------------!
  integer                      :: i
  type(Cell_Type), allocatable :: work(:)
!==============================================================================!

  allocate(work(n))

  do i = 1, n
    work( indx(i) ) = cells(i)
  end do

  do i=1,N
    cells(i) = work(i)
  end do

  deallocate(work)

  end subroutine

!==============================================================================!
  subroutine Allocate_Grid_Cells(grid, nb, ni)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nb        ! number of cells on the bounday
  integer         :: ni        ! number of cells inside
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  allocate(grid % cells(-nb:ni))

  ! Initialize nodes and cells
  do c = -nb, ni
    grid % cells(c) % n_nodes = 0 
    grid % cells(c) % n       = 0 
    grid % cells(c) % c       = 0 
  end do

  ! Initialize coordinates
  grid % cells(-nb:ni) % x = 0 
  grid % cells(-nb:ni) % y = 0 
  grid % cells(-nb:ni) % z = 0 

  end subroutine

!==============================================================================!
  subroutine Allocate_Grid_Nodes(grid, n)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: n
!==============================================================================!

  allocate(grid % nodes(n))

  grid % nodes(1:n) % x = 0.0 
  grid % nodes(1:n) % y = 0.0 
  grid % nodes(1:n) % z = 0.0 

  end subroutine

  end module Grid_Mod
