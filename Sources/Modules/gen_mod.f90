!+++++++++++++++++++++++++++++++++++++!
!                                     !
!     Global variable definitions     !
!       for the mesh generator        !
!                                     !
!+++++++++++++++++++++++++++++++++++++!
module gen_mod

  use allp_mod

  implicit none

  integer, parameter :: ELIPSOID  = 3
  integer, parameter :: RECTANGLE = 4
  integer, parameter :: PLANE     = 5

  real,    allocatable :: walln(:)           ! node distance from the wall 
  integer, allocatable :: SideN(:,:)         ! numb, n1, n2, n3, n4
  integer, allocatable :: SideCc(:,:)
                                                
  integer, allocatable :: TwinN(:,:)

  integer, allocatable :: NewN(:)    ! new number for the nodes and cells
  integer, allocatable :: NewC(:)    ! new number for cells
  integer, allocatable :: NewS(:)    ! new number for sides 
  integer, allocatable :: CelMar(:)  ! cell marker

  integer, allocatable :: level(:)   ! refinement level

  integer :: MAXN, MAXB, MAXS

  ! variables for refinement
  integer, allocatable :: n_refined_regions(:)    ! number of refin. regions
  real,    allocatable :: refined_regions(:,:,:)  ! levels, regions

  ! variables for smoothing
  integer              :: n_smoothing_regions   
  logical, allocatable :: smooth_in_x(:), smooth_in_y(:), smooth_in_z(:)   
  integer, allocatable :: smooth_iters(:)   
  real,    allocatable :: smooth_regions(:,:), smooth_relax(:)

  integer, allocatable :: b_cond(:,:),             &  
                          periodic_cond(:,:),      &                         
                          copy_cond(:,:)

  integer :: Nsurf, n_b_cond, n_periodic_cond, n_copy_cond, n_refine_levels
  integer :: NN
  integer :: NSsh                  ! number of shadow faces

  integer :: WallFacFst, WallFacLst 

  character(len=4), allocatable :: BndFac(:)

end module
