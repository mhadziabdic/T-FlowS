!+++++++++++++++++++++++++++++++++++++!
!                                     !
!     Global variable definitions     !
!       for the mesh generator        !
!                                     !
!+++++++++++++++++++++++++++++++++++++!
module gen_mod

  implicit none


  real,    allocatable :: walln(:)           ! node distance from the wall 
  integer, allocatable :: SideN(:,:)         ! numb, n1, n2, n3, n4
  integer, allocatable :: SideCc(:,:)
                                                
  integer, allocatable :: TwinN(:,:)

  integer, allocatable :: NewN(:)    ! new number for the nodes and cells
  integer, allocatable :: NewC(:)    ! new number for cells
  integer, allocatable :: NewS(:)    ! new number for sides 
  integer, allocatable :: CelMar(:)  ! cell marker

  ! variables for refinement
  integer, parameter   :: ELIPSOID  = 3
  integer, parameter   :: RECTANGLE = 4
  integer, parameter   :: PLANE     = 5
  integer              :: n_refine_levels
  integer, allocatable :: level(:)   ! refinement level
  integer, allocatable :: n_refined_regions(:)    ! number of refin. regions
  real,    allocatable :: refined_regions(:,:,:)  ! levels, regions

  ! variables for smoothing
  integer              :: n_smoothing_regions   
  logical, allocatable :: smooth_in_x(:), smooth_in_y(:), smooth_in_z(:)   
  integer, allocatable :: smooth_iters(:)   
  real,    allocatable :: smooth_regions(:,:), smooth_relax(:)

  integer              :: n_periodic_cond,     &
                          n_copy_cond
  integer, allocatable :: periodic_cond(:,:),  &                         
                          copy_cond(:,:)

  integer :: Nsurf
  integer :: NN
  integer :: NSsh                  ! number of shadow faces

  integer :: WallFacFst, WallFacLst 

end module
