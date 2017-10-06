!+++++++++++++++++++++++++++++++++++++!
!                                     !
!     Global variable definitions     !
!       for the mesh generator        !
!                                     !
!+++++++++++++++++++++++++++++++++++++!
module gen_mod

  use allp_mod

  implicit none

  integer, parameter :: YES       = 1
  integer, parameter :: NO        = 2
  integer, parameter :: ELIPSOID  = 3
  integer, parameter :: RECTANGLE = 4
  integer, parameter :: PLANE     = 5

  real, allocatable :: x_node(:),  &  ! node coordinates
                       y_node(:),  &
                       z_node(:) 

  real,    allocatable :: walln(:)           ! node distance from the wall 
  integer, allocatable :: SideN(:,:)         ! numb, n1, n2, n3, n4
  integer, allocatable :: SideCc(:,:)
                                                
  integer, allocatable :: CellC(:,:)        ! cell's neighbours
  integer, allocatable :: CellN(:,:)        ! cell nodes

  integer, allocatable :: TwinN(:,:)

  integer, allocatable :: NewN(:)    ! new number for the nodes and cells
  integer, allocatable :: NewC(:)    ! new number for cells
  integer, allocatable :: NewS(:)    ! new number for sides 
  integer, allocatable :: CelMar(:)  ! cell marker

  integer, allocatable :: NodeN2(:,:)    
  integer, allocatable :: NodeN4(:,:)  
  integer, allocatable :: NodeN8(:,:)    

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

  real, allocatable :: x_point(:), &  ! point coordinates
                       y_point(:), &
                       z_point(:)      

  real, allocatable :: xl(:,:), yl(:,:),zl(:,:), LinWgt(:)
  real              :: BlkWgt(MAXL,3), BlFaWt(MAXL,3)     ! leave this 

  integer, allocatable :: block_points(:,:),       &  ! 0 for orientation                
                          block_resolutions(:,:),  &  ! NI,NJ,NK,NI*NJ*NK,NNo,NV
                          block_faces(:,:,:),      & 
                          face_laplace(:),         &  ! never used 
                          b_cond(:,:),             &  
                          periodic_cond(:,:),      &                         
                          copy_cond(:,:)

  integer :: LinPnt(MAXL,2), LinRes(MAXL)
  integer :: Nbloc, Nline, Nsurf, n_b_cond, n_periodic_cond, n_copy_cond, n_refine_levels
  integer :: NN
  integer :: NSsh                  ! number of shadow faces

  integer :: WallFacFst, WallFacLst 

  character(len=4), allocatable :: BndFac(:)

end module
