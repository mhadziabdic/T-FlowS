!==============================================================================!
!   Module for global variable definitions.  Variables defined here are        !
!   supposed to be accessible to all routines in all programs.                 !
!==============================================================================!
module all_mod

  use allp_mod

  implicit none

  !----------------------------------------------------!
  !   Geometrical quantities for describing the grid   !
  !----------------------------------------------------!
  real,allocatable :: xc(:),yc(:),zc(:) 
  real,allocatable :: Sx(:),Sy(:),Sz(:)
  real,allocatable :: volume(:)            ! cell's volume
  real,allocatable :: delta(:)             ! delta (max(dx,dy,dz))
  real,allocatable :: a1(:)                ! scotti
  real,allocatable :: a2(:)                ! scotti
  real,allocatable :: Dx(:),Dy(:),Dz(:)
  real,allocatable :: xsp(:),ysp(:),zsp(:) ! face coordinates    
  real,allocatable :: WallDs(:), f(:)

  !----------------------------------------!
  !   Variables for ease of input/output   !
  !----------------------------------------!
  character(len=80)  :: name
  character(len=300) :: inp*300
  integer            :: tn, ts(300), te(300)
  integer            :: cmn_line_count

  !-------------------------------------------!
  !   Logical quantities desribing the grid   !
  !-------------------------------------------!
  integer   :: NC, NS                    ! num. of nodes and cells 
  integer   :: NbC
  integer   :: MNBS
  integer   :: NRL
  integer   :: n_copy                    ! number of copy cells/faces
  integer   :: Nmat                      ! number of materials
  logical   :: Mater(1024)               ! is the material present ?

  integer,allocatable :: material(:)     ! material markers
  integer,allocatable :: SideC(:,:)      !  c0, c1, c2

  integer,allocatable :: TypeBC(:)       ! type of boundary condition
  integer,allocatable :: bcmark(:)

  integer,allocatable :: CopyC(:)        !  might be shorter
  integer,allocatable :: CopyS(:,:)      !  similar to SideC 

  integer,allocatable :: SideC1C2(:,:)   !  similar to SideC 

  real, allocatable   :: Dxsp(:,:)       !  similar to SideC 
  real, allocatable   :: Dysp(:,:)       !  similar to SideC 
  real, allocatable   :: Dzsp(:,:)       !  similar to SideC 

end module 
