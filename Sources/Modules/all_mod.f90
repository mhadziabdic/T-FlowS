!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!         for all modules         !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!   Note: cell_n, parent, A_row, A_col, A_dia, side_c, side_cc, 
!         sideAij, are for all grids
!======================================================================!
module all_mod

  implicit none

  real,allocatable :: xc(:),yc(:),zc(:) 
  real,allocatable :: Sx(:),Sy(:),Sz(:)
  real,allocatable :: volume(:)            ! cell's volume
  real,allocatable :: delta(:)             ! delta (max(dx,dy,dz))
  real,allocatable :: a1(:)                ! scotti
  real,allocatable :: a2(:)                ! scotti
  real,allocatable :: Dx(:),Dy(:),Dz(:)
  real,allocatable :: xsp(:),ysp(:),zsp(:) ! face coordinates    
  real,allocatable :: WallDs(:), f(:)

  character :: name*80
  character :: inp*300
  integer   :: tn, ts(300), te(300)

  integer   :: cmn_line_count

  integer   :: NC, NS                    ! num. of nodes and cells 
  integer   :: NbC
  integer   :: MNBS
  integer   :: NRL
  integer   :: Ncopy
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
