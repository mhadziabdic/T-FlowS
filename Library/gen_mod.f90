!+++++++++++++++++++++++++++++++++++++!
!                                     !
!     Global variable definitions     !
!       for the mesh generator        !
!                                     !
!+++++++++++++++++++++++++++++++++++++!
MODULE gen_mod

  use allp_mod

  implicit none

  real,allocatable :: x(:),  y(:),  z(:)   ! node coordinates
  real,allocatable :: walln(:)             ! node distance from the wall 
  integer,allocatable :: SideN(:,:)        ! numb, n1, n2, n3, n4
  integer,allocatable :: SideCc(:,:)
                                                
  integer,allocatable :: CellC(:,:)        ! cell's neighbours
  integer,allocatable :: CellN(:,:)        ! cell nodes

  integer,allocatable :: TwinN(:,:)

  integer,allocatable :: NewN(:)    ! new number for the nodes and cells
  integer,allocatable :: NewC(:)    ! new number for cells
  integer,allocatable :: NewS(:)    ! new number for sides 
  integer,allocatable :: CelMar(:)  ! cell marker

  integer,allocatable :: NodeN2(:,:)    
  integer,allocatable :: NodeN4(:,:)  
  integer,allocatable :: NodeN8(:,:)    

  integer,allocatable :: level(:)   ! refinement level

  integer :: MAXN, MAXB, MAXS

  integer :: NR(MAXP)               ! refin. levels, refin. regions

  real    :: xp(MAXP), yp(MAXP), zp(MAXP)       ! point coordinates
  real    :: xl(MAXP,MAXL),yl(MAXP,MAXL),zl(MAXP,MAXL),LinWgt(MAXP)
  real    :: BlkWgt(MAXL,3), BlFaWt(MAXL,3)     ! leave this 
  real    :: FRegio(MAXP,MAXP,0:6)              ! levels, regions

  real    :: SRegio(MAXP,0:6), Srelax(MAXP)  ! levels, regions
  logical :: SdirX(MAXP), SdirY(MAXP), SdirZ(MAXP)   
  integer :: Siter(MAXP)   

  integer :: BlkPnt(MAXP,0:8),  & ! 0 for orientation                
             BlkRes(MAXP,6),    & ! NI,NJ,NK,NI*NJ*NK,NNo,NVo       
             BlkFac(MAXP,6,4),  &                                    
             BlFaLa(MAXP),      &                                   
             Bound(MAXP,8),     &                                  
             Period(MAXP,8),    &                                 
             Copy(MAXP,0:8)

  integer :: LinPnt(MAXL,2), LinRes(MAXL)
  integer :: Nbloc, NP, Nline, Nsurf, Nboun, Nperi
  integer :: NN, NN2, NN4, NN8
  integer :: NSR                   ! smoothing regions
  integer :: NSsh                  ! number of shadow faces

  integer :: WallFacFst, WallFacLst 

  integer :: ELIPSO, RECTAN, PLANE,YES,NO

  character*4 :: BndFac(MAXP)

end MODULE
