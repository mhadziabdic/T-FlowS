!======================================================================!
  subroutine GenLoa
!----------------------------------------------------------------------!
! Reads: NAME.d                                                        !
! ~~~~~~                                                               !
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod
  use par_mod
!----------------------------------------------------------------------! 
  implicit none
!------------------------------[Calling]-------------------------------!
  real :: TetVol   
!-------------------------------[Locals]-------------------------------!
  integer   :: b, i, l, s, fc, n, n1,n2,n3,n4
  integer   :: n_faces_check, n_nodes_check
  integer   :: ni, nj, nk
  integer   :: dum
  character :: domain_name*80
  character :: answer*12 
  real      :: xt(8), yt(8), zt(8)

  integer   :: face_nodes(6,4)
  integer   :: n_points
!======================================================================!
  data face_nodes / 1, 1, 2, 4, 3, 5,                               &
                    2, 5, 6, 8, 7, 7,                               &
                    4, 6, 8, 7, 5, 8,                               &
                    3, 2, 4, 3, 1, 6  /
!----------------------------------------------------------------------!

  write(6,'(A41)') '# Input problem name: (without extension)'
  call ReadC(5,inp,tn,ts,te) 
  read(inp, '(A80)')  name

  domain_name = name
  domain_name(len_trim(name)+1:len_trim(name)+2) = '.d'
  write(6, '(A24,A)') '# Now reading the file: ', domain_name
  open(9, FILE=domain_name)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Max. number of nodes (cells), boundary faces and cell faces     ! 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) MAXN, MAXB, MAXS  

!/////////////////////////!
!     Allocate memory     !
!/////////////////////////!

  write(6,'(A25)')       '# Allocating memory for: ' 
  write(6,'(A1,I8,A16)') '#', MAXN, ' nodes and cells' 
  write(6,'(A1,I8,A15)') '#', MAXB, ' boundary cells'         
  write(6,'(A1,I8,A11)') '#', MAXS, ' cell faces' 

!---- variables declared in all_mod.h90:
!---- (these are stored in .geo)
  allocate (xc(-MAXB:MAXN)); xc=0.0
  allocate (yc(-MAXB:MAXN)); yc=0.0
  allocate (zc(-MAXB:MAXN)); zc=0.0

  allocate (Sx(MAXS)); Sx=0.0
  allocate (Sy(MAXS)); Sy=0.0
  allocate (Sz(MAXS)); Sz=0.0

  allocate (Dx(MAXS)); Dx=0.0
  allocate (Dy(MAXS)); Dy=0.0
  allocate (Dz(MAXS)); Dz=0.0

  allocate (volume(-MAXB:MAXN)); volume=0.0
  allocate (delta(-MAXB:MAXN));  delta=0.0

  allocate (WallDs(MAXN)); WallDs=0.0
  allocate (f(MAXS)); f=0.0

!---- ()
  allocate (material(-MAXB:MAXN));  material=0
  allocate (SideC(0:2,MAXS)); SideC   =0

  allocate (CopyC(-MAXB:MAXN)); CopyC=0
  allocate (CopyS(2,MAXB));     CopyS=0    

  allocate (BCmark(-MAXB-1:-1)); BCmark=0;

!---- variables declared in gen_mod.h90:
  allocate (x_node(MAXN));     x_node=0 
  allocate (y_node(MAXN));     y_node=0
  allocate (z_node(MAXN));     z_node=0
  allocate (walln(MAXN)); walln=0
  allocate (xsp(MAXS));   xsp=0
  allocate (ysp(MAXS));   ysp=0
  allocate (zsp(MAXS));   zsp=0

  allocate (SideN(MAXS,0:4));       SideN =0 
  allocate (SideCc(MAXS,2));        SideCc=0 
  allocate (CellC(-MAXB:MAXN,24));  CellC =0

  allocate (NewN(-MAXB:MAXN));   NewN=0
  allocate (NewC(-MAXB:MAXN));   NewC=0
  allocate (NewS(MAXS));         NewS=0
  allocate (CelMar(-MAXB:MAXN)); CelMar=0

  allocate (CellN(MAXN,0:8));    CellN=0
  allocate (TwinN(MAXN,0:8));    TwinN=0

  allocate (NodeN2(MAXN,0:2));   NodeN2=0 
  allocate (NodeN4(MAXN,0:4));   NodeN4=0
  allocate (NodeN8(MAXN,0:8));   NodeN8=0

  allocate (level(MAXN)); level=0

!---- variables declared in pro_mod.h90:
  allocate (proces(MAXN)); proces=0
  allocate (BuSeIn(MAXS)); BuSeIn=0
  allocate (BuReIn(MAXS)); BuReIn=0
  allocate (BufPos(MAXS)); BufPos=0

  write(6,'(A26)') '# Allocation successfull !'

!////////////////////////!
!     Initialization     !
!////////////////////////!
  do i=1,120 
    do n=1,3
      BlkWgt(i,n) = 1.0
      BlFaWt(i,n) = 1.0
    end do
  end do

!>>>>>>>>>>>>>>>>>!
!     Corners     !
!>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) n_points  ! number of points

  allocate (x_point(n_points))
  allocate (y_point(n_points))
  allocate (z_point(n_points))

  do i = 1, n_points
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(2):te(4)),*) x_point(i), y_point(i), z_point(i)
  end do

!>>>>>>>>>>>>>>>>!
!     Blocks     !
!>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) Nbloc  ! number of blocks 

  allocate (block_points(Nbloc,0:8))
  allocate (block_resolutions(Nbloc,6))
  allocate (block_faces(Nbloc,6,4))
  allocate (face_laplace(Nbloc*6))

  do b=1,Nbloc
    block_points(b,0)=1       ! suppose it is properly oriented

    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(2):te(4)),*)  &  ! ni, nj, nk  for a block
         block_resolutions(b, 1), block_resolutions(b, 2), block_resolutions(b, 3) 
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)  &            ! Block weights 
         BlkWgt(b,1),BlkWgt(b,2),BlkWgt(b,3)
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)                       &
         block_points(b, 1), block_points(b, 2),  &
         block_points(b, 3), block_points(b, 4),  &
         block_points(b, 5), block_points(b, 6),  &
         block_points(b, 7), block_points(b, 8)

!-------------------------------!
!     Check if the block is     ! 
!       properly oriented       !
!-------------------------------!
    do n=1,8
      xt(n)=x_point(block_points(b, n))
      yt(n)=y_point(block_points(b, n))
      zt(n)=z_point(block_points(b, n))
!->>>     write(6,'(3F8.4)') xt(n),yt(n),zt(n)
    end do

    if(tetvol( xt(2),yt(2),zt(2), xt(5),yt(5),zt(5),  &
               xt(3),yt(3),zt(3), xt(1),yt(1),zt(1) )  < 0) then
      block_points(b,0)=-1            !  It's nor properly oriented
      call swapi(block_points(b,2),block_points(b,3))
      call swapi(block_points(b,6),block_points(b,7))
      call swapr(BlkWgt(b,1),BlkWgt(b,2))
      BlkWgt(b,1)=1.0/BlkWgt(b,1)
      BlkWgt(b,2)=1.0/BlkWgt(b,2)
      call swapi(block_resolutions(b,1),block_resolutions(b,2))
      write(6,*) 'Warning: Block ',b,' was not properly oriented'
    end if
  end do                 ! through blocks

!---------------------------------!
!     Set the corners of each     !
!        face of the block        !
!---------------------------------!
  do b=1,Nbloc
    do fc=1,6
      do n=1,4
        block_faces(b, fc, n)=block_points(b, face_nodes(fc,n))
      end do
    end do
  end do

!>>>>>>>>>>>>>>>!
!     Lines     !   -->> Under construction
!>>>>>>>>>>>>>>>!--------------------------------!
!     Linije mogu biti zadan tocka po tocka,     ! 
!             ili samo tezinski faktor           !
!------------------------------------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)       Nline     ! number of defined lines
  do l=1, Nline
    call ReadC(9,inp,tn,ts,te)

    read(inp(ts(1):te(3)),*) dum, LinPnt(l,1), LinPnt(l,2)

    call FinLin(LinPnt(l,1),LinPnt(l,2),LinRes(l))

    if(LinRes(l)  > MAXL) then
      write(6,*) 'ERROR MESSAGE FROM TFlowS:'
      write(6,*) 'You tried to define ', LinRes(l), ' points on'
      write(6,*) 'the line ', l, ' and the limit is: ', MAXL
      write(6,*) 'Increase the parameter MAXL in the file param.all'
      write(6,*) 'and recompile the code. Good Luck !'
      stop
    end if 

!----- zadana tocka po tocka
    if(dum  > 0) then
      do n=1,LinRes(l)
        call ReadC(9,inp,tn,ts,te)
        read(inp(ts(2):te(4)),*) xl(l,n), yl(l,n), zl(l,n)
!->>>
        write(6,*)  xl(l,n), yl(l,n), zl(l,n)
      end do
!----- zadan tezinski faktor
    else
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) LinWgt(l)
    endif 

  end do

!--------------------------------!
!                                !
!                                !
!--------------------------------!
  do b=1,Nbloc
    do fc=1,6                          !  face of the block
      n = (b-1)*6 + fc                 !  surface number
      BlFaWt(n,1)=BlkWgt(b,1)
      BlFaWt(n,2)=BlkWgt(b,2)
      BlFaWt(n,3)=BlkWgt(b,3)
      face_laplace(n)  = NO
    end do
  end do

!>>>>>>>>>>>>>>>>>>!
!     Surfaces     !  -->> Under construction
!>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)     Nsurf     ! number of defined surfaces

  do s=1,Nsurf
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) dum, n1,n2,n3,n4
    call FinSur(n1,n2,n3,n4,b,fc)
    write(6,*) 'block: ', b, ' surf: ', fc
    n = (b-1)*6 + fc         ! surface number
    face_laplace(n) = YES          ! perform Laplace

    call ReadC(9,inp,tn,ts,te)
    read(inp,*)  BlFaWt(n,1),BlFaWt(n,2),BlFaWt(n,3)
  end do

!??????????????????????????????????????????!
!     Is there enough allocated memory     !
!??????????????????????????????????????????!

!----- Nodes & Sides
  n_nodes_check = 0
  n_faces_check = 0
  do b=1,Nbloc
    ni=block_resolutions(b,1)
    nj=block_resolutions(b,2)
    nk=block_resolutions(b,3)
    n_nodes_check=n_nodes_check + ni*nj*nk
    n_faces_check=n_faces_check + ni*nj*nk + 2*( (ni*nj)+(nj*nk)+(ni*nk) )
  end do

  if( (n_faces_check  > MAXS).or.(n_nodes_check  > MAXN) ) then
    write(6,*) 'ERROR MESSAGE FROM TFlowS:'
  end if

  if( n_faces_check  > MAXS ) then
    write(6,*) 'The estimated number of sides is :', n_faces_check
    write(6,*) 'There is space available only for:', MAXS
    write(6,*) 'Increase the parameter MAXS in the input file'
    write(6,*) 'and re-run the code !'
  end if

  if( n_nodes_check  > MAXN ) then
    write(6,*) 'The estimated number of nodes is :', n_nodes_check
    write(6,*) 'There is space available only for:', MAXN
    write(6,*) 'Increase the parameter MAXN in the input file'
    write(6,*) 'and re-run the code !'
  end if 

  if( (n_faces_check  > MAXS).or.(n_nodes_check  > MAXN) ) then
    stop
  end if

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Boundary conditions and materials     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)     n_b_cond      ! number of boundary conditions 

  allocate (b_cond(n_b_cond,8))
  allocate (BndFac(n_b_cond))

  do n=1,n_b_cond
    BndFac(n)=''
    call ReadC(9,inp,tn,ts,te)
    if(tn == 7) then
      read(inp,*)  dum,                         &  
           b_cond(n,1), b_cond(n,2), b_cond(n,3),  &  ! is, js, ks 
           b_cond(n,4), b_cond(n,5), b_cond(n,6)      ! ie, je, ke
    else if(tn == 2) then
      read(inp(ts(1):te(1)),*)       dum           
      read(inp(ts(2):te(2)),'(A4)')  BndFac(n)
      call ToUppr(BndFac(n))           
    end if
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)       &  
         b_cond(n,7),  &  ! block,  
         b_cond(n,8)      ! mark         
  if( block_points(b_cond(n,7),0) == -1 ) then
    call swapi( b_cond(n,1),b_cond(n,2) )
    call swapi( b_cond(n,4),b_cond(n,5) )
  end if

  end do

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     periodic boundaries     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)  n_periodic_cond      ! number of periodic boundaries
  write(*,*) 'Number of periodic boundaries: ', n_periodic_cond 

  allocate (periodic_cond(n_periodic_cond,8))

  do n=1,n_periodic_cond
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) dum, periodic_cond(n,1), periodic_cond(n,2),  &
                     periodic_cond(n,3), periodic_cond(n,4) 
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)      periodic_cond(n,5), periodic_cond(n,6),  &
                     periodic_cond(n,7), periodic_cond(n,8)
  end do

!>>>>>>>>>>>>>>>>>>>>>>>>>!
!     copy boundaries     !
!>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)  n_copy_cond      ! number of copy boundaries
  write(*,*) 'Number of copy boundaries: ', n_copy_cond

  allocate (copy_cond(n_copy_cond,8))

  do n=1,n_copy_cond
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) dum, copy_cond(n,1), copy_cond(n,2),  &
                     copy_cond(n,3), copy_cond(n,4) 
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)      copy_cond(n,5), copy_cond(n,6),  &
                     copy_cond(n,7),copy_cond(n,8)
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) copy_cond(n,0)
  end do

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Refinement levels and regions     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) n_refine_levels ! number of refinement levels

  write(*,*) 'Number of refinement levels: ', n_refine_levels

  allocate (refined_regions(n_refine_levels, 1024, 0:6))
  allocate (n_refined_regions(n_refine_levels))

  do l=1,n_refine_levels
    write(*,*) 'Level: ', l
    call ReadC(9,inp,tn,ts,te)
    inp = inp(ts(2):te(2)) ! this is the only way Microsoft Fortran
                           ! compiles this part (two lines) of the code
    read(inp,*) n_refined_regions(l)      ! number of regions in level n

    do n=1, n_refined_regions(l)
      call ReadC(9,inp,tn,ts,te)
      read(inp(ts(3):te(3)),*) answer
      call ToUppr(answer)
      if(answer == 'RECTANGLE') then
        refined_regions(l,n,0) = RECTANGLE
      elseif(answer == 'ELIPSOID') then
        refined_regions(l,n,0) = ELIPSOID 
      elseif(answer == 'PLANE') then
        refined_regions(l,n,0) = PLANE
      else
        write(*,*) 'Error in input file: ', answer 
        stop
      endif 

      call ReadC(9,inp,tn,ts,te)
      read(inp,*)                                                   &
                  refined_regions(l,n,1),refined_regions(l,n,2),    &
                  refined_regions(l,n,3),refined_regions(l,n,4),    &
                  refined_regions(l,n,5),refined_regions(l,n,6)   
    end do
  end do

!>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Smoothing regions     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) n_smoothing_regions  ! number of smoothing regions 

  write(*,*) 'Number of (non)smoothing regions: ', n_smoothing_regions 

  allocate (smooth_in_x(n_smoothing_regions))
  allocate (smooth_in_y(n_smoothing_regions))
  allocate (smooth_in_z(n_smoothing_regions))
  allocate (smooth_iters(n_smoothing_regions))
  allocate (smooth_relax(n_smoothing_regions))
  allocate (smooth_regions(n_smoothing_regions,0:6))

  do n=1, n_smoothing_regions
    smooth_in_x(n) = .FALSE.
    smooth_in_y(n) = .FALSE.
    smooth_in_z(n) = .FALSE.
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) smooth_regions(n,0)  
    if(tn == 4) then   ! smoothing in three directions
      smooth_in_x(n) = .TRUE.
      smooth_in_y(n) = .TRUE.
      smooth_in_z(n) = .TRUE.
    else if(tn == 3) then
      call ToUppr(inp(ts(2):te(2)))
      call ToUppr(inp(ts(3):te(3)))
      if( inp(ts(2):te(2))  ==  'X' ) smooth_in_x(n) = .TRUE.
      if( inp(ts(3):te(3))  ==  'X' ) smooth_in_x(n) = .TRUE.
      if( inp(ts(2):te(2))  ==  'Y' ) smooth_in_y(n) = .TRUE.
      if( inp(ts(3):te(3))  ==  'Y' ) smooth_in_y(n) = .TRUE.
      if( inp(ts(2):te(2))  ==  'Z' ) smooth_in_z(n) = .TRUE.
      if( inp(ts(3):te(3))  ==  'Z' ) smooth_in_z(n) = .TRUE.
    else if(tn == 2) then
      call ToUppr(inp(ts(2):te(2)))
      if( inp(ts(2):te(2))  ==  'X' ) smooth_in_x(n) = .TRUE.
      if( inp(ts(2):te(2))  ==  'Y' ) smooth_in_y(n) = .TRUE.
      if( inp(ts(2):te(2))  ==  'Z' ) smooth_in_z(n) = .TRUE.
    end if 

!---- read the coordinates of the (non)smoothed region
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) smooth_iters(n), smooth_relax(n)
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) smooth_regions(n,1),smooth_regions(n,2),smooth_regions(n,3), &
                smooth_regions(n,4),smooth_regions(n,5),smooth_regions(n,6)   
  end do

  close(9)

  end subroutine GenLoa
