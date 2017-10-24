!==============================================================================!
  subroutine Load_Domain
!------------------------------------------------------------------------------!
!   Reads: name.d                                                              !
!----------------------------------[Modules]-----------------------------------!
  use Block_Mod
  use all_mod
  use gen_mod
  use par_mod
!------------------------------------------------------------------------------! 
  implicit none
!----------------------------------[Calling]-----------------------------------!
  real :: Tet_Volume   
!-----------------------------------[Locals]-----------------------------------!
  integer           :: b, i, l, s, fc, n, n1,n2,n3,n4
  integer           :: n_faces_check, n_nodes_check
  integer           :: ni, nj, nk
  integer           :: dum
  character(len=80) :: domain_name
  character(len=12) :: answer
  real              :: xt(8), yt(8), zt(8)
  integer           :: face_nodes(6,4)
  integer           :: n_points
!==============================================================================!
  data face_nodes / 1, 1, 2, 4, 3, 5,                               &
                    2, 5, 6, 8, 7, 7,                               &
                    4, 6, 8, 7, 5, 8,                               &
                    3, 2, 4, 3, 1, 6  /
!------------------------------------------------------------------------------!

  write(*,*) '#========================================'
  write(*,*) '# Input problem name: (without extension)'
  write(*,*) '#----------------------------------------'
  call ReadC(5,inp,tn,ts,te) 
  read(inp, '(A80)')  name

  domain_name = name
  domain_name(len_trim(name)+1:len_trim(name)+2) = '.d'
  write(*, *) '# Now reading the file: ', domain_name
  open(9, file=domain_name)

  !-----------------------------------------------------------------!
  !   Max. number of nodes (cells), boundary faces and cell faces   ! 
  !-----------------------------------------------------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) MAXN, MAXB, MAXS  

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  write(*,*) '# Allocating memory for: ' 
  write(*,*) '#', MAXN, ' nodes and cells' 
  write(*,*) '#', MAXB, ' boundary cells'         
  write(*,*) '#', MAXS, ' cell faces' 

  ! Variables declared in all_mod.h90:
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

  allocate (material(-MAXB:MAXN));  material=0
  allocate (SideC(0:2,MAXS)); SideC   =0

  allocate (CopyC(-MAXB:MAXN)); CopyC=0
  allocate (CopyS(2,MAXB));     CopyS=0    

  allocate (BCmark(-MAXB-1:-1)); BCmark=0;

  ! Variables declared in gen_mod.h90:
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

  ! Variables declared in pro_mod.h90:
  allocate (proces(MAXN)); proces=0
  allocate (BuSeIn(MAXS)); BuSeIn=0
  allocate (BuReIn(MAXS)); BuReIn=0
  allocate (BufPos(MAXS)); BufPos=0

  write(*,*) '# Allocation successfull !'

  !--------------------!
  !   Initialization   !
  !--------------------!
  do i=1,120 
    do n=1,3
      BlkWgt(i,n) = 1.0
      BlFaWt(i,n) = 1.0
    end do
  end do

  !-------------!
  !   Corners   !
  !-------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) n_points  ! number of points

  allocate (x_point(n_points))
  allocate (y_point(n_points))
  allocate (z_point(n_points))

  do i = 1, n_points
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(2):te(4)),*) x_point(i), y_point(i), z_point(i)
  end do

  !------------!
  !   Blocks   !
  !------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) Nbloc  ! number of blocks 

  allocate (blocks(Nbloc))

  do b=1,Nbloc
    blocks(b) % points(0)=1       ! suppose it is properly oriented

    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(2):te(4)),*)          &  ! ni, nj, nk  for a block
         blocks(b) % resolutions(1),  &
         blocks(b) % resolutions(2),  &
         blocks(b) % resolutions(3) 
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)  &            ! Block weights 
         BlkWgt(b,1),BlkWgt(b,2),BlkWgt(b,3)
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)                       &
         blocks(b) % points( 1), blocks(b) % points( 2),  &
         blocks(b) % points( 3), blocks(b) % points( 4),  &
         blocks(b) % points( 5), blocks(b) % points( 6),  &
         blocks(b) % points( 7), blocks(b) % points( 8)

    !---------------------------!
    !   Check if the block is   ! 
    !     properly oriented     !
    !---------------------------!
    do n=1,8
      xt(n)=x_point(blocks(b) % points( n))
      yt(n)=y_point(blocks(b) % points( n))
      zt(n)=z_point(blocks(b) % points( n))
    end do

    if(Tet_Volume( xt(2),yt(2),zt(2), xt(5),yt(5),zt(5),  &
                   xt(3),yt(3),zt(3), xt(1),yt(1),zt(1) )  < 0) then
      blocks(b) % points(0)=-1            !  It's nor properly oriented
      call Swap_Integers(blocks(b) % points(2),blocks(b) % points(3))
      call Swap_Integers(blocks(b) % points(6),blocks(b) % points(7))
      call Swap_Reals(BlkWgt(b,1),BlkWgt(b,2))
      BlkWgt(b,1)=1.0/BlkWgt(b,1)
      BlkWgt(b,2)=1.0/BlkWgt(b,2)
      call Swap_Integers(blocks(b) % resolutions(1),blocks(b) % resolutions(2))
      write(*,*) 'Warning: Block ',b,' was not properly oriented'
    end if
  end do                 ! through blocks

  !-----------------------------!
  !   Set the corners of each   !
  !      face of the block      !
  !-----------------------------!
  do b=1,Nbloc
    do fc=1,6
      do n=1,4
        blocks(b) % faces(fc, n)=blocks(b) % points(face_nodes(fc,n))
      end do
    end do
  end do

  !----------------------------------------------!
  !   Lines                                      !
  !----------------------------------------------!
  !   Lines can be prescribed point by point     ! 
  !   or with just a weighting factor.           !
  !----------------------------------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)       Nline     ! number of defined lines
  do l=1, Nline
    call ReadC(9,inp,tn,ts,te)

    read(inp(ts(1):te(3)),*) dum, LinPnt(l,1), LinPnt(l,2)

    call Find_Line(LinPnt(l,1),LinPnt(l,2),LinRes(l))

    if(LinRes(l)  > MAXL) then
      write(*,*) 'ERROR MESSAGE FROM TFlowS:'
      write(*,*) 'You tried to define ', LinRes(l), ' points on'
      write(*,*) 'the line ', l, ' and the limit is: ', MAXL
      write(*,*) 'Increase the parameter MAXL in the file param.all'
      write(*,*) 'and recompile the code. Good Luck !'
      stop
    end if 

    ! Point by point
    if(dum  > 0) then
      do n=1,LinRes(l)
        call ReadC(9,inp,tn,ts,te)
        read(inp(ts(2):te(4)),*) xl(l,n), yl(l,n), zl(l,n)
        write(*,*)  xl(l,n), yl(l,n), zl(l,n)
      end do

    ! Weight factor
    else
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) LinWgt(l)
    endif 

  end do

  !----------------------------------------!
  !   Copy block weights to face weights   !
  !----------------------------------------!
  do b=1,Nbloc
    do fc=1,6                          !  face of the block
      n = (b-1)*6 + fc                 !  surface number
      BlFaWt(n,1)=BlkWgt(b,1)
      BlFaWt(n,2)=BlkWgt(b,2)
      BlFaWt(n,3)=BlkWgt(b,3)
    end do
  end do

  !--------------!
  !   Surfaces   !
  !--------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)     Nsurf     ! number of defined surfaces

  do s=1,Nsurf
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) dum, n1,n2,n3,n4
    call Find_Surface(n1,n2,n3,n4,b,fc)
    write(*,*) 'block: ', b, ' surf: ', fc
    n = (b-1)*6 + fc         ! surface number

    call ReadC(9,inp,tn,ts,te)
    read(inp,*)  BlFaWt(n,1),BlFaWt(n,2),BlFaWt(n,3)
  end do

  !---------------------------------------!
  !   Is there enough allocated memory?   !
  !---------------------------------------!

  ! Nodes & faces
  n_nodes_check = 0
  n_faces_check = 0
  do b=1,Nbloc
    ni = blocks(b) % resolutions(1)
    nj = blocks(b) % resolutions(2)
    nk = blocks(b) % resolutions(3)
    n_nodes_check=n_nodes_check + ni*nj*nk
    n_faces_check=n_faces_check + ni*nj*nk + 2*( (ni*nj)+(nj*nk)+(ni*nk) )
  end do

  if( (n_faces_check  > MAXS).or.(n_nodes_check  > MAXN) ) then
    write(*,*) '# Error message from TFlowS:'
  end if

  if( n_faces_check  > MAXS ) then
    write(*,*) '# The estimated number of sides is :', n_faces_check
    write(*,*) '# There is space available only for:', MAXS
    write(*,*) '# Increase the parameter MAXS in the input file'
    write(*,*) '# and re-run the code !'
  end if

  if( n_nodes_check  > MAXN ) then
    write(*,*) '# The estimated number of nodes is :', n_nodes_check
    write(*,*) '# There is space available only for:', MAXN
    write(*,*) '# Increase the parameter MAXN in the input file'
    write(*,*) '# and re-run the code !'
  end if 

  if( (n_faces_check  > MAXS).or.(n_nodes_check  > MAXN) ) then
    stop
  end if

  !---------------------------------------!
  !   Boundary conditions and materials   !
  !---------------------------------------!
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
      call To_Upper_Case(BndFac(n))           
    end if
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)       &  
         b_cond(n,7),  &  ! block,  
         b_cond(n,8)      ! mark         
  if( blocks(b_cond(n,7)) % points(0) == -1 ) then
    call Swap_Integers( b_cond(n,1),b_cond(n,2) )
    call Swap_Integers( b_cond(n,4),b_cond(n,5) )
  end if

  end do

  !-------------------------!
  !   Periodic boundaries   !
  !-------------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)  n_periodic_cond      ! number of periodic boundaries
  write(*,*) '# Number of periodic boundaries: ', n_periodic_cond 

  allocate (periodic_cond(n_periodic_cond,8))

  do n=1,n_periodic_cond
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) dum, periodic_cond(n,1), periodic_cond(n,2),  &
                     periodic_cond(n,3), periodic_cond(n,4) 
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)      periodic_cond(n,5), periodic_cond(n,6),  &
                     periodic_cond(n,7), periodic_cond(n,8)
  end do

  !---------------------!
  !   Copy boundaries   !
  !---------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*)  n_copy_cond      ! number of copy boundaries
  write(*,*) '# Number of copy boundaries: ', n_copy_cond

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

  !-----------------------------------!
  !   Refinement levels and regions   !
  !-----------------------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) n_refine_levels ! number of refinement levels

  write(*,*) '# Number of refinement levels: ', n_refine_levels

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
      call To_Upper_Case(answer)
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

  !-----------------------!
  !   Smoothing regions   !
  !-----------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) n_smoothing_regions  ! number of smoothing regions 

  write(*,*) '# Number of (non)smoothing regions: ', n_smoothing_regions 

  allocate (smooth_in_x(n_smoothing_regions))
  allocate (smooth_in_y(n_smoothing_regions))
  allocate (smooth_in_z(n_smoothing_regions))
  allocate (smooth_iters(n_smoothing_regions))
  allocate (smooth_relax(n_smoothing_regions))
  allocate (smooth_regions(n_smoothing_regions,0:6))

  do n=1, n_smoothing_regions
    smooth_in_x(n) = .false.
    smooth_in_y(n) = .false.
    smooth_in_z(n) = .false.
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) smooth_regions(n,0)  
    if(tn == 4) then   ! smoothing in three directions
      smooth_in_x(n) = .true.
      smooth_in_y(n) = .true.
      smooth_in_z(n) = .true.
    else if(tn == 3) then
      call To_Upper_Case(inp(ts(2):te(2)))
      call To_Upper_Case(inp(ts(3):te(3)))
      if( inp(ts(2):te(2))  ==  'X' ) smooth_in_x(n) = .true.
      if( inp(ts(3):te(3))  ==  'X' ) smooth_in_x(n) = .true.
      if( inp(ts(2):te(2))  ==  'Y' ) smooth_in_y(n) = .true.
      if( inp(ts(3):te(3))  ==  'Y' ) smooth_in_y(n) = .true.
      if( inp(ts(2):te(2))  ==  'Z' ) smooth_in_z(n) = .true.
      if( inp(ts(3):te(3))  ==  'Z' ) smooth_in_z(n) = .true.
    else if(tn == 2) then
      call To_Upper_Case(inp(ts(2):te(2)))
      if( inp(ts(2):te(2))  ==  'X' ) smooth_in_x(n) = .true.
      if( inp(ts(2):te(2))  ==  'Y' ) smooth_in_y(n) = .true.
      if( inp(ts(2):te(2))  ==  'Z' ) smooth_in_z(n) = .true.
    end if 

    ! Read the coordinates of the (non)smoothed region
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) smooth_iters(n), smooth_relax(n)
    call ReadC(9,inp,tn,ts,te)
    read(inp,*) smooth_regions(n,1),smooth_regions(n,2),smooth_regions(n,3), &
                smooth_regions(n,4),smooth_regions(n,5),smooth_regions(n,6)   
  end do

  close(9)

  end subroutine Load_Domain
