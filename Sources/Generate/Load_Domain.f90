!==============================================================================!
  subroutine Load_Domain(dom, grid)
!------------------------------------------------------------------------------!
!   Reads: name.d                                                              !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use par_mod
  use Domain_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Domain_Type) :: dom
  type(Grid_Type)   :: grid
!----------------------------------[Calling]-----------------------------------!
  real :: Tet_Volume   
!-----------------------------------[Locals]-----------------------------------!
  integer           :: b, i, l, s, fc, n, n1,n2,n3,n4
  integer           :: n_faces_check, n_nodes_check
  integer           :: ni, nj, nk, npnt
  character(len=12) :: dum
  character(len=80) :: domain_name
  character(len=12) :: answer
  real              :: xt(8), yt(8), zt(8)
  integer           :: face_nodes(6,4)
!==============================================================================!
  data face_nodes / 1, 1, 2, 4, 3, 5,         &
                    2, 5, 6, 8, 7, 7,         &
                    4, 6, 8, 7, 5, 8,         &
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
  read(inp,*) grid % max_n_nodes,           &
              grid % max_n_bnd_cells,  &
              grid % max_n_faces  

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  write(*,*) '# Allocating memory for: ' 
  write(*,*) '#', grid % max_n_nodes, ' nodes and cells' 
  write(*,*) '#', grid % max_n_bnd_cells, ' boundary cells'         
  write(*,*) '#', grid % max_n_faces, ' cell faces' 

  ! Variables declared in all_mod.h90:
  allocate (delta (-grid % max_n_bnd_cells:grid % max_n_nodes))
  delta=0.0

  allocate (WallDs(grid % max_n_nodes))
  WallDs=0.0
  allocate (f(grid % max_n_faces))
  f=0.0

  allocate (material(-grid % max_n_bnd_cells:grid % max_n_nodes)) 
  material=0
  allocate (SideC(0:2,grid % max_n_faces))
  SideC = 0

  allocate (CopyC(-grid % max_n_bnd_cells:grid % max_n_nodes))
  CopyC=0
  allocate (CopyS(2,grid % max_n_bnd_cells))
  CopyS=0    

  allocate (BCmark(-grid % max_n_bnd_cells-1:-1))
  BCmark=0;

  ! Variables in Grid_Mod
  call Grid_Mod_Allocate_Nodes(grid,  &
                               grid % max_n_nodes) 
  call Grid_Mod_Allocate_Cells(grid,                         &
                               grid % max_n_bnd_cells,  &
                               grid % max_n_nodes) 
  call Grid_Mod_Allocate_Faces(grid,  &
                               grid % max_n_faces) 

  ! Variables declared in gen_mod.h90:
  allocate (walln(grid % max_n_nodes))
  walln=0

  allocate (SideCc(grid % max_n_faces,2))
  SideCc=0 

  ! Variables for renumbering
  allocate (NewN(-grid % max_n_bnd_cells:grid % max_n_nodes))
  NewN=0
  allocate (NewC(-grid % max_n_bnd_cells:grid % max_n_nodes))
  NewC=0
  allocate (NewS( grid % max_n_faces))
  NewS=0
  allocate (CelMar(-grid % max_n_bnd_cells:grid % max_n_nodes))
  CelMar=0

  allocate (TwinN(grid % max_n_nodes,0:8))
  TwinN=0

  allocate (level(grid % max_n_nodes))
  level=0

  ! Variables declared in pro_mod.h90:
  allocate (proces(grid % max_n_nodes))
  proces=0
  allocate (BuSeIn(grid % max_n_faces))
  BuSeIn=0
  allocate (BuReIn(grid % max_n_faces))
  BuReIn=0
  allocate (BufPos(grid % max_n_faces))
  BufPos=0

  write(*,*) '# Allocation successfull !'

  !-------------!
  !   Corners   !
  !-------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) dom % n_points  ! number of points

  call Domain_Mod_Allocate_Points(dom, dom % n_points)

  do i = 1, dom % n_points
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(2):te(4)),*) dom % points(i) % x,  &
                             dom % points(i) % y,  &
                             dom % points(i) % z
  end do

  !------------!
  !   Blocks   !
  !------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) dom % n_blocks  ! number of blocks 

  call Domain_Mod_Allocate_Blocks(dom, dom % n_blocks)

  ! Initialize weights 
  do b=1, dom % n_blocks
    dom % blocks(b) % weights      = 1.0
    dom % blocks(b) % face_weights = 1.0
  end do

  do b = 1, dom % n_blocks
    dom % blocks(b) % corners(0)=1       ! suppose it is properly oriented

    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(2):te(4)),*)                &  ! ni, nj, nk  for a block
         dom % blocks(b) % resolutions(1),  &
         dom % blocks(b) % resolutions(2),  &
         dom % blocks(b) % resolutions(3) 
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)                         &  ! Block weights 
         dom % blocks(b) % weights(1),  &
         dom % blocks(b) % weights(2),  &
         dom % blocks(b) % weights(3)
    call ReadC(9,inp,tn,ts,te)
    read(inp,*)                                                       &
         dom % blocks(b) % corners(1), dom % blocks(b) % corners(2),  &
         dom % blocks(b) % corners(3), dom % blocks(b) % corners(4),  &
         dom % blocks(b) % corners(5), dom % blocks(b) % corners(6),  &
         dom % blocks(b) % corners(7), dom % blocks(b) % corners(8)

    !---------------------------!
    !   Check if the block is   ! 
    !     properly oriented     !
    !---------------------------!
    do n=1,8
      xt(n) = dom % points(dom % blocks(b) % corners(n)) % x
      yt(n) = dom % points(dom % blocks(b) % corners(n)) % y
      zt(n) = dom % points(dom % blocks(b) % corners(n)) % z
    end do

    if(Tet_Volume( xt(2),yt(2),zt(2), xt(5),yt(5),zt(5),  &
                   xt(3),yt(3),zt(3), xt(1),yt(1),zt(1) )  < 0) then
      dom % blocks(b) % corners(0)=-1            !  It's nor properly oriented
      call Swap_Integers(dom % blocks(b) % corners(2),  &
                         dom % blocks(b) % corners(3))
      call Swap_Integers(dom % blocks(b) % corners(6),  &
                         dom % blocks(b) % corners(7))
      call Swap_Reals(dom % blocks(b) % weights(1),  &
                      dom % blocks(b) % weights(2))
      dom % blocks(b) % weights(1) = 1.0 / dom % blocks(b) % weights(1)
      dom % blocks(b) % weights(2) = 1.0 / dom % blocks(b) % weights(2)
      call Swap_Integers(dom % blocks(b) % resolutions(1),  &
                         dom % blocks(b) % resolutions(2))
      write(*,*) 'Warning: Block ',b,' was not properly oriented'
    end if
  end do                 ! through dom % blocks

  !-----------------------------!
  !   Set the corners of each   !
  !      face of the block      !
  !-----------------------------!
  do b = 1, dom % n_blocks
    do fc = 1,6
      do n = 1,4
        dom % blocks(b) % faces(fc, n) =  &
        dom % blocks(b) % corners(face_nodes(fc,n))
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
  read(inp,*) dom % n_lines     ! number of defined dom % lines

  call Domain_Mod_Allocate_Lines(dom, dom % n_lines)

  do l=1, dom % n_lines
    call ReadC(9,inp,tn,ts,te)

    read(inp(ts(1):te(3)),*) npnt, dom % lines(l) % points(1),  &
                                   dom % lines(l) % points(2)

    call Find_Line(dom,                         &
                   dom % lines(l) % points(1),  &
                   dom % lines(l) % points(2),  &
                   dom % lines(l) % resolution)

    ! Does this need a more elegant solution?
    allocate(dom % lines(l) % x( dom % lines(l) % resolution ))
    allocate(dom % lines(l) % y( dom % lines(l) % resolution ))
    allocate(dom % lines(l) % z( dom % lines(l) % resolution ))

    ! Point by point
    if(npnt > 0) then
      do n=1,dom % lines(l) % resolution
        call ReadC(9,inp,tn,ts,te)
        read(inp(ts(2):te(4)),*) dom % lines(l) % x(n),  &
                                 dom % lines(l) % y(n),  &
                                 dom % lines(l) % z(n)
      end do

    ! Weight factor
    else
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) dom % lines(l) % weight
    endif 

  end do

  !----------------------------------------!
  !   Copy block weights to face weights   !
  !----------------------------------------!
  do b = 1, dom % n_blocks
    do fc = 1,6                          !  face of the block
      dom % blocks(b) % face_weights(fc, 1) = dom % blocks(b) % weights(1)
      dom % blocks(b) % face_weights(fc, 2) = dom % blocks(b) % weights(2)
      dom % blocks(b) % face_weights(fc, 3) = dom % blocks(b) % weights(3)
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
    call Find_Surface(dom, n1, n2, n3, n4, b, fc)
    write(*,*) '# block: ', b, ' surf: ', fc
    n = (b-1)*6 + fc         ! surface number

    call ReadC(9,inp,tn,ts,te)
    read(inp,*) dom % blocks(b) % face_weights(fc,1), &
                dom % blocks(b) % face_weights(fc,2), &
                dom % blocks(b) % face_weights(fc,2)
  end do

  !---------------------------------------!
  !   Is there enough allocated memory?   !
  !---------------------------------------!

  ! Nodes & faces
  n_nodes_check = 0
  n_faces_check = 0
  do b = 1, dom % n_blocks
    ni = dom % blocks(b) % resolutions(1)
    nj = dom % blocks(b) % resolutions(2)
    nk = dom % blocks(b) % resolutions(3)
    n_nodes_check=n_nodes_check + ni*nj*nk
    n_faces_check=n_faces_check + ni*nj*nk + 2*( (ni*nj)+(nj*nk)+(ni*nk) )
  end do

  if( (n_faces_check > grid % max_n_faces) .or.  &
      (n_nodes_check > grid % max_n_nodes) ) then
    write(*,*) '# Error message from TFlowS:'
  end if

  if( n_faces_check  > grid % max_n_faces ) then
    write(*,*) '# The estimated number of sides is :', n_faces_check
    write(*,*) '# There is space available only for:', grid % max_n_faces
    write(*,*) '# Increase the parameter grid % max_n_faces in the input file'
    write(*,*) '# and re-run the code !'
  end if

  if( n_nodes_check  > grid % max_n_nodes ) then
    write(*,*) '# The estimated number of nodes is :', n_nodes_check
    write(*,*) '# There is space available only for:', grid % max_n_nodes
    write(*,*) '# Increase the parameter grid % max_n_nodes in the input file'
    write(*,*) '# and re-run the code !'
  end if 

  if( (n_faces_check > grid % max_n_faces) .or.  &
      (n_nodes_check > grid % max_n_nodes) ) then
    stop
  end if

  !---------------------------------------!
  !   Boundary conditions and materials   !
  !---------------------------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) dom % n_regions  ! number of regions (can be boundary 
                               ! conditions or materials)

  call Domain_Mod_Allocate_Regions(dom, dom % n_regions)

  do n=1, dom % n_regions
    dom % regions(n) % face=''

    call ReadC(9,inp,tn,ts,te)
    if(tn == 7) then
      read(inp,*)  dum,                     &  
                   dom % regions(n) % is,   &
                   dom % regions(n) % js,   &
                   dom % regions(n) % ks,   & 
                   dom % regions(n) % ie,   &
                   dom % regions(n) % je,   &
                   dom % regions(n) % ke   
    else if(tn == 2) then
      read(inp(ts(1):te(1)),*)       dum           
      read(inp(ts(2):te(2)),'(A4)')  dom % regions(n) % face
      call To_Upper_Case(dom % regions(n) % face)           
    end if

    call ReadC(9,inp,tn,ts,te)
    read(inp,*) dom % regions(n) % block,  & 
                dom % regions(n) % name    
    call To_Upper_Case(dom % regions(n) % name)           

    ! if( dom % blocks(b_cond(n,7)) % points(0) == -1 ) then
    !   call Swap_Integers( dom % regions(n) % is,dom % regions(n) % js )
    !   call Swap_Integers( dom % regions(n) % ie,dom % regions(n) % je )
    ! end if

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
    inp = inp(ts(2):te(2))  ! this is the only way Microsoft Fortran
                            ! compiles this part (two lines) of the code
    read(inp,*) n_refined_regions(l)  ! number of regions in level n

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

  end subroutine
