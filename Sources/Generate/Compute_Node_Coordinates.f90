!==============================================================================!
  subroutine Compute_Node_Coordinates
!------------------------------------------------------------------------------!
!   Calculate node coordinates inside the domain, block by block.              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Block_Mod
!------------------------------------------------------------------------------! 
  implicit none
!----------------------------------[Calling]-----------------------------------!
  integer :: Is_Line_in_Block
!---------------------------------[Interface]----------------------------------!
  include "../Shared/Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer :: fc, b, bl, i, j, k, n, c, ig
  integer :: l, l1, l2
  integer :: is, js, ks, ie, je, ke, face 
  integer :: ni, nj, nk, ci, cj, ck
  integer :: trans(3,2)
!==============================================================================!

  do n=1,MAXN
    x_node(n) = HUGE
  end do 

  !------------------------------------!
  !   Calculate the node coordinates   !
  !------------------------------------!
  nn = 0  ! initialize n.o.n.
  nc = 0  ! initialize n.o.v.

  do b = 1, size(blocks)

    write(*,*) '# Generating block: ', b
    ni=blocks(b) % resolutions( 1)
    nj=blocks(b) % resolutions( 2)
    nk=blocks(b) % resolutions( 3)   

    ! ( 1 )
    n = nn + ( 1-1)*ni*nj + ( 1-1)*ni +  1
    x_node(n) = x_point(blocks(b) % points(1)) 
    y_node(n) = y_point(blocks(b) % points(1)) 
    z_node(n) = z_point(blocks(b) % points(1)) 

    ! ( 2 )
    n = nn + ( 1-1)*ni*nj + ( 1-1)*ni + ni
    x_node(n) = x_point(blocks(b) % points(2)) 
    y_node(n) = y_point(blocks(b) % points(2)) 
    z_node(n) = z_point(blocks(b) % points(2)) 

    ! ( 3 )
    n = nn + ( 1-1)*ni*nj + (nj-1)*ni +  1
    x_node(n) = x_point(blocks(b) % points(3)) 
    y_node(n) = y_point(blocks(b) % points(3)) 
    z_node(n) = z_point(blocks(b) % points(3)) 

    ! ( 4 )
    n = nn + ( 1-1)*ni*nj + (nj-1)*ni + ni
    x_node(n) = x_point(blocks(b) % points(4)) 
    y_node(n) = y_point(blocks(b) % points(4)) 
    z_node(n) = z_point(blocks(b) % points(4)) 

    ! ( 5 ) !
    n = nn + (nk-1)*ni*nj + ( 1-1)*ni +  1
    x_node(n) = x_point(blocks(b) % points(5)) 
    y_node(n) = y_point(blocks(b) % points(5)) 
    z_node(n) = z_point(blocks(b) % points(5)) 

    ! ( 6 ) !
    n = nn + (nk-1)*ni*nj + ( 1-1)*ni + ni
    x_node(n) = x_point(blocks(b) % points(6)) 
    y_node(n) = y_point(blocks(b) % points(6)) 
    z_node(n) = z_point(blocks(b) % points(6)) 

    ! ( 7 ) !
    n = nn + (nk-1)*ni*nj + (nj-1)*ni +  1
    x_node(n) = x_point(blocks(b) % points(7)) 
    y_node(n) = y_point(blocks(b) % points(7)) 
    z_node(n) = z_point(blocks(b) % points(7)) 

    ! ( 8 ) !
    n = nn + (nk-1)*ni*nj + (nj-1)*ni + ni
    x_node(n) = x_point(blocks(b) % points(8)) 
    y_node(n) = y_point(blocks(b) % points(8)) 
    z_node(n) = z_point(blocks(b) % points(8)) 

    !----------------------------!
    !     First on the lines     !
    !   defined point by point   !
    !----------------------------!
    do l=1, Nline    

      bl = Is_Line_in_Block(LinPnt(l,1),LinPnt(l,2),b)

      if(bl == b) then

        do n=1,8
          if( LinPnt(l,1) == blocks(b) % points(n)  ) l1=n
          if( LinPnt(l,2) == blocks(b) % points(n)  ) l2=n
        end do

        ! Line is defined in the +i direction
        if     ( (l2-l1) == +1 ) then
          trans(1,1) = 0
          trans(1,2) =+1
          trans(2,2) = 0
          trans(3,2) = 0
          if( (l1 == 1).or.(l1 == 5) ) then 
            trans(2,1) =1
          else                         
            trans(2,1) =blocks(b) % resolutions(2)
          endif
          if( (l1 == 1).or.(l1 == 3) ) then
            trans(3,1) =1
          else                         
            trans(3,1) =blocks(b) % resolutions(3)
          endif

        ! Line is defined in the -i direction
        else if( (l2-l1) == -1 ) then
          trans(1,1) =blocks(b) % resolutions(1)+1   ! ni from the block + 1
          trans(1,2) =-1
          trans(2,2) = 0
          trans(3,2) = 0
          if( (l1 == 2).or.(l1 == 6) ) then
            trans(2,1) =1
          else                         
            trans(2,1) =blocks(b) % resolutions(2)
          endif
          if( (l1 == 2).or.(l1 == 4) ) then
            trans(3,1) =1
          else                         
            trans(3,1) =blocks(b) % resolutions(3)
          endif

        ! Line is defined in the +j direction
        else if( (l2-l1) == +2 ) then
          trans(2,1) = 0
          trans(2,2) =+1
          trans(1,2) = 0
          trans(3,2) = 0 
          if( (l1 == 1).or.(l1 == 5) ) then
            trans(1,1) =1
          else                         
            trans(1,1) =blocks(b) % resolutions(1)
          endif
          if( (l1 == 1).or.(l1 == 2) ) then
            trans(3,1) =1
          else                         
            trans(3,1) =blocks(b) % resolutions(3)
          endif

        ! Line is defined in the -j direction
        else if( (l2-l1) == -2 ) then
          trans(2,1) =blocks(b) % resolutions(2)+1   ! nj from the block + 1
          trans(2,2) =-1
          trans(1,2) = 0
          trans(3,2) = 0 
          if( (l1 == 3).or.(l1 == 7) ) then
            trans(1,1) =1
          else                         
            trans(1,1) =blocks(b) % resolutions(1)
          endif
          if( (l1 == 3).or.(l1 == 4) ) then
            trans(3,1) =1
          else                         
            trans(3,1) =blocks(b) % resolutions(3)
          endif

        ! Line is defined in the +k direction
        else if( (l2-l1) == +4 ) then
          trans(3,1) = 0
          trans(3,2) =+1
          trans(1,2) = 0
          trans(2,2) = 0 
          if( (l1 == 1).or.(l1 == 3) ) then
            trans(1,1) =1
          else                         
            trans(1,1) =blocks(b) % resolutions(1)
          end if
          if( (l1 == 1).or.(l1 == 2) ) then
            trans(2,1) =1
          else                         
            trans(2,1) =blocks(b) % resolutions(2)
          endif

        ! Line is defined in the -k direction
        else if( (l2-l1) == -4 ) then
          trans(3,1) =blocks(b) % resolutions(3) + 1  ! nk from the block + 1
          trans(3,2) =-1
          trans(1,2) = 0
          trans(2,2) = 0 
          if( (l1 == 5).or.(l1 == 7) ) then
            trans(1,1) =1
          else                         
            trans(1,1) =blocks(b) % resolutions(1)
          endif
          if( (l1 == 5).or.(l1 == 6) ) then
            trans(2,1) =1
          else                         
            trans(2,1) =blocks(b) % resolutions(2)
          endif

        endif ! l1-l2

        ! Line is defined point by point
        if(LinWgt(l) ==  0.0) then
          write(*,*) '# Line: ', l
          write(*,*) '# l1= ', l1
          write(*,*) '# l2= ', l2
          do ig=1,LinRes(l)
            i=trans(1,1)+trans(1,2)*ig
            j=trans(2,1)+trans(2,2)*ig
            k=trans(3,1)+trans(3,2)*ig

            n = nn + (k-1)*ni*nj + (j-1)*ni + i
            x_node(n) = xl(l,ig)
            y_node(n) = yl(l,ig)
            z_node(n) = zl(l,ig)
          end do

        ! Line is defined with a weight factor
        else
          is=trans(1,1)+trans(1,2)
          js=trans(2,1)+trans(2,2)
          ks=trans(3,1)+trans(3,2)
          ie=trans(1,1)+trans(1,2)*LinRes(l)
          je=trans(2,1)+trans(2,2)*LinRes(l)
          ke=trans(3,1)+trans(3,2)*LinRes(l)
          call Distribute_Nodes(b, LinWgt(l), is, js, ks, ie, je, ke)
        endif  

      endif ! if the block contains

    end do ! for the lines

    !-----------!
    !   Lines   !
    !-----------!
    do k=1,nk,nk-1
      do j=1,nj,nj-1
        call Distribute_Nodes(b, blocks(b) % weights(1), 1,j,k,ni,j,k)
      end do
    end do

    do k=1,nk,nk-1
      do i=1,ni,ni-1
        call Distribute_Nodes(b, blocks(b) % weights(2), i,1,k,i,nj,k)
      end do
    end do

    do j=1,nj,nj-1
      do i=1,ni,ni-1
        call Distribute_Nodes(b, blocks(b) % weights(3), i,j,1,i,j,nk)
      end do
    end do

    !------------------------------------------------------------!
    !   Surfaces...                                              !
    !                                                            !
    !   I think this is the propper way to calculate surfaces:   !
    !   it spans the lines in the direction of higher weigh      !
    !------------------------------------------------------------!

    ! I (k=1) 
    fc = 1   ! face index
    k = 1
    if( .not. Approx(blocks(b) % face_weights(fc,1),1.0 ) ) then
      do j=1,nj
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,1), 1,j,k,ni,j,k)
      end do
    else ! lines in the j direction
      do i=1,ni
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,2), i,1,k,i,nj,k)
      end do
    endif

    ! VI (k=nk)
    fc = 6   ! face index
    k = nk
    if( .not. Approx(blocks(b) % face_weights(fc,1),1.0 ) ) then
     do j=1,nj
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,1), 1,j,k,ni,j,k)
      end do
    else ! lines in the j direction
      do i=1,ni
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,2), i,1,k,i,nj,k)
      end do
    endif

    ! V (i=1)
    fc = 5   ! face index
    i = 1
    if( .not. Approx(blocks(b) % face_weights(fc,3),1.0 ) ) then
      do j=1,nj
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,3), i,j,1,i,j,nk)
      end do
    else ! lines in the j direction
      do k=1,nk
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,2), i,1,k,i,nj,k)
      end do
    end if 

    ! III (i=ni)
    fc = 3   ! face index
    i = ni
    if( .not. Approx(blocks(b) % face_weights(fc,3),1.0 ) ) then
      do j=1,nj
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,3), i,j,1,i,j,nk)
      end do
    else ! lines in the j direction
      do k=1,nk
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,2), i,1,k,i,nj,k)
      end do
    end if 

    ! II (j=1)       
    fc = 2   ! face index
    j = 1
    if( .not. Approx(blocks(b) % face_weights(fc,3),1.0 ) ) then
      do i=1,ni
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,3), i,j,1,i,j,nk)
      end do
    else ! lines in the i direction
      do k=1,nk
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,1), 1,j,k,ni,j,k)
      end do
    endif

    ! IV (j=nj)       
    fc = 4   ! face index
    j = nj
    if( .not. Approx(blocks(b) % face_weights(fc,3),1.0 ) ) then
      do i=1,ni
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,3), i,j,1,i,j,nk)
      end do
    else ! lines in the i direction
      do k=1,nk
        call Distribute_Nodes(b, blocks(b) % face_weights(fc,1), 1,j,k,ni,j,k)
      end do
    endif

    !-------------!
    !   Volumes   !
    !-------------!
    if( .not. Approx( blocks(b) % weights(3), 1.0 ) ) then
      do i=1,ni
        do j=1,nj
          call Distribute_Nodes(b, blocks(b) % weights(3), i,j,1,i,j,nk)
        end do
      end do
    else if( .not. Approx( blocks(b) % weights(1), 1.0 ) ) then
      do k=1,nk
        do j=1,nj
          call Distribute_Nodes(b, blocks(b) % weights(1), 1,j,k,ni,j,k)
        end do
      end do
    else if( .not. Approx( blocks(b) % weights(2), 1.0 ) ) then
      do k=1,nk
        do i=1,ni
          call Distribute_Nodes(b, blocks(b) % weights(2), i,1,k,i,nj,k)
        end do
      end do
    else

      do i=1,ni
        do j=1,nj
          do k=1,nk
            n = nn+(k-1)*ni*nj + (j-1)*ni + i
            call Laplac(b, i, j, k, 0.333, 0.333, 0.334,           &
                                    0.333, 0.333, 0.334,           &
                                    0.333, 0.333, 0.334)
          end do
        end do
      end do
    endif

    !-----------------------------------------!
    !   Set the control volume nodes (CellN)  !
    !    and  the control volume neighbours   !
    !-----------------------------------------!
    ci = ni-1     
    cj = nj-1
    ck = nk-1

    do k=1,ck
      do j=1,cj
        do i=1,ci
          c = nc + (k-1)*ci*cj + (j-1)*ci + i ! cell 
          n = nn + (k-1)*ni*nj + (j-1)*ni + i ! 1st node

          ! Nodes
          CellN(c,1) =n
          CellN(c,2) =n+1
          CellN(c,3) =n+ni
          CellN(c,4) =n+ni+1
          CellN(c,5) =CellN(c,1)+ni*nj
          CellN(c,6) =CellN(c,2)+ni*nj
          CellN(c,7) =CellN(c,3)+ni*nj
          CellN(c,8) =CellN(c,4)+ni*nj

          ! Neighbours
          CellC(c,1) = c-ci*cj
          CellC(c,2) = c-ci 
          CellC(c,3) = c+1
          CellC(c,4) = c+ci
          CellC(c,5) = c-1
          CellC(c,6) = c+ci*cj

          ! This value (-1) is also the default boundary marker
          if(i == 1)  CellC(c,5) =-1
          if(i == ci) CellC(c,3) =-1
          if(j == 1)  CellC(c,2) =-1
          if(j == cj) CellC(c,4) =-1
          if(k == 1)  CellC(c,1) =-1
          if(k == ck) CellC(c,6) =-1

        end do
      end do
    end do

    blocks(b) % n_nodes = nn       ! old number of nodes, for fusion 
    blocks(b) % n_cells = nc       ! old number of volumes, for fusion
    nn = nn + ni*nj*nk
    nc = nc + ci*cj*ck

  end do   ! through blocks 

  !-----------------------------------------!
  !   Insertion of the boundary condition   ! 
  !        and materials information        !
  !-----------------------------------------!

  ! Initialize all the material markers to 1
  do c=1,nc
    material(c) =1
  end do

  do n=1,n_b_cond

    b  = abs(b_cond(n,7))   ! block

    ! Block resolution
    ci=blocks(b) % resolutions(1)-1
    cj=blocks(b) % resolutions(2)-1
    ck=blocks(b) % resolutions(3)-1

    ! Default values
    is=1
    ie=ci
    js=1
    je=cj
    ks=1
    ke=ck

    ! Boundary conditions prescribed with mnemonics
    if(BndFac(n) == 'IMIN') then
      ie=1 
      face = 5
    else if(BndFac(n) == 'IMAX') then 
      is=ci
      face = 3
    else if(BndFac(n) == 'JMIN') then 
      je=1
      face = 2
    else if(BndFac(n) == 'JMAX') then 
      js=cj
      face = 4
    else if(BndFac(n) == 'KMIN') then 
      ke=1
      face = 1
    else if(BndFac(n) == 'KMAX') then 
      ks=ck
      face = 6

    ! Boundary conditions (materials) prescribed explicitly
    !  (error prone and difficult, but might be usefull)
    else   
      is = b_cond(n,1)
      js = b_cond(n,2)
      ks = b_cond(n,3)
      ie = b_cond(n,4)
      je = b_cond(n,5)
      ke = b_cond(n,6)
      face = 0
      if( (is == ie).and.(is ==  1) ) face=5
      if( (is == ie).and.(is == ci) ) face=3
      if( (js == je).and.(js ==  1) ) face=2
      if( (js == je).and.(js == cj) ) face=4
      if( (ks == ke).and.(ks ==  1) ) face=1
      if( (ks == ke).and.(ks == ck) ) face=6
    end if

    do i=is,ie
      do j=js,je
        do k=ks,ke
          c = blocks(b) % n_cells + (k-1)*ci*cj + (j-1)*ci + i   
          if(face /= 0) then 
            CellC(c,face) = -b_cond(n,8) ! marker
          else
            material(c) = b_cond(n,8)    ! material
          end if
        end do
      end do
    end do

  end do  !  n_b_cond

  end subroutine Compute_Node_Coordinates
