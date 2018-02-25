!==============================================================================!
  subroutine Connect_Periodicity(dom, grid)
!------------------------------------------------------------------------------!
!   Solve the cell connectivity for periodic boundary conditions.              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use gen_mod
  use Domain_Mod
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Domain_Type) :: dom
  type(Grid_Type)   :: grid
!----------------------------------[Calling]-----------------------------------!
  logical :: Are_Nodes_Twins
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, n, p                      ! counters
  integer :: b1, b2                          ! block 1 and 2
  integer :: f1, f2                          ! faces of block 1 and 2
  integer :: n11,n12,n13,n14,n21,n22,n23,n24 ! global node numbers
  integer :: p11,p12,p13,p14,p21,p22,p23,p24 ! global node numbers
  integer :: l11,l12,l13,l14,l21,l22,l23,l24 ! local  node numbers
  integer :: i1, j1, i2, j2, k1, k2          ! directions in dom % blocks
  integer :: ig, jg, nig, njg                ! generic plane 
  integer :: ci1, cj1, ck1, ci2, cj2, ck2    ! resolution of dom % blocks
  integer :: c1, c2                          ! cells from block 1, 2
  integer :: ni1, nj1, nk1, ni2, nj2, nk2    ! resolution of dom % blocks
  integer :: n1, n2                          !  from block 1, 2
  integer :: n3, i3, new
  integer :: trans1(3,3), trans2(3,3)
!==============================================================================!

  ! Initialise TwinN.
  do n = 1, grid % max_n_nodes
    TwinN(n,0) = 0
  end do

  !-----------------------------------------------------!
  !   Search through all block and all of their faces   !
  !-----------------------------------------------------!
  do p = 1, n_periodic_cond    
    do b2 = 1, size(dom % blocks)
      do b1 = 1, size(dom % blocks)
        do f2 = 1, 6    ! faces of the second block
          do f1 = 1, 6  ! faces of the first block

           ! Initialize the transformation matrixes             
            do i=1,3
              do j=1,3
                trans1(i,j)=0
                trans2(i,j)=0
              end do
            end do

            n11 = dom % blocks(b1) % faces(f1, 1)
            n12 = dom % blocks(b1) % faces(f1, 2)
            n13 = dom % blocks(b1) % faces(f1, 3)
            n14 = dom % blocks(b1) % faces(f1, 4) 
            n21 = dom % blocks(b2) % faces(f2, 1)
            n22 = dom % blocks(b2) % faces(f2, 2)
            n23 = dom % blocks(b2) % faces(f2, 3)
            n24 = dom % blocks(b2) % faces(f2, 4)

            p11=periodic_cond(p, 1)
            p12=periodic_cond(p, 2)
            p13=periodic_cond(p, 3)
            p14=periodic_cond(p, 4) 
            p21=periodic_cond(p, 5)
            p22=periodic_cond(p, 6)
            p23=periodic_cond(p, 7)
            p24=periodic_cond(p, 8)

          ! Check if they are connected 
          if( ( ((n11 == p11).and.(n13 == p13)) .or.                &
                ((n11 == p14).and.(n13 == p12)) .or.                &
                ((n11 == p13).and.(n13 == p11)) .or.                &
                ((n11 == p12).and.(n13 == p14)) )                   &
                             .and.                                  &
              ( ((n21 == p21).and.(n23 == p23)) .or.                &
                ((n21 == p24).and.(n23 == p22)) .or.                &
                ((n21 == p23).and.(n23 == p21)) .or.                &
                ((n21 == p22).and.(n23 == p24)) ) ) then

              ! Find local nodes (1-8) from blocks 1 and 2 on generic surface
              do n=1,8
                if(dom % blocks(b1) % corners(n) == p11) l11=n
                if(dom % blocks(b1) % corners(n) == p12) l12=n
                if(dom % blocks(b1) % corners(n) == p13) l13=n
                if(dom % blocks(b1) % corners(n) == p14) l14=n
                if(dom % blocks(b2) % corners(n) == p21) l21=n
                if(dom % blocks(b2) % corners(n) == p22) l22=n
                if(dom % blocks(b2) % corners(n) == p23) l23=n
                if(dom % blocks(b2) % corners(n) == p24) l24=n
              end do

               write(*, *) '# Periodicity between:', b1, b2

              ! Direction ig, block 1
              if((l14-l11) == +1) then
                nig = dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,2)=+1
              elseif((l14-l11) == +2) then
                nig = dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,2)=+1
              elseif((l14-l11) == +4) then 
                nig = dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,2)=+1
              elseif((l14-l11) == -1) then 
                nig = dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,1)=nig
                trans1(1,2)=-1
              elseif((l14-l11) == -2) then 
                nig = dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,1)=nig
                trans1(2,2)=-1
              elseif((l14-l11) == -4) then 
                nig = dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,1)=nig
                trans1(3,2)=-1
              endif

              ! Direction jg, block 1 
              if((l12-l11) == +1) then 
                njg = dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,3)=+1
              elseif((l12-l11) == +2) then
                njg = dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,3)=+1
              elseif((l12-l11) == +4) then
                njg = dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,3)=+1
              elseif((l12-l11) == -1) then
                njg = dom % blocks(b1) % resolutions(1)       ! ni from block 1
                trans1(1,1)=njg
                trans1(1,3)=-1
              elseif((l12-l11) == -2) then
                njg = dom % blocks(b1) % resolutions(2)       ! nj from block 1
                trans1(2,1)=njg
                trans1(2,3)=-1
              elseif((l12-l11) == -4) then
                njg = dom % blocks(b1) % resolutions(3)       ! nk from block 1
                trans1(3,1)=njg
                trans1(3,3)=-1
              endif

              ! Direction ig, block 2
              if((l24-l21) == +1) then
                nig = dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,2)=+1
              elseif((l24-l21) == +2) then
                nig = dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,2)=+1
              elseif((l24-l21) == +4) then 
                nig = dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,2)=+1
              elseif((l24-l21) == -1) then 
                nig = dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,1)=nig
                trans2(1,2)=-1
              elseif((l24-l21) == -2) then 
                nig = dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,1)=nig
                trans2(2,2)=-1
              elseif((l24-l21) == -4) then 
                nig = dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,1)=nig
                trans2(3,2)=-1
              endif

              ! Direction jg, block 2 
              if((l22-l21) == +1) then 
                njg = dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,3)=+1
              elseif((l22-l21) == +2) then
                njg = dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,3)=+1
              elseif((l22-l21) == +4) then
                njg = dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,3)=+1
              elseif((l22-l21) == -1) then
                njg = dom % blocks(b2) % resolutions(1)       ! ni from block 2
                trans2(1,1)=njg
                trans2(1,3)=-1
              elseif((l22-l21) == -2) then
                njg = dom % blocks(b2) % resolutions(2)       ! nj from block 2
                trans2(2,1)=njg
                trans2(2,3)=-1
              elseif((l22-l21) == -4) then
                njg = dom % blocks(b2) % resolutions(3)       ! nk from block 2
                trans2(3,1)=njg
                trans2(3,3)=-1
              endif

              ! Set the constant directions
              if(f1 == 1) trans1(3,1)=1
              if(f1 == 2) trans1(2,1)=1
              if(f1 == 3) trans1(1,1)=dom % blocks(b1) % resolutions(1)-1
              if(f1 == 4) trans1(2,1)=dom % blocks(b1) % resolutions(2)-1
              if(f1 == 5) trans1(1,1)=1
              if(f1 == 6) trans1(3,1)=dom % blocks(b1) % resolutions(3)-1

              if(f2 == 1) trans2(3,1)=1
              if(f2 == 2) trans2(2,1)=1
              if(f2 == 3) trans2(1,1)=dom % blocks(b2) % resolutions(1)-1
              if(f2 == 4) trans2(2,1)=dom % blocks(b2) % resolutions(2)-1
              if(f2 == 5) trans2(1,1)=1
              if(f2 == 6) trans2(3,1)=dom % blocks(b2) % resolutions(3)-1

              ! Finally conect the two periodic boundaries
              do jg=1,njg-1              ! through cells only
                do ig=1,nig-1            ! through cells only
                  ci1=dom % blocks(b1) % resolutions(1)-1
                  cj1=dom % blocks(b1) % resolutions(2)-1
                  ck1=dom % blocks(b1) % resolutions(3)-1
                  ci2=dom % blocks(b2) % resolutions(1)-1
                  cj2=dom % blocks(b2) % resolutions(2)-1
                  ck2=dom % blocks(b2) % resolutions(3)-1
                  i1 = trans1(1,1)+trans1(1,2)*ig+trans1(1,3)*jg
                  j1 = trans1(2,1)+trans1(2,2)*ig+trans1(2,3)*jg
                  k1 = trans1(3,1)+trans1(3,2)*ig+trans1(3,3)*jg 
                  i2 = trans2(1,1)+trans2(1,2)*ig+trans2(1,3)*jg
                  j2 = trans2(2,1)+trans2(2,2)*ig+trans2(2,3)*jg
                  k2 = trans2(3,1)+trans2(3,2)*ig+trans2(3,3)*jg
                  c1 = dom % blocks(b1) % n_cells  &
                     + (k1-1)*ci1*cj1 + (j1-1)*ci1 + i1
                  c2 = dom % blocks(b2) % n_cells  &
                     + (k2-1)*ci2*cj2 + (j2-1)*ci2 + i2
                  grid % cells_c(f1,c1) = c2
                  grid % cells_c(f2,c2) = c1
                end do
              end do

              ! Modify the transformation matrices for nodal connection
              if(trans1(1,1)  > 1) trans1(1,1)=trans1(1,1)+1
              if(trans1(2,1)  > 1) trans1(2,1)=trans1(2,1)+1
              if(trans1(3,1)  > 1) trans1(3,1)=trans1(3,1)+1
              if(trans2(1,1)  > 1) trans2(1,1)=trans2(1,1)+1
              if(trans2(2,1)  > 1) trans2(2,1)=trans2(2,1)+1
              if(trans2(3,1)  > 1) trans2(3,1)=trans2(3,1)+1   

              ! Conect the nodes 
              do jg=1,njg                ! through nodes 
                do ig=1,nig              ! through nodes
                  ni1=dom % blocks(b1) % resolutions(1)
                  nj1=dom % blocks(b1) % resolutions(2)
                  nk1=dom % blocks(b1) % resolutions(3)
                  ni2=dom % blocks(b2) % resolutions(1)
                  nj2=dom % blocks(b2) % resolutions(2)
                  nk2=dom % blocks(b2) % resolutions(3)
                  i1 = trans1(1,1) + trans1(1,2)*ig + trans1(1,3)*jg
                  j1 = trans1(2,1) + trans1(2,2)*ig + trans1(2,3)*jg
                  k1 = trans1(3,1) + trans1(3,2)*ig + trans1(3,3)*jg
                  i2 = trans2(1,1) + trans2(1,2)*ig + trans2(1,3)*jg
                  j2 = trans2(2,1) + trans2(2,2)*ig + trans2(2,3)*jg
                  k2 = trans2(3,1) + trans2(3,2)*ig + trans2(3,3)*jg
                  n1 = dom % blocks(b1) % n_nodes  &
                     + (k1-1)*ni1*nj1 + (j1-1)*ni1 + i1
                  n2 = dom % blocks(b2) % n_nodes  &
                     + (k2-1)*ni2*nj2 + (j2-1)*ni2 + i2
                  n1 = NewN(n1)
                  n2 = NewN(n2)

                  ! Check if they are already connected
                  do n=1, TwinN(n1,0)
                    if(Are_Nodes_Twins(n1,n2)) goto 10
                  end do

                  ! If they were not, connect them
                  TwinN(n1,0)=TwinN(n1,0)+1
                  TwinN(n1,TwinN(n1,0))=n2
                  TwinN(n2,0)=TwinN(n2,0)+1
                  TwinN(n2,TwinN(n2,0))=n1

10                end do       ! jg
              end do    ! ig

            end if  ! are they connected ? 

          end do    ! f1
        end do      ! f2
      end do        ! b1
    end do          ! b2 
  end do            ! p periods

  !---------------------!
  !   Twin of my twin   ! 
  !   is also my twin   !
  !---------------------!
  do n1=1,grid % n_nodes
    do i1=1,TwinN(n1,0)
      n2=TwinN(n1,i1) 
      do i2=1,TwinN(n2,0)
        n3=TwinN(n2,i2)   ! blizanci od n2
        new=n3
        do i3=1,TwinN(n1,0)  
          if( (TwinN(n1,i3) == n3) .or. (n3 == n1) ) new=0 
        end do
        if(new == n3) then 
          TwinN(n1,0)=TwinN(n1,0)+1
          TwinN(n1,TwinN(n1,0))=n3
        end if 
      end do
    end do
  end do

  end subroutine
