!==============================================================================!
  subroutine Connect_Blocks
!------------------------------------------------------------------------------!
!   Solve the cell connectivity after block by block grid generation           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, n                         ! counters
  integer :: b1, b2                          ! block 1 and 2
  integer :: f1, f2                          ! faces of block 1 and 2
  integer :: n11,n12,n13,n14,n21,n22,n23,n24 ! global node numbers
  integer :: l11,l12,l13,l14,l21,l22,l23,l24 ! local  node numbers
  integer :: g1, g2, g3, g4                  ! generic points
  integer :: i1, j1, i2, j2, k1, k2          ! directions in blocks
  integer :: ig, jg, nig, njg                ! generic plane 
  integer :: ci1, cj1, ck1, ci2, cj2, ck2    ! resolution of blocks
  integer :: c1, c2                          ! cells from block 1, 2
  integer :: ni1, nj1, nk1, ni2, nj2, nk2    ! resolution of blocks
  integer :: n1, n2                          ! from block 1, 2
  integer :: trans1(3,3), trans2(3,3)
  integer :: del                             ! number of deleted nodes
!==============================================================================!

  ! Initialize the NewN array
  do n=1,NN
    NewN(n)=n
  end do

  if(Nbloc == 1) return 

  ! Initialize the number of deleted nodes
  del=0

  ni1=0; nj1=0; nk1=0; ni2=0; nj2=0; nk2=0 

  !-----------------------------------------------------!
  !   Search through all block and all of their faces   !
  !-----------------------------------------------------!
  do b2=2,Nbloc
    do b1=1,b2-1
      do f2=1,6    ! faces of the second block
        do f1=1,6  ! faces of the first block

          ! Initialize the transformation matrixes             
          do i=1,3
            do j=1,3
              trans1(i,j)=0
              trans2(i,j)=0
            end do
          end do

          n11=block_faces(b1, f1, 1)
          n12=block_faces(b1, f1, 2)
          n13=block_faces(b1, f1, 3)
          n14=block_faces(b1, f1, 4) 
          n21=block_faces(b2, f2, 1)
          n22=block_faces(b2, f2, 2)
          n23=block_faces(b2, f2, 3)
          n24=block_faces(b2, f2, 4)

          ! Check if they are connected 
          if( ((n11 == n21).and.(n13 == n23)) .or.  &
              ((n11 == n24).and.(n13 == n22)) .or.  &
              ((n11 == n23).and.(n13 == n21)) .or.  &
              ((n11 == n22).and.(n13 == n24)) ) then

            ! Define generic surface (g1-g4 are in essence not needed)
            g1=n11
            g2=n12
            g3=n13
            g4=n14

            ! Find local nodes (1-8) from blocks 1 and 2 on generic surface
            do n=1,8
              if(block_points(b1,n) == g1) l11=n
              if(block_points(b2,n) == g1) l21=n
              if(block_points(b1,n) == g2) l12=n
              if(block_points(b2,n) == g2) l22=n
              if(block_points(b1,n) == g3) l13=n
              if(block_points(b2,n) == g3) l23=n
              if(block_points(b1,n) == g4) l14=n
              if(block_points(b2,n) == g4) l24=n
            end do

            ! Direction ig, block 1
            if((l14-l11) == +1) then
              nig = block_resolutions(b1,1)       ! ni from block 1
              trans1(1,2)=+1
            elseif((l14-l11) == +2) then
              nig = block_resolutions(b1,2)       ! nj from block 1
              trans1(2,2)=+1
            elseif((l14-l11) == +4) then 
              nig = block_resolutions(b1,3)       ! nk from block 1
              trans1(3,2)=+1
            elseif((l14-l11) == -1) then 
              nig = block_resolutions(b1,1)       ! ni from block 1
              trans1(1,1)=nig
              trans1(1,2)=-1
            elseif((l14-l11) == -2) then 
              nig = block_resolutions(b1,2)       ! nj from block 1
              trans1(2,1)=nig
              trans1(2,2)=-1
            elseif((l14-l11) == -4) then 
              nig = block_resolutions(b1,3)       ! nk from block 1
              trans1(3,1)=nig
              trans1(3,2)=-1
            endif

            ! Direction jg, block 1 
            if((l12-l11) == +1) then 
              njg = block_resolutions(b1,1)       ! ni from block 1
              trans1(1,3)=+1
            elseif((l12-l11) == +2) then
              njg = block_resolutions(b1,2)       ! nj from block 1
              trans1(2,3)=+1
            elseif((l12-l11) == +4) then
              njg = block_resolutions(b1,3)       ! nk from block 1
              trans1(3,3)=+1
            elseif((l12-l11) == -1) then
              njg = block_resolutions(b1,1)       ! ni from block 1
              trans1(1,1)=njg
              trans1(1,3)=-1
            elseif((l12-l11) == -2) then
              njg = block_resolutions(b1,2)       ! nj from block 1
              trans1(2,1)=njg
              trans1(2,3)=-1
            elseif((l12-l11) == -4) then
              njg = block_resolutions(b1,3)       ! nk from block 1
              trans1(3,1)=njg
              trans1(3,3)=-1
            endif

            ! Direction ig, block 2
            if((l24-l21) == +1) then
              nig = block_resolutions(b2,1)       ! ni from block 2
              trans2(1,2)=+1
            elseif((l24-l21) == +2) then
              nig = block_resolutions(b2,2)       ! nj from block 2
              trans2(2,2)=+1
            elseif((l24-l21) == +4) then 
              nig = block_resolutions(b2,3)       ! nk from block 2
              trans2(3,2)=+1
            elseif((l24-l21) == -1) then 
              nig = block_resolutions(b2,1)       ! ni from block 2
              trans2(1,1)=nig
              trans2(1,2)=-1
            elseif((l24-l21) == -2) then 
              nig = block_resolutions(b2,2)       ! nj from block 2
              trans2(2,1)=nig
              trans2(2,2)=-1
            elseif((l24-l21) == -4) then 
              nig = block_resolutions(b2,3)       ! nk from block 2
              trans2(3,1)=nig
              trans2(3,2)=-1
            endif

            ! Direction jg, block 2 
            if((l22-l21) == +1) then 
              njg = block_resolutions(b2,1)       ! ni from block 2
              trans2(1,3)=+1
            elseif((l22-l21) == +2) then
              njg = block_resolutions(b2,2)       ! nj from block 2
              trans2(2,3)=+1
            elseif((l22-l21) == +4) then
              njg = block_resolutions(b2,3)       ! nk from block 2
              trans2(3,3)=+1
            elseif((l22-l21) == -1) then
              njg = block_resolutions(b2,1)       ! ni from block 2
              trans2(1,1)=njg
              trans2(1,3)=-1
            elseif((l22-l21) == -2) then
              njg = block_resolutions(b2,2)       ! nj from block 2
              trans2(2,1)=njg
              trans2(2,3)=-1
            elseif((l22-l21) == -4) then
              njg = block_resolutions(b2,3)       ! nk from block 2
              trans2(3,1)=njg
              trans2(3,3)=-1
            endif

            ! Set the constant directions
            if(f1 == 1) trans1(3,1)=1
            if(f1 == 2) trans1(2,1)=1
            if(f1 == 3) trans1(1,1)=block_resolutions(b1,1)-1
            if(f1 == 4) trans1(2,1)=block_resolutions(b1,2)-1
            if(f1 == 5) trans1(1,1)=1
            if(f1 == 6) trans1(3,1)=block_resolutions(b1,3)-1

            if(f2 == 1) trans2(3,1)=1
            if(f2 == 2) trans2(2,1)=1
            if(f2 == 3) trans2(1,1)=block_resolutions(b2,1)-1
            if(f2 == 4) trans2(2,1)=block_resolutions(b2,2)-1
            if(f2 == 5) trans2(1,1)=1
            if(f2 == 6) trans2(3,1)=block_resolutions(b2,3)-1

            ! Finally conect the two blocks
            do jg=1,njg-1              ! through volumes only
              do ig=1,nig-1            ! through volumes only
                ci1=block_resolutions(b1,1)-1
                cj1=block_resolutions(b1,2)-1
                ck1=block_resolutions(b1,3)-1
                ci2=block_resolutions(b2,1)-1
                cj2=block_resolutions(b2,2)-1
                ck2=block_resolutions(b2,3)-1
                i1 = trans1(1,1) + trans1(1,2)*ig + trans1(1,3)*jg
                j1 = trans1(2,1) + trans1(2,2)*ig + trans1(2,3)*jg
                k1 = trans1(3,1) + trans1(3,2)*ig + trans1(3,3)*jg 
                i2 = trans2(1,1) + trans2(1,2)*ig + trans2(1,3)*jg
                j2 = trans2(2,1) + trans2(2,2)*ig + trans2(2,3)*jg
                k2 = trans2(3,1) + trans2(3,2)*ig + trans2(3,3)*jg
                c1 = block_resolutions(b1,6)  &
                     + (k1-1)*ci1*cj1 + (j1-1)*ci1 + i1
                c2 = block_resolutions(b2,6)  &
                     + (k2-1)*ci2*cj2 + (j2-1)*ci2 + i2
                CellC(c1, f1) = c2
                CellC(c2, f2) = c1
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
                ni1=block_resolutions(b1,1)
                nj1=block_resolutions(b1,2)
                nk1=block_resolutions(b1,3)
                ni2=block_resolutions(b2,1)
                nj2=block_resolutions(b2,2)
                nk2=block_resolutions(b2,3)
                i1 = trans1(1,1) + trans1(1,2)*ig + trans1(1,3)*jg
                j1 = trans1(2,1) + trans1(2,2)*ig + trans1(2,3)*jg
                k1 = trans1(3,1) + trans1(3,2)*ig + trans1(3,3)*jg
                i2 = trans2(1,1) + trans2(1,2)*ig + trans2(1,3)*jg
                j2 = trans2(2,1) + trans2(2,2)*ig + trans2(2,3)*jg
                k2 = trans2(3,1) + trans2(3,2)*ig + trans2(3,3)*jg
                n1 = block_resolutions(b1,5)  &
                     + (k1-1)*ni1*nj1 + (j1-1)*ni1 + i1
                n2 = block_resolutions(b2,5)  &
                     + (k2-1)*ni2*nj2 + (j2-1)*ni2 + i2
                NewN(n2) = NewN(n1)
              end do
            end do

          end if  ! are they connected ? 

        end do    ! f1
      end do      ! f2
    end do        ! b1

    ! Update node numbers
    do n=block_resolutions(b2,5)+1, block_resolutions(b2,5)+ni2*nj2*nk2
      if(NewN(n) /= n) del=del+1
      if(NewN(n) == n) NewN(n)=NewN(n)-del
    end do 

  end do          ! b2 


  do n=1,NN
    x_node(NewN(n))=x_node(n)
    y_node(NewN(n))=y_node(n)
    z_node(NewN(n))=z_node(n)
  end do

  NN=NN-del

  ! Skip the merged points in the node() structure
  do i=1,NC
    do n=1,8
      CellN(i,n)=NewN(CellN(i,n))
    end do
  end do

  write(6, '(I8)') del       

  end subroutine Connect_Blocks   
