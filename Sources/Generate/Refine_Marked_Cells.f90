!==============================================================================!
  subroutine Refine_Marked_Cells(lev)
!------------------------------------------------------------------------------!
!   Refine the marked cells.                                                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Parameters]---------------------------------!
  integer :: lev 
!----------------------------------[Calling]-----------------------------------!
  logical :: Are_Nodes_Twins
  integer :: Which_Node
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n_cells_old, c1, c2, c3, c4, c5, c6
  integer :: cr1, cr2, cr3, cr4, cr5, cr6, cr7, cr8
  integer :: n, n_nodes_old, n1, n2, n3, n4, n5, n6, n7, n8
  integer :: n12,n13,n24,n34,n15,n26,n37,n48,n56,n57,n68,n78
  integer :: nF1, nF2, nF3, nF4, nF5, nF6, n0
  integer :: del   ! number of deleted nodes 
  integer :: nA, nA0, nA1, nA2, nB, nB0, nB1, nB2
  integer :: NN2, NN4, NN8
!==============================================================================!
!                                                                              !
!                               c6      c3                                     !
!                               |      /                                       !
!                         8-----|---------6                                    !
!                        /  cr8 |/  cr6  /|                                    !
!                       /-------+-------/ |                                    !
!                      /       /       /| |                                    !
!                     7-------+-------5 | |                                    !
!                     |       |       | |/|                                    !
!                c4---|  cr7  |  cr5  | +-----c2                               !
!                     |       |       |/| |                                    !
!                     +-------+-------+ | 2                                    !
!                     |      /|       | |/                                     !
!                     |  cr3/ |  cr1  | /                                      !
!                     |    /  |       |/                                       !
!                     3---c5--+-------1                                        !
!                               |                                              !
!                               c1                                             !
!                                                                              !
!------------------------------------------------------------------------------!

  write(*,*) 'Refine: Number of nodes: ', NN 
  write(*,*) '        Number of cells: ', NC 

  n_cells_old = NC 
  n_nodes_old = NN
  NN2 = 0 
  NN4 = 0
  NN8 = 0

  !---------------------!
  !                     !
  !   Count new celss   !
  !                     !
  !---------------------!
  do c=1,n_cells_old
    if(CelMar(c)  ==  -1) then 
      NC = NC + 8
      CelMar(c) = NC   ! now points to cr8
    end if
  end do

  do c=1,n_cells_old

    c1=CellC(c,1)
    c2=CellC(c,2)
    c3=CellC(c,3)
    c4=CellC(c,4)
    c5=CellC(c,5)
    c6=CellC(c,6)

    n1=CellN(c,1)
    n2=CellN(c,2)
    n3=CellN(c,3)
    n4=CellN(c,4)
    n5=CellN(c,5)
    n6=CellN(c,6)
    n7=CellN(c,7)
    n8=CellN(c,8)

    !-------------------!
    !                   !
    !   Refined cells   !
    !                   !
    !-------------------!
    if(CelMar(c)   >  0) then  ! only refined

      !------------------------------------!
      !   Take care of neighboring cells   !
      !------------------------------------!
      cr1=CelMar(c)-7
      cr2=CelMar(c)-6
      cr3=CelMar(c)-5
      cr4=CelMar(c)-4
      cr5=CelMar(c)-3
      cr6=CelMar(c)-2
      cr7=CelMar(c)-1
      cr8=CelMar(c)

      material(cr1) = material(c) 
      material(cr2) = material(c) 
      material(cr3) = material(c) 
      material(cr4) = material(c) 
      material(cr5) = material(c) 
      material(cr6) = material(c) 
      material(cr7) = material(c) 
      material(cr8) = material(c) 

      !-----------------------------------------------!
      !   Internal links do not depend on neighbors   !
      !-----------------------------------------------!

      ! 6
      CellC(cr1,6) = cr5
      CellC(cr2,6) = cr6
      CellC(cr3,6) = cr7
      CellC(cr4,6) = cr8

      ! 5
      CellC(cr2,5) = cr1
      CellC(cr4,5) = cr3
      CellC(cr6,5) = cr5
      CellC(cr8,5) = cr7

      ! 4
      CellC(cr1,4) = cr3
      CellC(cr2,4) = cr4
      CellC(cr5,4) = cr7
      CellC(cr6,4) = cr8

      ! 3
      CellC(cr1,3) = cr2
      CellC(cr3,3) = cr4
      CellC(cr5,3) = cr6
      CellC(cr7,3) = cr8

      ! 2
      CellC(cr3,2) = cr1
      CellC(cr4,2) = cr2
      CellC(cr7,2) = cr5
      CellC(cr8,2) = cr6

      ! 1
      CellC(cr5,1) = cr1
      CellC(cr6,1) = cr2
      CellC(cr7,1) = cr3
      CellC(cr8,1) = cr4

      !-------------------------!
      !   Level of refinement   !
      !-------------------------!
      level(cr1) = lev
      level(cr2) = lev
      level(cr3) = lev
      level(cr4) = lev
      level(cr5) = lev
      level(cr6) = lev
      level(cr7) = lev
      level(cr8) = lev

      !----------------------------------------!         
      !   External links depend on neighbors   !
      !----------------------------------------!
      if(CelMar(c1)  ==  0) then  ! neighbor 1 not refined
        CellC(cr1,1) = c1
        CellC(cr2,1) = c1
        CellC(cr3,1) = c1
        CellC(cr4,1) = c1
      else                        ! neighbor 1 refined
        CellC(cr1,1) = CelMar(c1) - 8 + Which_Node(c1,n1)
        CellC(cr2,1) = CelMar(c1) - 8 + Which_Node(c1,n2)
        CellC(cr3,1) = CelMar(c1) - 8 + Which_Node(c1,n3)
        CellC(cr4,1) = CelMar(c1) - 8 + Which_Node(c1,n4)
      endif           

      if(CelMar(c2)  ==  0) then  ! neighbor 2 not refined
        CellC(cr1,2) = c2
        CellC(cr2,2) = c2
        CellC(cr5,2) = c2
        CellC(cr6,2) = c2
      else                        ! neighbor 2 refined
        CellC(cr1,2) = CelMar(c2) - 8 + Which_Node(c2, n1)
        CellC(cr2,2) = CelMar(c2) - 8 + Which_Node(c2, n2)
        CellC(cr5,2) = CelMar(c2) - 8 + Which_Node(c2, n5)
        CellC(cr6,2) = CelMar(c2) - 8 + Which_Node(c2, n6)
      endif           

      if(CelMar(c3)  ==  0) then  ! neighbor 3 not refined
        CellC(cr2,3) = c3
        CellC(cr4,3) = c3
        CellC(cr6,3) = c3
        CellC(cr8,3) = c3
      else                        ! neighbor 3 refined
        CellC(cr2,3) = CelMar(c3) - 8 + Which_Node(c3, n2)
        CellC(cr4,3) = CelMar(c3) - 8 + Which_Node(c3, n4)
        CellC(cr6,3) = CelMar(c3) - 8 + Which_Node(c3, n6)
        CellC(cr8,3) = CelMar(c3) - 8 + Which_Node(c3, n8)
      endif           

      if(CelMar(c4)  ==  0) then  ! neighbor 4 not refine
        CellC(cr3,4) = c4
        CellC(cr4,4) = c4
        CellC(cr7,4) = c4
        CellC(cr8,4) = c4
      else                        ! neighbor 4 refine
        CellC(cr3,4) = CelMar(c4) - 8 + Which_Node(c4, n3)
        CellC(cr4,4) = CelMar(c4) - 8 + Which_Node(c4, n4)
        CellC(cr7,4) = CelMar(c4) - 8 + Which_Node(c4, n7)
        CellC(cr8,4) = CelMar(c4) - 8 + Which_Node(c4, n8)
      endif           

      if(CelMar(c5)  ==  0) then  ! neighbor 5 not refined
        CellC(cr1,5) = c5
        CellC(cr3,5) = c5
        CellC(cr5,5) = c5
        CellC(cr7,5) = c5
      else                        ! neighbor 5 refined
        CellC(cr1,5) = CelMar(c5) - 8 + Which_Node(c5, n1)
        CellC(cr3,5) = CelMar(c5) - 8 + Which_Node(c5, n3)
        CellC(cr5,5) = CelMar(c5) - 8 + Which_Node(c5, n5)
        CellC(cr7,5) = CelMar(c5) - 8 + Which_Node(c5, n7)
      endif

      if(CelMar(c6)  ==  0) then  ! neighbor 6 not refined
        CellC(cr5,6) = c6
        CellC(cr6,6) = c6
        CellC(cr7,6) = c6
        CellC(cr8,6) = c6
      else                        ! neighbor 6 refined
        CellC(cr5,6) = CelMar(c6) - 8 + Which_Node(c6, n5)
        CellC(cr6,6) = CelMar(c6) - 8 + Which_Node(c6, n6)
        CellC(cr7,6) = CelMar(c6) - 8 + Which_Node(c6, n7)
        CellC(cr8,6) = CelMar(c6) - 8 + Which_Node(c6, n8)
      endif           

    !------------------------!
    !   Take care of nodes   !
    !------------------------!

    !------------------------!
    !   Nodes on the edges   !
    !------------------------!
    n12 = 0     !
    n13 = 0     !         8-----n68-----6
    n24 = 0     !        /|            /|
    n34 = 0     !      n78|          n56|
    n15 = 0     !      / n48         / n26
    n26 = 0     !     7-----n57-----5   |
    n37 = 0     !     |   |         |   |
    n48 = 0     !     |   4- - -n24-| - 2
    n56 = 0     !    n37 /         n15 / 
    n57 = 0     !     |n34          |n12
    n68 = 0     !     |/            |/
    n78 = 0     !     3-----n13-----1

    ! n12
    do n=1,NN2
      if( ( NodeN2(n,1) == n1 .and. NodeN2(n,2) == n2 ) .or.  &
          ( NodeN2(n,1) == n2 .and. NodeN2(n,2) == n1) ) then
        n12 = NodeN2(n,0) 
      end if
    end do
    if (n12 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n12 = NN
      NodeN2(NN2,0) = n12
      NodeN2(NN2,1) = n1
      NodeN2(NN2,2) = n2
      x_node(n12) = 0.5 * (x_node(n1)+x_node(n2))
      y_node(n12) = 0.5 * (y_node(n1)+y_node(n2))
      z_node(n12) = 0.5 * (z_node(n1)+z_node(n2))
    end if 

    ! n13
    do n=1,NN2
      if( ( NodeN2(n,1) == n1 .and. NodeN2(n,2) == n3 ) .or.  &
          ( NodeN2(n,1) == n3 .and. NodeN2(n,2) == n1) ) then
        n13 = NodeN2(n,0) 
      end if
    end do
    if (n13 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n13 = NN
      NodeN2(NN2,0) = n13
      NodeN2(NN2,1) = n1
      NodeN2(NN2,2) = n3
      x_node(n13) = 0.5 * (x_node(n1)+x_node(n3))
      y_node(n13) = 0.5 * (y_node(n1)+y_node(n3))
      z_node(n13) = 0.5 * (z_node(n1)+z_node(n3))
    end if 

    ! n24
    do n=1,NN2
      if( ( NodeN2(n,1) == n2 .and. NodeN2(n,2) == n4 ) .or.  &
          ( NodeN2(n,1) == n4 .and. NodeN2(n,2) == n2) ) then
        n24 = NodeN2(n,0) 
      end if
    end do
    if (n24 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n24 = NN
      NodeN2(NN2,0) = n24
      NodeN2(NN2,1) = n2
      NodeN2(NN2,2) = n4
      x_node(n24) = 0.5 * (x_node(n2)+x_node(n4))
      y_node(n24) = 0.5 * (y_node(n2)+y_node(n4))
      z_node(n24) = 0.5 * (z_node(n2)+z_node(n4))
    end if 

    ! n34
    do n=1,NN2
      if( ( NodeN2(n,1) == n3 .and. NodeN2(n,2) == n4 ) .or.  &
          ( NodeN2(n,1) == n4 .and. NodeN2(n,2) == n3) ) then
        n34 = NodeN2(n,0) 
      end if
    end do
    if (n34 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n34 = NN
      NodeN2(NN2,0) = n34
      NodeN2(NN2,1) = n3
      NodeN2(NN2,2) = n4
      x_node(n34) = 0.5 * (x_node(n3)+x_node(n4))
      y_node(n34) = 0.5 * (y_node(n3)+y_node(n4))
      z_node(n34) = 0.5 * (z_node(n3)+z_node(n4))
    end if 

    ! n15
    do n=1,NN2
      if( ( NodeN2(n,1) == n1 .and. NodeN2(n,2) == n5 ) .or.  &
          ( NodeN2(n,1) == n5 .and. NodeN2(n,2) == n1) ) then
        n15 = NodeN2(n,0) 
      end if
    end do
    if (n15 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n15 = NN
      NodeN2(NN2,0) = n15
      NodeN2(NN2,1) = n1
      NodeN2(NN2,2) = n5
      x_node(n15) = 0.5 * (x_node(n1)+x_node(n5))
      y_node(n15) = 0.5 * (y_node(n1)+y_node(n5))
      z_node(n15) = 0.5 * (z_node(n1)+z_node(n5))
    end if 

    ! n26
    do n=1,NN2
      if( ( NodeN2(n,1) == n2 .and. NodeN2(n,2) == n6 ) .or.  &
          ( NodeN2(n,1) == n6 .and. NodeN2(n,2) == n2) ) then
        n26 = NodeN2(n,0) 
      end if
    end do
    if (n26 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n26 = NN
      NodeN2(NN2,0) = n26
      NodeN2(NN2,1) = n2
      NodeN2(NN2,2) = n6
      x_node(n26) = 0.5 * (x_node(n2)+x_node(n6))
      y_node(n26) = 0.5 * (y_node(n2)+y_node(n6))
      z_node(n26) = 0.5 * (z_node(n2)+z_node(n6))
    end if 

    ! n37
    do n=1,NN2
      if( ( NodeN2(n,1) == n3 .and. NodeN2(n,2) == n7 ) .or.  &
          ( NodeN2(n,1) == n7 .and. NodeN2(n,2) == n3) ) then
        n37 = NodeN2(n,0) 
      end if
    end do
    if (n37 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n37 = NN
      NodeN2(NN2,0) = n37
      NodeN2(NN2,1) = n3
      NodeN2(NN2,2) = n7
      x_node(n37) = 0.5 * (x_node(n3)+x_node(n7))
      y_node(n37) = 0.5 * (y_node(n3)+y_node(n7))
      z_node(n37) = 0.5 * (z_node(n3)+z_node(n7))
    end if 

    ! n48
    do n=1,NN2
      if( ( NodeN2(n,1) == n4 .and. NodeN2(n,2) == n8 ) .or.  &
          ( NodeN2(n,1) == n8 .and. NodeN2(n,2) == n4) ) then
        n48 = NodeN2(n,0) 
      end if
    end do
    if (n48 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n48 = NN
      NodeN2(NN2,0) = n48
      NodeN2(NN2,1) = n4
      NodeN2(NN2,2) = n8
      x_node(n48) = 0.5 * (x_node(n4)+x_node(n8))
      y_node(n48) = 0.5 * (y_node(n4)+y_node(n8))
      z_node(n48) = 0.5 * (z_node(n4)+z_node(n8))
    end if 

    ! n56
    do n=1,NN2
      if( ( NodeN2(n,1) == n5 .and. NodeN2(n,2) == n6 ) .or.  &
          ( NodeN2(n,1) == n6 .and. NodeN2(n,2) == n5) ) then
        n56 = NodeN2(n,0) 
      end if
    end do
    if (n56 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n56 = NN
      NodeN2(NN2,0) = n56
      NodeN2(NN2,1) = n5
      NodeN2(NN2,2) = n6
      x_node(n56) = 0.5 * (x_node(n5)+x_node(n6))
      y_node(n56) = 0.5 * (y_node(n5)+y_node(n6))
      z_node(n56) = 0.5 * (z_node(n5)+z_node(n6))
    end if 

    ! n57
    do n=1,NN2
      if( ( NodeN2(n,1) == n5 .and. NodeN2(n,2) == n7 ) .or.  &
          ( NodeN2(n,1) == n7 .and. NodeN2(n,2) == n5) ) then
        n57 = NodeN2(n,0) 
      end if
    end do
    if (n57 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n57 = NN
      NodeN2(NN2,0) = n57
      NodeN2(NN2,1) = n5
      NodeN2(NN2,2) = n7
      x_node(n57) = 0.5 * (x_node(n5)+x_node(n7))
      y_node(n57) = 0.5 * (y_node(n5)+y_node(n7))
      z_node(n57) = 0.5 * (z_node(n5)+z_node(n7))
    end if 

    ! n68 
    do n=1,NN2
      if( ( NodeN2(n,1) == n6 .and. NodeN2(n,2) == n8 ) .or.  &
          ( NodeN2(n,1) == n8 .and. NodeN2(n,2) == n6) ) then
        n68 = NodeN2(n,0) 
      end if
    end do
    if (n68 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n68 = NN
      NodeN2(NN2,0) = n68
      NodeN2(NN2,1) = n6
      NodeN2(NN2,2) = n8
      x_node(n68) = 0.5 * (x_node(n6)+x_node(n8))
      y_node(n68) = 0.5 * (y_node(n6)+y_node(n8))
      z_node(n68) = 0.5 * (z_node(n6)+z_node(n8))
    end if 

    ! n78
    do n=1,NN2
      if( ( NodeN2(n,1) == n7 .and. NodeN2(n,2) == n8 ) .or.  &
          ( NodeN2(n,1) == n8 .and. NodeN2(n,2) == n7) ) then
        n78 = NodeN2(n,0) 
      end if
    end do
    if (n78 == 0) then
      NN2 = NN2 + 1
      NN  = NN  + 1
      n78 = NN
      NodeN2(NN2,0) = n78
      NodeN2(NN2,1) = n7
      NodeN2(NN2,2) = n8
      x_node(n78) = 0.5 * (x_node(n7)+x_node(n8))
      y_node(n78) = 0.5 * (y_node(n7)+y_node(n8))
      z_node(n78) = 0.5 * (z_node(n7)+z_node(n8))
    end if 

    !-------------------------!
    !   Then nodes on faces   !
    !-------------------------!
    nF1 = 0
    nF2 = 0
    nF3 = 0
    nF4 = 0
    nF5 = 0
    nF6 = 0

    ! nF1
    do n=1,NN4
      if( ( NodeN4(n,1) == n1 .and. NodeN4(n,4) == n4 ) .or.  &
          ( NodeN4(n,1) == n4 .and. NodeN4(n,4) == n1 ) .or.  &
          ( NodeN4(n,1) == n2 .and. NodeN4(n,4) == n3 ) .or.  &
          ( NodeN4(n,1) == n3 .and. NodeN4(n,4) == n2 ) ) then
        nF1 = NodeN4(n,0) 
      end if
    end do
    if (nF1 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF1 = NN
      NodeN4(NN4,0) = nF1
      NodeN4(NN4,1) = n1
      NodeN4(NN4,2) = n2
      NodeN4(NN4,3) = n3
      NodeN4(NN4,4) = n4
      x_node(nF1) = 0.25 * (x_node(n1)+x_node(n2)+x_node(n3)+x_node(n4))
      y_node(nF1) = 0.25 * (y_node(n1)+y_node(n2)+y_node(n3)+y_node(n4))
      z_node(nF1) = 0.25 * (z_node(n1)+z_node(n2)+z_node(n3)+z_node(n4))
    end if 

    ! nF2
    do n=1,NN4
      if( ( NodeN4(n,1) == n1 .and. NodeN4(n,4) == n6 ) .or.  &
          ( NodeN4(n,1) == n6 .and. NodeN4(n,4) == n1 ) .or.  &
          ( NodeN4(n,1) == n2 .and. NodeN4(n,4) == n5 ) .or.  &
          ( NodeN4(n,1) == n5 .and. NodeN4(n,4) == n2 ) ) then
        nF2 = NodeN4(n,0) 
      end if
    end do
    if (nF2 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF2 = NN
      NodeN4(NN4,0) = nF2
      NodeN4(NN4,1) = n1
      NodeN4(NN4,2) = n2
      NodeN4(NN4,3) = n5
      NodeN4(NN4,4) = n6
      x_node(nF2) = 0.25 * (x_node(n1)+x_node(n2)+x_node(n5)+x_node(n6))
      y_node(nF2) = 0.25 * (y_node(n1)+y_node(n2)+y_node(n5)+y_node(n6))
      z_node(nF2) = 0.25 * (z_node(n1)+z_node(n2)+z_node(n5)+z_node(n6))
    end if 

    ! nF3
    do n=1,NN4
      if( ( NodeN4(n,1) == n2 .and. NodeN4(n,4) == n8 ) .or.  &
          ( NodeN4(n,1) == n8 .and. NodeN4(n,4) == n2 ) .or.  &
          ( NodeN4(n,1) == n4 .and. NodeN4(n,4) == n6 ) .or.  &
          ( NodeN4(n,1) == n6 .and. NodeN4(n,4) == n4 ) ) then
        nF3 = NodeN4(n,0) 
      end if
    end do
    if (nF3 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF3 = NN
      NodeN4(NN4,0) = nF3
      NodeN4(NN4,1) = n2
      NodeN4(NN4,2) = n4
      NodeN4(NN4,3) = n6
      NodeN4(NN4,4) = n8
      x_node(nF3) = 0.25 * (x_node(n2)+x_node(n4)+x_node(n6)+x_node(n8))
      y_node(nF3) = 0.25 * (y_node(n2)+y_node(n4)+y_node(n6)+y_node(n8))
      z_node(nF3) = 0.25 * (z_node(n2)+z_node(n4)+z_node(n6)+z_node(n8))
    end if 

    ! nF4
    do n=1,NN4
      if( ( NodeN4(n,1) == n3 .and. NodeN4(n,4) == n8 ) .or.  &
          ( NodeN4(n,1) == n8 .and. NodeN4(n,4) == n3 ) .or.  &
          ( NodeN4(n,1) == n4 .and. NodeN4(n,4) == n7 ) .or.  &
          ( NodeN4(n,1) == n7 .and. NodeN4(n,4) == n4 ) ) then
        nF4 = NodeN4(n,0) 
      end if
    end do
    if (nF4 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF4 = NN
      NodeN4(NN4,0) = nF4
      NodeN4(NN4,1) = n3
      NodeN4(NN4,2) = n4
      NodeN4(NN4,3) = n7
      NodeN4(NN4,4) = n8
      x_node(nF4) = 0.25 * (x_node(n3)+x_node(n4)+x_node(n7)+x_node(n8))
      y_node(nF4) = 0.25 * (y_node(n3)+y_node(n4)+y_node(n7)+y_node(n8))
      z_node(nF4) = 0.25 * (z_node(n3)+z_node(n4)+z_node(n7)+z_node(n8))
    end if 

    ! nF5
    do n=1,NN4
      if( ( NodeN4(n,1) == n1 .and. NodeN4(n,4) == n7 ) .or.  &
          ( NodeN4(n,1) == n7 .and. NodeN4(n,4) == n1 ) .or.  &
          ( NodeN4(n,1) == n3 .and. NodeN4(n,4) == n5 ) .or.  &
          ( NodeN4(n,1) == n5 .and. NodeN4(n,4) == n3 ) ) then
        nF5 = NodeN4(n,0) 
      end if
    end do
    if (nF5 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF5 = NN
      NodeN4(NN4,0) = nF5
      NodeN4(NN4,1) = n1
      NodeN4(NN4,2) = n3
      NodeN4(NN4,3) = n5
      NodeN4(NN4,4) = n7
      x_node(nF5) = 0.25 * (x_node(n1)+x_node(n3)+x_node(n5)+x_node(n7))
      y_node(nF5) = 0.25 * (y_node(n1)+y_node(n3)+y_node(n5)+y_node(n7))
      z_node(nF5) = 0.25 * (z_node(n1)+z_node(n3)+z_node(n5)+z_node(n7))
    end if 

    ! nF6
    do n=1,NN4
      if( ( NodeN4(n,1) == n5 .and. NodeN4(n,4) == n8 ) .or.  &
          ( NodeN4(n,1) == n8 .and. NodeN4(n,4) == n5 ) .or.  &
          ( NodeN4(n,1) == n6 .and. NodeN4(n,4) == n7 ) .or.  &
          ( NodeN4(n,1) == n7 .and. NodeN4(n,4) == n6 ) ) then
        nF6 = NodeN4(n,0) 
      end if
    end do
    if (nF6 == 0) then
      NN4 = NN4 + 1
      NN  = NN  + 1
      nF6 = NN
      NodeN4(NN4,0) = nF6
      NodeN4(NN4,1) = n5
      NodeN4(NN4,2) = n6
      NodeN4(NN4,3) = n7
      NodeN4(NN4,4) = n8
      x_node(nF6) = 0.25 * (x_node(n5)+x_node(n6)+x_node(n7)+x_node(n8))
      y_node(nF6) = 0.25 * (y_node(n5)+y_node(n6)+y_node(n7)+y_node(n8))
      z_node(nF6) = 0.25 * (z_node(n5)+z_node(n6)+z_node(n7)+z_node(n8))
    end if 

    !----------------------------------------!
    !   Eventually, the node in the middle   !
    !----------------------------------------!
    NN8 = NN8 + 1
    NN  = NN + 1
    n0  = NN
    NodeN8(NN8,0) = n0 
    NodeN8(NN8,1) = n1
    NodeN8(NN8,2) = n2
    NodeN8(NN8,3) = n3
    NodeN8(NN8,4) = n4
    NodeN8(NN8,5) = n5
    NodeN8(NN8,6) = n6
    NodeN8(NN8,7) = n7
    NodeN8(NN8,8) = n8
    x_node(n0) = 0.125*(x_node(n1)+x_node(n2)+x_node(n3)+x_node(n4)+  &
                        x_node(n5)+x_node(n6)+x_node(n7)+x_node(n8))
    y_node(n0) = 0.125*(y_node(n1)+y_node(n2)+y_node(n3)+y_node(n4)+  &
                        y_node(n5)+y_node(n6)+y_node(n7)+y_node(n8))
    z_node(n0) = 0.125*(z_node(n1)+z_node(n2)+z_node(n3)+z_node(n4)+  &
                        z_node(n5)+z_node(n6)+z_node(n7)+z_node(n8))

    !----------------------------!
    !   Set nodes to new cells   !
    !----------------------------!

    ! cr1 -!
    CellN(cr1,1) = n1
    CellN(cr1,2) = n12
    CellN(cr1,3) = n13
    CellN(cr1,4) = nF1
    CellN(cr1,5) = n15
    CellN(cr1,6) = nF2
    CellN(cr1,7) = nF5
    CellN(cr1,8) = n0 

    ! cr2 -!
    CellN(cr2,1) = n12
    CellN(cr2,2) = n2 
    CellN(cr2,3) = nF1
    CellN(cr2,4) = n24
    CellN(cr2,5) = nF2
    CellN(cr2,6) = n26 
    CellN(cr2,7) = n0 
    CellN(cr2,8) = nF3

    ! cr3 -!
    CellN(cr3,1) = n13
    CellN(cr3,2) = nF1
    CellN(cr3,3) = n3 
    CellN(cr3,4) = n34
    CellN(cr3,5) = nF5
    CellN(cr3,6) = n0  
    CellN(cr3,7) = n37
    CellN(cr3,8) = nF4

    ! cr4 -!
    CellN(cr4,1) = nF1
    CellN(cr4,2) = n24
    CellN(cr4,3) = n34
    CellN(cr4,4) = n4 
    CellN(cr4,5) = n0 
    CellN(cr4,6) = nF3 
    CellN(cr4,7) = nF4
    CellN(cr4,8) = n48

    ! cr5 -!
    CellN(cr5,1) = n15
    CellN(cr5,2) = nF2
    CellN(cr5,3) = nF5
    CellN(cr5,4) = n0 
    CellN(cr5,5) = n5 
    CellN(cr5,6) = n56 
    CellN(cr5,7) = n57
    CellN(cr5,8) = nF6

    ! cr6 -!
    CellN(cr6,1) = nF2
    CellN(cr6,2) = n26
    CellN(cr6,3) = n0 
    CellN(cr6,4) = nF3
    CellN(cr6,5) = n56
    CellN(cr6,6) = n6  
    CellN(cr6,7) = nF6
    CellN(cr6,8) = n68

    ! cr7 -!
    CellN(cr7,1) = nF5
    CellN(cr7,2) = n0 
    CellN(cr7,3) = n37
    CellN(cr7,4) = nF4
    CellN(cr7,5) = n57
    CellN(cr7,6) = nF6 
    CellN(cr7,7) = n7 
    CellN(cr7,8) = n78

    ! cr8 -!
    CellN(cr8,1) = n0 
    CellN(cr8,2) = nF3
    CellN(cr8,3) = nF4
    CellN(cr8,4) = n48
    CellN(cr8,5) = nF6
    CellN(cr8,6) = n68 
    CellN(cr8,7) = n78
    CellN(cr8,8) = n8 

    !-----------------------!
    !                       !
    !   Non-refined cells   !
    !                       !
    !-----------------------!
    else

      if(CelMar(c1)   >  0) then  ! neighbor 1 refined
        CellC(c, 1) = CelMar(c1) - 8 + Which_Node(c1,n1)
        CellC(c, 7) = CelMar(c1) - 8 + Which_Node(c1,n2)
        CellC(c,13) = CelMar(c1) - 8 + Which_Node(c1,n3)
        CellC(c,19) = CelMar(c1) - 8 + Which_Node(c1,n4)
      endif           

      if(CelMar(c2)   >  0) then  ! neighbor 2 refined
        CellC(c, 2) = CelMar(c2) - 8 + Which_Node(c2, n1)
        CellC(c, 8) = CelMar(c2) - 8 + Which_Node(c2, n2)
        CellC(c,14) = CelMar(c2) - 8 + Which_Node(c2, n5)
        CellC(c,20) = CelMar(c2) - 8 + Which_Node(c2, n6)
      endif           

      if(CelMar(c3)   >  0) then  ! neighbor 3 refined
        CellC(c, 3) = CelMar(c3) - 8 + Which_Node(c3, n2)
        CellC(c, 9) = CelMar(c3) - 8 + Which_Node(c3, n4)
        CellC(c,15) = CelMar(c3) - 8 + Which_Node(c3, n6)
        CellC(c,21) = CelMar(c3) - 8 + Which_Node(c3, n8)
      endif           

      if(CelMar(c4)   >  0) then  ! neighbor 4 refined
        CellC(c, 4) = CelMar(c4) - 8 + Which_Node(c4, n3)
        CellC(c,10) = CelMar(c4) - 8 + Which_Node(c4, n4)
        CellC(c,16) = CelMar(c4) - 8 + Which_Node(c4, n7)
        CellC(c,22) = CelMar(c4) - 8 + Which_Node(c4, n8)
      endif           

      if(CelMar(c5)   >  0) then  ! neighbor 5 refined
        CellC(c, 5) = CelMar(c5) - 8 + Which_Node(c5, n1)
        CellC(c,11) = CelMar(c5) - 8 + Which_Node(c5, n3)
        CellC(c,17) = CelMar(c5) - 8 + Which_Node(c5, n5)
        CellC(c,23) = CelMar(c5) - 8 + Which_Node(c5, n7)
      endif

      if(CelMar(c6)   >  0) then  ! neighbor 6 refined
        CellC(c, 6) = CelMar(c6) - 8 + Which_Node(c6, n5)
        CellC(c,12) = CelMar(c6) - 8 + Which_Node(c6, n6)
        CellC(c,18) = CelMar(c6) - 8 + Which_Node(c6, n7)
        CellC(c,24) = CelMar(c6) - 8 + Which_Node(c6, n8)
      endif           

    end if   
  end do

  write(*,*) 'Number of nodes after the refinement: ', NN 
  write(*,*) 'Number of cells after the refinement: ', NC 

  !------------------------------------------!
  !                                          !
  !   Connect the new twins, if they exist   !
  !                                          !
  !------------------------------------------!  

  do nA=1,NN2
    nA0=NodeN2(nA,0)
    nA1=NodeN2(nA,1)
    nA2=NodeN2(nA,2)

    if( (TwinN(nA1,0) /= 0).and.(TwinN(nA2,0) /= 0) ) then

      do nB=nA+1,NN2
        nB0=NodeN2(nB,0)
        nB1=NodeN2(nB,1)
        nB2=NodeN2(nB,2)

        if( (TwinN(nB1,0) /= 0).and.(TwinN(nB2,0) /= 0) ) then

          if( (Are_Nodes_Twins(nA1,nB1) .and. Are_Nodes_Twins(nA2,nB2)) .or.          &
              (Are_Nodes_Twins(nA1,nB2) .and. Are_Nodes_Twins(nA2,nB1))  ) then
            if (.not. Are_Nodes_Twins(nA0,nB0)) then
              TwinN(nA0,0)=TwinN(nA0,0)+1
              TwinN(nA0,TwinN(nA0,0))=nB0
              TwinN(nB0,0)=TwinN(nB0,0)+1
              TwinN(nB0,TwinN(nB0,0))=nA0
            end if
          end if
        end if
      end do
    end if
  end do

  do nA=1,NN4
    nA0=NodeN4(nA,0)
    nA1=NodeN4(nA,1)
    nA2=NodeN4(nA,4) ! diagonal

    if( (TwinN(nA1,0) /= 0).and.(TwinN(nA2,0) /= 0) ) then

      do nB=nA+1,NN4
        nB0=NodeN4(nB,0)
        nB1=NodeN4(nB,1)
        nB2=NodeN4(nB,4) ! diagonal

        if( (TwinN(nB1,0) /= 0).and.(TwinN(nB2,0) /= 0) ) then

          if( (Are_Nodes_Twins(nA1,nB1) .and. Are_Nodes_Twins(nA2,nB2)) .or.          &
              (Are_Nodes_Twins(nA1,nB2) .and. Are_Nodes_Twins(nA2,nB1))  ) then
            if (.not. Are_Nodes_Twins(nA0,nB0)) then
              TwinN(nA0,0)=TwinN(nA0,0)+1
              TwinN(nA0,TwinN(nA0,0))=nB0
              TwinN(nB0,0)=TwinN(nB0,0)+1
              TwinN(nB0,TwinN(nB0,0))=nA0
            end if
          end if
        end if
      end do
    end if
  end do

  !----------------------------!
  !                            !
  !   Delete redundant cells   !
  !                            !
  !----------------------------!  

  ! Initialize the new numbers for the cells
  do c=-NbC,NC
    NewN(c)=c
  end do

  del=0 
  do c=1,NC
    if(CelMar(c) /= 0) then
      NewN(c) = -1
      del = del+1
    else
      NewN(c) = c - del 
    endif 
  end do
  write(*,*) 'Deleted cells:', del

  do c=1,NC
    if(NewN(c) /= -1) then

      ! Update the cell numbers. Watch out ! The numbers you are
      ! updating are old, so double indexing is needed
      do n=1,24  ! n is neighbour now
        CellC( NewN(c),n ) = NewN(CellC( c,n ))
      end do

      ! Update the node numbers
      do n=1,8   ! n is node now
        CellN( NewN(c),n ) = CellN( c,n ) 
      end do

      material( NewN(c) ) = material( c )  ! -> never checked !
      level( NewN(c) )    = level( c )     ! -> never checked !
    end if
  end do

  do c=NC-del+1, MAXN   ! erase old data
    do n=1,24           ! n is neighbour now
      CellC(c,n ) = 0
    end do
  end do

  NC = NC - del    

  write(*,*) 'Number of cells after the renumeration: ', NC 

  end subroutine Refine_Marked_Cells