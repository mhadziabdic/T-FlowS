!==============================================================================!
  subroutine Refine_Marked_Cells(lev)
!------------------------------------------------------------------------------!
!   Refine the marked cells.                                                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!---------------------------------[Arguments]----------------------------------!
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
  integer :: nn_2, nn_4, nn_8
  integer, allocatable :: node_n2(:,:)    
  integer, allocatable :: node_n4(:,:)  
  integer, allocatable :: node_n8(:,:)    
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

  write(*,*) '# Refine: Number of nodes: ', NN 
  write(*,*) '#         Number of cells: ', NC 

  n_cells_old = NC 
  n_nodes_old = NN
  nn_2 = 0 
  nn_4 = 0
  nn_8 = 0

  allocate (node_n2(grid % max_n_nodes,0:2));   node_n2=0 
  allocate (node_n4(grid % max_n_nodes,0:4));   node_n4=0
  allocate (node_n8(grid % max_n_nodes,0:8));   node_n8=0

  !---------------------!
  !                     !
  !   Count new celss   !
  !                     !
  !---------------------!
  do c=1,n_cells_old
    if(CelMar(c) == -1) then 
      NC = NC + 8
      CelMar(c) = NC   ! now points to cr8
    end if
  end do

  do c=1,n_cells_old

    c1 = grid % cells_c(1,c)
    c2 = grid % cells_c(2,c)
    c3 = grid % cells_c(3,c)
    c4 = grid % cells_c(4,c)
    c5 = grid % cells_c(5,c)
    c6 = grid % cells_c(6,c)

    n1 = grid % cells_n(1,c)
    n2 = grid % cells_n(2,c)
    n3 = grid % cells_n(3,c)
    n4 = grid % cells_n(4,c)
    n5 = grid % cells_n(5,c)
    n6 = grid % cells_n(6,c)
    n7 = grid % cells_n(7,c)
    n8 = grid % cells_n(8,c)

    !-------------------!
    !                   !
    !   Refined cells   !
    !                   !
    !-------------------!
    if( CelMar(c) > 0 ) then  ! only refined

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
      grid % cells_c(6,cr1) = cr5
      grid % cells_c(6,cr2) = cr6
      grid % cells_c(6,cr3) = cr7
      grid % cells_c(6,cr4) = cr8

      ! 5
      grid % cells_c(5,cr2) = cr1
      grid % cells_c(5,cr4) = cr3
      grid % cells_c(5,cr6) = cr5
      grid % cells_c(5,cr8) = cr7

      ! 4
      grid % cells_c(4,cr1) = cr3
      grid % cells_c(4,cr2) = cr4
      grid % cells_c(4,cr5) = cr7
      grid % cells_c(4,cr6) = cr8

      ! 3
      grid % cells_c(3,cr1) = cr2
      grid % cells_c(3,cr3) = cr4
      grid % cells_c(3,cr5) = cr6
      grid % cells_c(3,cr7) = cr8

      ! 2
      grid % cells_c(2,cr3) = cr1
      grid % cells_c(2,cr4) = cr2
      grid % cells_c(2,cr7) = cr5
      grid % cells_c(2,cr8) = cr6

      ! 1
      grid % cells_c(1,cr5) = cr1
      grid % cells_c(1,cr6) = cr2
      grid % cells_c(1,cr7) = cr3
      grid % cells_c(1,cr8) = cr4

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
      if(CelMar(c1) == 0) then  ! neighbor 1 not refined
        grid % cells_c(1,cr1) = c1
        grid % cells_c(1,cr2) = c1
        grid % cells_c(1,cr3) = c1
        grid % cells_c(1,cr4) = c1
      else                        ! neighbor 1 refined
        grid % cells_c(1,cr1) = CelMar(c1) - 8 + Which_Node(c1,n1)
        grid % cells_c(1,cr2) = CelMar(c1) - 8 + Which_Node(c1,n2)
        grid % cells_c(1,cr3) = CelMar(c1) - 8 + Which_Node(c1,n3)
        grid % cells_c(1,cr4) = CelMar(c1) - 8 + Which_Node(c1,n4)
      endif           

      if(CelMar(c2) == 0) then  ! neighbor 2 not refined
        grid % cells_c(2,cr1) = c2
        grid % cells_c(2,cr2) = c2
        grid % cells_c(2,cr5) = c2
        grid % cells_c(2,cr6) = c2
      else                        ! neighbor 2 refined
        grid % cells_c(2,cr1) = CelMar(c2) - 8 + Which_Node(c2, n1)
        grid % cells_c(2,cr2) = CelMar(c2) - 8 + Which_Node(c2, n2)
        grid % cells_c(2,cr5) = CelMar(c2) - 8 + Which_Node(c2, n5)
        grid % cells_c(2,cr6) = CelMar(c2) - 8 + Which_Node(c2, n6)
      endif           

      if(CelMar(c3) == 0) then  ! neighbor 3 not refined
        grid % cells_c(3,cr2) = c3
        grid % cells_c(3,cr4) = c3
        grid % cells_c(3,cr6) = c3
        grid % cells_c(3,cr8) = c3
      else                        ! neighbor 3 refined
        grid % cells_c(3,cr2) = CelMar(c3) - 8 + Which_Node(c3, n2)
        grid % cells_c(3,cr4) = CelMar(c3) - 8 + Which_Node(c3, n4)
        grid % cells_c(3,cr6) = CelMar(c3) - 8 + Which_Node(c3, n6)
        grid % cells_c(3,cr8) = CelMar(c3) - 8 + Which_Node(c3, n8)
      endif           

      if(CelMar(c4) == 0) then  ! neighbor 4 not refine
        grid % cells_c(4,cr3) = c4
        grid % cells_c(4,cr4) = c4
        grid % cells_c(4,cr7) = c4
        grid % cells_c(4,cr8) = c4
      else                        ! neighbor 4 refine
        grid % cells_c(4,cr3) = CelMar(c4) - 8 + Which_Node(c4, n3)
        grid % cells_c(4,cr4) = CelMar(c4) - 8 + Which_Node(c4, n4)
        grid % cells_c(4,cr7) = CelMar(c4) - 8 + Which_Node(c4, n7)
        grid % cells_c(4,cr8) = CelMar(c4) - 8 + Which_Node(c4, n8)
      endif           

      if(CelMar(c5) == 0) then  ! neighbor 5 not refined
        grid % cells_c(5,cr1) = c5
        grid % cells_c(5,cr3) = c5
        grid % cells_c(5,cr5) = c5
        grid % cells_c(5,cr7) = c5
      else                        ! neighbor 5 refined
        grid % cells_c(5,cr1) = CelMar(c5) - 8 + Which_Node(c5, n1)
        grid % cells_c(5,cr3) = CelMar(c5) - 8 + Which_Node(c5, n3)
        grid % cells_c(5,cr5) = CelMar(c5) - 8 + Which_Node(c5, n5)
        grid % cells_c(5,cr7) = CelMar(c5) - 8 + Which_Node(c5, n7)
      endif

      if(CelMar(c6) == 0) then  ! neighbor 6 not refined
        grid % cells_c(6,cr5) = c6
        grid % cells_c(6,cr6) = c6
        grid % cells_c(6,cr7) = c6
        grid % cells_c(6,cr8) = c6
      else                        ! neighbor 6 refined
        grid % cells_c(6,cr5) = CelMar(c6) - 8 + Which_Node(c6, n5)
        grid % cells_c(6,cr6) = CelMar(c6) - 8 + Which_Node(c6, n6)
        grid % cells_c(6,cr7) = CelMar(c6) - 8 + Which_Node(c6, n7)
        grid % cells_c(6,cr8) = CelMar(c6) - 8 + Which_Node(c6, n8)
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
    do n=1,nn_2
      if( ( node_n2(n,1) == n1 .and. node_n2(n,2) == n2 ) .or.  &
          ( node_n2(n,1) == n2 .and. node_n2(n,2) == n1) ) then
        n12 = node_n2(n,0) 
      end if
    end do
    if (n12 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n12 = NN
      node_n2(nn_2,0) = n12
      node_n2(nn_2,1) = n1
      node_n2(nn_2,2) = n2
      grid % nodes(n12) % x = .5 * (grid % nodes(n1) % x + grid % nodes(n2) % x)
      grid % nodes(n12) % y = .5 * (grid % nodes(n1) % y + grid % nodes(n2) % y)
      grid % nodes(n12) % z = .5 * (grid % nodes(n1) % z + grid % nodes(n2) % z)
    end if 

    ! n13
    do n=1,nn_2
      if( ( node_n2(n,1) == n1 .and. node_n2(n,2) == n3 ) .or.  &
          ( node_n2(n,1) == n3 .and. node_n2(n,2) == n1) ) then
        n13 = node_n2(n,0) 
      end if
    end do
    if (n13 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n13 = NN
      node_n2(nn_2,0) = n13
      node_n2(nn_2,1) = n1
      node_n2(nn_2,2) = n3
      grid % nodes(n13) % x = .5 * (grid % nodes(n1) % x + grid % nodes(n3) % x)
      grid % nodes(n13) % y = .5 * (grid % nodes(n1) % y + grid % nodes(n3) % y)
      grid % nodes(n13) % z = .5 * (grid % nodes(n1) % z + grid % nodes(n3) % z)
    end if 

    ! n24
    do n=1,nn_2
      if( ( node_n2(n,1) == n2 .and. node_n2(n,2) == n4 ) .or.  &
          ( node_n2(n,1) == n4 .and. node_n2(n,2) == n2) ) then
        n24 = node_n2(n,0) 
      end if
    end do
    if (n24 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n24 = NN
      node_n2(nn_2,0) = n24
      node_n2(nn_2,1) = n2
      node_n2(nn_2,2) = n4
      grid % nodes(n24) % x = .5 * (grid % nodes(n2) % x + grid % nodes(n4) % x)
      grid % nodes(n24) % y = .5 * (grid % nodes(n2) % y + grid % nodes(n4) % y)
      grid % nodes(n24) % z = .5 * (grid % nodes(n2) % z + grid % nodes(n4) % z)
    end if 

    ! n34
    do n=1,nn_2
      if( ( node_n2(n,1) == n3 .and. node_n2(n,2) == n4 ) .or.  &
          ( node_n2(n,1) == n4 .and. node_n2(n,2) == n3) ) then
        n34 = node_n2(n,0) 
      end if
    end do
    if (n34 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n34 = NN
      node_n2(nn_2,0) = n34
      node_n2(nn_2,1) = n3
      node_n2(nn_2,2) = n4
      grid % nodes(n34) % x = .5 * (grid % nodes(n3) % x + grid % nodes(n4) % x)
      grid % nodes(n34) % y = .5 * (grid % nodes(n3) % y + grid % nodes(n4) % y)
      grid % nodes(n34) % z = .5 * (grid % nodes(n3) % z + grid % nodes(n4) % z)
    end if 

    ! n15
    do n=1,nn_2
      if( ( node_n2(n,1) == n1 .and. node_n2(n,2) == n5 ) .or.  &
          ( node_n2(n,1) == n5 .and. node_n2(n,2) == n1) ) then
        n15 = node_n2(n,0) 
      end if
    end do
    if (n15 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n15 = NN
      node_n2(nn_2,0) = n15
      node_n2(nn_2,1) = n1
      node_n2(nn_2,2) = n5
      grid % nodes(n15) % x = .5 * (grid % nodes(n1) % x + grid % nodes(n5) % x)
      grid % nodes(n15) % y = .5 * (grid % nodes(n1) % y + grid % nodes(n5) % y)
      grid % nodes(n15) % z = .5 * (grid % nodes(n1) % z + grid % nodes(n5) % z)
    end if 

    ! n26
    do n=1,nn_2
      if( ( node_n2(n,1) == n2 .and. node_n2(n,2) == n6 ) .or.  &
          ( node_n2(n,1) == n6 .and. node_n2(n,2) == n2) ) then
        n26 = node_n2(n,0) 
      end if
    end do
    if (n26 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n26 = NN
      node_n2(nn_2,0) = n26
      node_n2(nn_2,1) = n2
      node_n2(nn_2,2) = n6
      grid % nodes(n26) % x = .5 * (grid % nodes(n2) % x + grid % nodes(n6) % x)
      grid % nodes(n26) % y = .5 * (grid % nodes(n2) % y + grid % nodes(n6) % y)
      grid % nodes(n26) % z = .5 * (grid % nodes(n2) % z + grid % nodes(n6) % z)
    end if 

    ! n37
    do n=1,nn_2
      if( ( node_n2(n,1) == n3 .and. node_n2(n,2) == n7 ) .or.  &
          ( node_n2(n,1) == n7 .and. node_n2(n,2) == n3) ) then
        n37 = node_n2(n,0) 
      end if
    end do
    if (n37 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n37 = NN
      node_n2(nn_2,0) = n37
      node_n2(nn_2,1) = n3
      node_n2(nn_2,2) = n7
      grid % nodes(n37) % x = .5 * (grid % nodes(n3) % x + grid % nodes(n7) % x)
      grid % nodes(n37) % y = .5 * (grid % nodes(n3) % y + grid % nodes(n7) % y)
      grid % nodes(n37) % z = .5 * (grid % nodes(n3) % z + grid % nodes(n7) % z)
    end if 

    ! n48
    do n=1,nn_2
      if( ( node_n2(n,1) == n4 .and. node_n2(n,2) == n8 ) .or.  &
          ( node_n2(n,1) == n8 .and. node_n2(n,2) == n4) ) then
        n48 = node_n2(n,0) 
      end if
    end do
    if (n48 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n48 = NN
      node_n2(nn_2,0) = n48
      node_n2(nn_2,1) = n4
      node_n2(nn_2,2) = n8
      grid % nodes(n48) % x = .5 * (grid % nodes(n4) % x + grid % nodes(n8) % x)
      grid % nodes(n48) % y = .5 * (grid % nodes(n4) % y + grid % nodes(n8) % y)
      grid % nodes(n48) % z = .5 * (grid % nodes(n4) % z + grid % nodes(n8) % z)
    end if 

    ! n56
    do n=1,nn_2
      if( ( node_n2(n,1) == n5 .and. node_n2(n,2) == n6 ) .or.  &
          ( node_n2(n,1) == n6 .and. node_n2(n,2) == n5) ) then
        n56 = node_n2(n,0) 
      end if
    end do
    if (n56 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n56 = NN
      node_n2(nn_2,0) = n56
      node_n2(nn_2,1) = n5
      node_n2(nn_2,2) = n6
      grid % nodes(n56) % x = .5 * (grid % nodes(n5) % x + grid % nodes(n6) % x)
      grid % nodes(n56) % y = .5 * (grid % nodes(n5) % y + grid % nodes(n6) % y)
      grid % nodes(n56) % z = .5 * (grid % nodes(n5) % z + grid % nodes(n6) % z)
    end if 

    ! n57
    do n=1,nn_2
      if( ( node_n2(n,1) == n5 .and. node_n2(n,2) == n7 ) .or.  &
          ( node_n2(n,1) == n7 .and. node_n2(n,2) == n5) ) then
        n57 = node_n2(n,0) 
      end if
    end do
    if (n57 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n57 = NN
      node_n2(nn_2,0) = n57
      node_n2(nn_2,1) = n5
      node_n2(nn_2,2) = n7
      grid % nodes(n57) % x = .5 * (grid % nodes(n5) % x + grid % nodes(n7) % x)
      grid % nodes(n57) % y = .5 * (grid % nodes(n5) % y + grid % nodes(n7) % y)
      grid % nodes(n57) % z = .5 * (grid % nodes(n5) % z + grid % nodes(n7) % z)
    end if 

    ! n68 
    do n=1,nn_2
      if( ( node_n2(n,1) == n6 .and. node_n2(n,2) == n8 ) .or.  &
          ( node_n2(n,1) == n8 .and. node_n2(n,2) == n6) ) then
        n68 = node_n2(n,0) 
      end if
    end do
    if (n68 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n68 = NN
      node_n2(nn_2,0) = n68
      node_n2(nn_2,1) = n6
      node_n2(nn_2,2) = n8
      grid % nodes(n68) % x = .5 * (grid % nodes(n6) % x + grid % nodes(n8) % x)
      grid % nodes(n68) % y = .5 * (grid % nodes(n6) % y + grid % nodes(n8) % y)
      grid % nodes(n68) % z = .5 * (grid % nodes(n6) % z + grid % nodes(n8) % z)
    end if 

    ! n78
    do n=1,nn_2
      if( ( node_n2(n,1) == n7 .and. node_n2(n,2) == n8 ) .or.  &
          ( node_n2(n,1) == n8 .and. node_n2(n,2) == n7) ) then
        n78 = node_n2(n,0) 
      end if
    end do
    if (n78 == 0) then
      nn_2 = nn_2 + 1
      NN  = NN  + 1
      n78 = NN
      node_n2(nn_2,0) = n78
      node_n2(nn_2,1) = n7
      node_n2(nn_2,2) = n8
      grid % nodes(n78) % x = .5 * (grid % nodes(n7) % x + grid % nodes(n8) % x)
      grid % nodes(n78) % y = .5 * (grid % nodes(n7) % y + grid % nodes(n8) % y)
      grid % nodes(n78) % z = .5 * (grid % nodes(n7) % z + grid % nodes(n8) % z)
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
    do n=1,nn_4
      if( ( node_n4(n,1) == n1 .and. node_n4(n,4) == n4 ) .or.  &
          ( node_n4(n,1) == n4 .and. node_n4(n,4) == n1 ) .or.  &
          ( node_n4(n,1) == n2 .and. node_n4(n,4) == n3 ) .or.  &
          ( node_n4(n,1) == n3 .and. node_n4(n,4) == n2 ) ) then
        nF1 = node_n4(n,0) 
      end if
    end do
    if (nF1 == 0) then
      nn_4 = nn_4 + 1
      NN  = NN  + 1
      nF1 = NN
      node_n4(nn_4,0) = nF1
      node_n4(nn_4,1) = n1
      node_n4(nn_4,2) = n2
      node_n4(nn_4,3) = n3
      node_n4(nn_4,4) = n4
      grid % nodes(nf1) % x = 0.25 * (                                        &
                               grid % nodes(n1) % x + grid % nodes(n2) % x +  &
                               grid % nodes(n3) % x + grid % nodes(n4) % x)
      grid % nodes(nf1) % y = 0.25 * (                                        &
                               grid % nodes(n1) % y + grid % nodes(n2) % y +  &
                               grid % nodes(n3) % y + grid % nodes(n4) % y)
      grid % nodes(nf1) % z = 0.25 * (                                        &
                               grid % nodes(n1) % z + grid % nodes(n2) % z +  &
                               grid % nodes(n3) % z + grid % nodes(n4) % z)
    end if 

    ! nF2
    do n=1,nn_4
      if( ( node_n4(n,1) == n1 .and. node_n4(n,4) == n6 ) .or.  &
          ( node_n4(n,1) == n6 .and. node_n4(n,4) == n1 ) .or.  &
          ( node_n4(n,1) == n2 .and. node_n4(n,4) == n5 ) .or.  &
          ( node_n4(n,1) == n5 .and. node_n4(n,4) == n2 ) ) then
        nF2 = node_n4(n,0) 
      end if
    end do
    if (nF2 == 0) then
      nn_4 = nn_4 + 1
      NN  = NN  + 1
      nF2 = NN
      node_n4(nn_4,0) = nF2
      node_n4(nn_4,1) = n1
      node_n4(nn_4,2) = n2
      node_n4(nn_4,3) = n5
      node_n4(nn_4,4) = n6
      grid % nodes(nf2) % x = 0.25 * (                                        &
                               grid % nodes(n1) % x + grid % nodes(n2) % x +  &
                               grid % nodes(n5) % x + grid % nodes(n6) % x)
      grid % nodes(nf2) % y = 0.25 * (                                        &
                               grid % nodes(n1) % y + grid % nodes(n2) % y +  &
                               grid % nodes(n5) % y + grid % nodes(n6) % y)
      grid % nodes(nf2) % z = 0.25 * (                                        &
                               grid % nodes(n1) % z + grid % nodes(n2) % z +  &
                               grid % nodes(n5) % z + grid % nodes(n6) % z)
    end if 

    ! nF3
    do n=1,nn_4
      if( ( node_n4(n,1) == n2 .and. node_n4(n,4) == n8 ) .or.  &
          ( node_n4(n,1) == n8 .and. node_n4(n,4) == n2 ) .or.  &
          ( node_n4(n,1) == n4 .and. node_n4(n,4) == n6 ) .or.  &
          ( node_n4(n,1) == n6 .and. node_n4(n,4) == n4 ) ) then
        nF3 = node_n4(n,0) 
      end if
    end do
    if (nF3 == 0) then
      nn_4 = nn_4 + 1
      NN  = NN  + 1
      nF3 = NN
      node_n4(nn_4,0) = nF3
      node_n4(nn_4,1) = n2
      node_n4(nn_4,2) = n4
      node_n4(nn_4,3) = n6
      node_n4(nn_4,4) = n8
      grid % nodes(nf3) % x = 0.25 * (                                        &
                               grid % nodes(n2) % x + grid % nodes(n4) % x +  &
                               grid % nodes(n6) % x + grid % nodes(n8) % x)
      grid % nodes(nf3) % y = 0.25 * (                                        &
                               grid % nodes(n2) % y + grid % nodes(n4) % y +  &
                               grid % nodes(n6) % y + grid % nodes(n8) % y)
      grid % nodes(nf3) % z = 0.25 * (                                        &
                               grid % nodes(n2) % z + grid % nodes(n4) % z +  &
                               grid % nodes(n6) % z + grid % nodes(n8) % z)
    end if 

    ! nF4
    do n=1,nn_4
      if( ( node_n4(n,1) == n3 .and. node_n4(n,4) == n8 ) .or.  &
          ( node_n4(n,1) == n8 .and. node_n4(n,4) == n3 ) .or.  &
          ( node_n4(n,1) == n4 .and. node_n4(n,4) == n7 ) .or.  &
          ( node_n4(n,1) == n7 .and. node_n4(n,4) == n4 ) ) then
        nF4 = node_n4(n,0) 
      end if
    end do
    if (nF4 == 0) then
      nn_4 = nn_4 + 1
      NN  = NN  + 1
      nF4 = NN
      node_n4(nn_4,0) = nF4
      node_n4(nn_4,1) = n3
      node_n4(nn_4,2) = n4
      node_n4(nn_4,3) = n7
      node_n4(nn_4,4) = n8
      grid % nodes(nf4) % x = 0.25 * (                                        &
                               grid % nodes(n3) % x + grid % nodes(n4) % x +  &
                               grid % nodes(n7) % x + grid % nodes(n8) % x)
      grid % nodes(nf4) % y = 0.25 * (                                        &
                               grid % nodes(n3) % y + grid % nodes(n4) % y +  &
                               grid % nodes(n7) % y + grid % nodes(n8) % y)
      grid % nodes(nf4) % z = 0.25 * (                                        &
                               grid % nodes(n3) % z + grid % nodes(n4) % z +  &
                               grid % nodes(n7) % z + grid % nodes(n8) % z)
    end if 

    ! nF5
    do n=1,nn_4
      if( ( node_n4(n,1) == n1 .and. node_n4(n,4) == n7 ) .or.  &
          ( node_n4(n,1) == n7 .and. node_n4(n,4) == n1 ) .or.  &
          ( node_n4(n,1) == n3 .and. node_n4(n,4) == n5 ) .or.  &
          ( node_n4(n,1) == n5 .and. node_n4(n,4) == n3 ) ) then
        nF5 = node_n4(n,0) 
      end if
    end do
    if (nF5 == 0) then
      nn_4 = nn_4 + 1
      NN  = NN  + 1
      nF5 = NN
      node_n4(nn_4,0) = nF5
      node_n4(nn_4,1) = n1
      node_n4(nn_4,2) = n3
      node_n4(nn_4,3) = n5
      node_n4(nn_4,4) = n7
      grid % nodes(nf5) % x = 0.25 * (                                        &
                               grid % nodes(n1) % x + grid % nodes(n3) % x +  &
                               grid % nodes(n5) % x + grid % nodes(n7) % x)
      grid % nodes(nf5) % y = 0.25 * (                                        &
                               grid % nodes(n1) % y + grid % nodes(n3) % y +  &
                               grid % nodes(n5) % y + grid % nodes(n7) % y)
      grid % nodes(nf5) % z = 0.25 * (                                        &
                               grid % nodes(n1) % z + grid % nodes(n3) % z +  &
                               grid % nodes(n5) % z + grid % nodes(n7) % z)
    end if 

    ! nF6
    do n=1,nn_4
      if( ( node_n4(n,1) == n5 .and. node_n4(n,4) == n8 ) .or.  &
          ( node_n4(n,1) == n8 .and. node_n4(n,4) == n5 ) .or.  &
          ( node_n4(n,1) == n6 .and. node_n4(n,4) == n7 ) .or.  &
          ( node_n4(n,1) == n7 .and. node_n4(n,4) == n6 ) ) then
        nF6 = node_n4(n,0) 
      end if
    end do
    if (nF6 == 0) then
      nn_4 = nn_4 + 1
      NN  = NN  + 1
      nF6 = NN
      node_n4(nn_4,0) = nF6
      node_n4(nn_4,1) = n5
      node_n4(nn_4,2) = n6
      node_n4(nn_4,3) = n7
      node_n4(nn_4,4) = n8
      grid % nodes(nf6) % x = 0.25 * (                                        &
                               grid % nodes(n5) % x + grid % nodes(n6) % x +  &
                               grid % nodes(n7) % x + grid % nodes(n8) % x)
      grid % nodes(nf6) % y = 0.25 * (                                        &
                               grid % nodes(n5) % y + grid % nodes(n6) % y +  &
                               grid % nodes(n7) % y + grid % nodes(n8) % y)
      grid % nodes(nf6) % z = 0.25 * (                                        &
                               grid % nodes(n5) % z + grid % nodes(n6) % z +  &
                               grid % nodes(n7) % z + grid % nodes(n8) % z)
    end if 

    !----------------------------------------!
    !   Eventually, the node in the middle   !
    !----------------------------------------!
    nn_8 = nn_8 + 1
    NN  = NN + 1
    n0  = NN
    node_n8(nn_8,0) = n0 
    node_n8(nn_8,1) = n1
    node_n8(nn_8,2) = n2
    node_n8(nn_8,3) = n3
    node_n8(nn_8,4) = n4
    node_n8(nn_8,5) = n5
    node_n8(nn_8,6) = n6
    node_n8(nn_8,7) = n7
    node_n8(nn_8,8) = n8
    grid % nodes(n0) % x = .125*(grid % nodes(n1) % x + grid % nodes(n2) % x + &
                                 grid % nodes(n3) % x + grid % nodes(n4) % x + &
                                 grid % nodes(n5) % x + grid % nodes(n6) % x + &
                                 grid % nodes(n7) % x + grid % nodes(n8) % x)
    grid % nodes(n0) % y = .125*(grid % nodes(n1) % y + grid % nodes(n2) % y + &
                                 grid % nodes(n3) % y + grid % nodes(n4) % y + &
                                 grid % nodes(n5) % y + grid % nodes(n6) % y + &
                                 grid % nodes(n7) % y + grid % nodes(n8) % y)
    grid % nodes(n0) % z = .125*(grid % nodes(n1) % z + grid % nodes(n2) % z + &
                                 grid % nodes(n3) % z + grid % nodes(n4) % z + &
                                 grid % nodes(n5) % z + grid % nodes(n6) % z + &
                                 grid % nodes(n7) % z + grid % nodes(n8) % z)

    !----------------------------!
    !   Set nodes to new cells   !
    !----------------------------!

    ! cr1 -!
    grid % cells_n(1,cr1) = n1
    grid % cells_n(2,cr1) = n12
    grid % cells_n(3,cr1) = n13
    grid % cells_n(4,cr1) = nF1
    grid % cells_n(5,cr1) = n15
    grid % cells_n(6,cr1) = nF2
    grid % cells_n(7,cr1) = nF5
    grid % cells_n(8,cr1) = n0 

    ! cr2 -!
    grid % cells_n(1,cr2) = n12
    grid % cells_n(2,cr2) = n2 
    grid % cells_n(3,cr2) = nF1
    grid % cells_n(4,cr2) = n24
    grid % cells_n(5,cr2) = nF2
    grid % cells_n(6,cr2) = n26 
    grid % cells_n(7,cr2) = n0 
    grid % cells_n(8,cr2) = nF3

    ! cr3 -!
    grid % cells_n(1,cr3) = n13
    grid % cells_n(2,cr3) = nF1
    grid % cells_n(3,cr3) = n3 
    grid % cells_n(4,cr3) = n34
    grid % cells_n(5,cr3) = nF5
    grid % cells_n(6,cr3) = n0  
    grid % cells_n(7,cr3) = n37
    grid % cells_n(8,cr3) = nF4

    ! cr4 -!
    grid % cells_n(1,cr4) = nF1
    grid % cells_n(2,cr4) = n24
    grid % cells_n(3,cr4) = n34
    grid % cells_n(4,cr4) = n4 
    grid % cells_n(5,cr4) = n0 
    grid % cells_n(6,cr4) = nF3 
    grid % cells_n(7,cr4) = nF4
    grid % cells_n(8,cr4) = n48

    ! cr5 -!
    grid % cells_n(1,cr5) = n15
    grid % cells_n(2,cr5) = nF2
    grid % cells_n(3,cr5) = nF5
    grid % cells_n(4,cr5) = n0 
    grid % cells_n(5,cr5) = n5 
    grid % cells_n(6,cr5) = n56 
    grid % cells_n(7,cr5) = n57
    grid % cells_n(8,cr5) = nF6

    ! cr6 -!
    grid % cells_n(1,cr6) = nF2
    grid % cells_n(2,cr6) = n26
    grid % cells_n(3,cr6) = n0 
    grid % cells_n(4,cr6) = nF3
    grid % cells_n(5,cr6) = n56
    grid % cells_n(6,cr6) = n6  
    grid % cells_n(7,cr6) = nF6
    grid % cells_n(8,cr6) = n68

    ! cr7 -!
    grid % cells_n(1,cr7) = nF5
    grid % cells_n(2,cr7) = n0 
    grid % cells_n(3,cr7) = n37
    grid % cells_n(4,cr7) = nF4
    grid % cells_n(5,cr7) = n57
    grid % cells_n(6,cr7) = nF6 
    grid % cells_n(7,cr7) = n7 
    grid % cells_n(8,cr7) = n78

    ! cr8 -!
    grid % cells_n(1,cr8) = n0 
    grid % cells_n(2,cr8) = nF3
    grid % cells_n(3,cr8) = nF4
    grid % cells_n(4,cr8) = n48
    grid % cells_n(5,cr8) = nF6
    grid % cells_n(6,cr8) = n68 
    grid % cells_n(7,cr8) = n78
    grid % cells_n(8,cr8) = n8 

    !-----------------------!
    !                       !
    !   Non-refined cells   !
    !                       !
    !-----------------------!
    else

      if(CelMar(c1) > 0) then  ! neighbor 1 refined
        grid % cells_c( 1,c) = CelMar(c1) - 8 + Which_Node(c1,n1)
        grid % cells_c( 7,c) = CelMar(c1) - 8 + Which_Node(c1,n2)
        grid % cells_c(13,c) = CelMar(c1) - 8 + Which_Node(c1,n3)
        grid % cells_c(19,c) = CelMar(c1) - 8 + Which_Node(c1,n4)
      endif         

      if(CelMar(c2) > 0) then  ! neighbor 2 refined
        grid % cells_c( 2,c) = CelMar(c2) - 8 + Which_Node(c2, n1)
        grid % cells_c( 8,c) = CelMar(c2) - 8 + Which_Node(c2, n2)
        grid % cells_c(14,c) = CelMar(c2) - 8 + Which_Node(c2, n5)
        grid % cells_c(20,c) = CelMar(c2) - 8 + Which_Node(c2, n6)
      endif           

      if(CelMar(c3) > 0) then  ! neighbor 3 refined
        grid % cells_c( 3,c) = CelMar(c3) - 8 + Which_Node(c3, n2)
        grid % cells_c( 9,c) = CelMar(c3) - 8 + Which_Node(c3, n4)
        grid % cells_c(15,c) = CelMar(c3) - 8 + Which_Node(c3, n6)
        grid % cells_c(21,c) = CelMar(c3) - 8 + Which_Node(c3, n8)
      endif           

      if(CelMar(c4) > 0) then  ! neighbor 4 refined
        grid % cells_c( 4,c) = CelMar(c4) - 8 + Which_Node(c4, n3)
        grid % cells_c(10,c) = CelMar(c4) - 8 + Which_Node(c4, n4)
        grid % cells_c(16,c) = CelMar(c4) - 8 + Which_Node(c4, n7)
        grid % cells_c(22,c) = CelMar(c4) - 8 + Which_Node(c4, n8)
      endif           

      if(CelMar(c5) > 0) then  ! neighbor 5 refined
        grid % cells_c( 5,c) = CelMar(c5) - 8 + Which_Node(c5, n1)
        grid % cells_c(11,c) = CelMar(c5) - 8 + Which_Node(c5, n3)
        grid % cells_c(17,c) = CelMar(c5) - 8 + Which_Node(c5, n5)
        grid % cells_c(23,c) = CelMar(c5) - 8 + Which_Node(c5, n7)
      endif

      if(CelMar(c6) > 0) then  ! neighbor 6 refined
        grid % cells_c( 6,c) = CelMar(c6) - 8 + Which_Node(c6, n5)
        grid % cells_c(12,c) = CelMar(c6) - 8 + Which_Node(c6, n6)
        grid % cells_c(18,c) = CelMar(c6) - 8 + Which_Node(c6, n7)
        grid % cells_c(24,c) = CelMar(c6) - 8 + Which_Node(c6, n8)
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

  do nA=1,nn_2
    nA0=node_n2(nA,0)
    nA1=node_n2(nA,1)
    nA2=node_n2(nA,2)

    if( (TwinN(nA1,0) /= 0).and.(TwinN(nA2,0) /= 0) ) then

      do nB=nA+1,nn_2
        nB0=node_n2(nB,0)
        nB1=node_n2(nB,1)
        nB2=node_n2(nB,2)

        if( (TwinN(nB1,0) /= 0).and.(TwinN(nB2,0) /= 0) ) then

          if( (Are_Nodes_Twins(nA1,nB1) .and.   &
               Are_Nodes_Twins(nA2,nB2)) .or.   &
              (Are_Nodes_Twins(nA1,nB2) .and.   &
               Are_Nodes_Twins(nA2,nB1))  ) then
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

  do nA=1,nn_4
    nA0=node_n4(nA,0)
    nA1=node_n4(nA,1)
    nA2=node_n4(nA,4) ! diagonal

    if( (TwinN(nA1,0) /= 0).and.(TwinN(nA2,0) /= 0) ) then

      do nB=nA+1,nn_4
        nB0=node_n4(nB,0)
        nB1=node_n4(nB,1)
        nB2=node_n4(nB,4) ! diagonal

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
        grid % cells_c( n, NewN(c) ) = NewN( grid % cells_c( n, c ) )
      end do

      ! Update the node numbers
      do n=1,8   ! n is node now
        grid % cells_n( n, NewN(c) ) = grid % cells_n( n, c )
      end do

      material( NewN(c) ) = material( c )  ! -> never checked !
      level( NewN(c) )    = level( c )     ! -> never checked !
    end if
  end do

  do c=NC-del+1, grid % max_n_nodes   ! erase old data
    do n=1,24                         ! n is neighbour now
      grid % cells_c( n, c ) = 0
    end do
  end do

  NC = NC - del    

  write(*,*) 'Number of cells after the renumeration: ', NC 

  deallocate(node_n2)
  deallocate(node_n4)
  deallocate(node_n8)

  end subroutine Refine_Marked_Cells
