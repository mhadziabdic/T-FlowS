!==============================================================================!
  subroutine Grid_Topology
!------------------------------------------------------------------------------!
!   Determines the topology of the grid.                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use gen_mod 
  use neu_mod 
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j
!==============================================================================!

!------------------------------------------------------------------------------!
!
!    *.NEU nodes are numbered as this:
!
!        7-----------8
!       /|          /|
!      /           / |
!     /  |        /  |
!    5-----------6   |
!    |   |       |   |
!    |   3- - - -|- -4
!    |  /        |  /
!    |           | /
!    |/          |/
!    1-----------2                 
!
!    Figure 1: Numbering of cell nodes 
!
!        7-----------8             7-----------8 
!       /|          /|            /|          /| 
!      /           / |           /    (6)    / | 
!     /  |    (3) /  |          /  |        /  | 
!    5-----------6   |         5-----------6   | 
!    |(4)|       |(2)|         |   |       |   | 
!    |   3- - - -|- -4         |   3- - - -|- -4 
!    |  / (1)    |  /          |  /        |  / 
!    |           | /           |      (5)  | / 
!    |/          |/            |/          |/ 
!    1-----------2             1-----------2 
!
!    Figure 2: Numbering of directions of boundary  
!
! The directions for boundary conditions are:
!
!   direction:    nodes (not in anticlockwise order):
!   ------------------------------------------------- 
!      (1)   ->   1,2,5,6
!      (2)   ->   2,4,6,8
!      (3)   ->   3,4,7,8
!      (4)   ->   1,3,5,7
!      (5)   ->   1,2,3,4
!      (6)   ->   5,6,7,8
!
!------------------------------------------------------------------------------!

  !------------------------------!
  !   Count the boundary cells   !
  !------------------------------!
  NBc = 0
  NS  = 0
  do i=1,NC
    do j=1,6
      if(BCtype(i,j) /= 0) then
        NBc = NBc + 1 

        ! BCmark
        BCmark(-NBc) = BCtype(i,j)

        ! Material:
        material(-NbC) = material(i)

        ! Sides
        NS  = NS  + 1
        SideC(1,NS) = i
        SideC(2,NS) = -NbC

        ! Hexahedra:
        if(grid % cells_n_nodes(i) == 8) then
          grid % cells_n_nodes(-NbC) = 4 
          grid % cells_n(1,-NbC) = grid % cells_n(f8n(j,1), i)
          grid % cells_n(2,-NbC) = grid % cells_n(f8n(j,2), i)
          grid % cells_n(3,-NbC) = grid % cells_n(f8n(j,3), i)
          grid % cells_n(4,-NbC) = grid % cells_n(f8n(j,4), i)

          grid % faces_n_nodes(NS) = 4
          grid % faces_n(1,NS) = grid % cells_n(1,-NbC)
          grid % faces_n(2,NS) = grid % cells_n(2,-NbC)
          grid % faces_n(3,NS) = grid % cells_n(3,-NbC)
          grid % faces_n(4,NS) = grid % cells_n(4,-NbC)
        end if

        ! Prisms:
        if(grid % cells_n_nodes(i) == 6) then
          if(j <= 3) then    ! faces (1), (2) and (3)
            grid % cells_n_nodes(-NbC) = 4 
            grid % cells_n(1,-NbC) = grid % cells_n(f6n(j,1), i)
            grid % cells_n(2,-NbC) = grid % cells_n(f6n(j,2), i)
            grid % cells_n(3,-NbC) = grid % cells_n(f6n(j,3), i)
            grid % cells_n(4,-NbC) = grid % cells_n(f6n(j,4), i)

            grid % faces_n_nodes(NS) = 4
            grid % faces_n(1,NS) = grid % cells_n(1,-NbC)
            grid % faces_n(2,NS) = grid % cells_n(2,-NbC)
            grid % faces_n(3,NS) = grid % cells_n(3,-NbC)
            grid % faces_n(4,NS) = grid % cells_n(4,-NbC)
          else if(j <= 5) then
            grid % cells_n_nodes(-NbC) = 3 
            grid % cells_n(1,-NbC) = grid % cells_n(f6n(j,1), i)
            grid % cells_n(2,-NbC) = grid % cells_n(f6n(j,2), i)
            grid % cells_n(3,-NbC) = grid % cells_n(f6n(j,3), i)

            grid % faces_n_nodes(NS) = 3
            grid % faces_n(1,NS) = grid % cells_n(1,-NbC)
            grid % faces_n(2,NS) = grid % cells_n(2,-NbC)
            grid % faces_n(3,NS) = grid % cells_n(3,-NbC)
          end if
        end if

        ! Tetrahedra:
        if(grid % cells_n_nodes(i) == 4) then
          if(j <= 4) then
            grid % cells_n_nodes(-NbC) = 3 
            grid % cells_n(1,-NbC) = grid % cells_n(f4n(j,1), i)
            grid % cells_n(2,-NbC) = grid % cells_n(f4n(j,2), i)
            grid % cells_n(3,-NbC) = grid % cells_n(f4n(j,3), i)

            grid % faces_n_nodes(NS) = 3
            grid % faces_n(1,NS) = grid % cells_n(1,-NbC)
            grid % faces_n(2,NS) = grid % cells_n(2,-NbC)
            grid % faces_n(3,NS) = grid % cells_n(3,-NbC)
          end if
        end if

        ! Pyramides:
        if(grid % cells_n_nodes(i) == 5) then
          if(j == 1) then    ! face (1)
            grid % cells_n_nodes(-NbC) = 4
            grid % cells_n(1,-NbC) = grid % cells_n(f5n(j,1), i)
            grid % cells_n(2,-NbC) = grid % cells_n(f5n(j,2), i)
            grid % cells_n(3,-NbC) = grid % cells_n(f5n(j,3), i)
            grid % cells_n(4,-NbC) = grid % cells_n(f5n(j,4), i)
 
            grid % faces_n_nodes(NS) = 4
            grid % faces_n(1,NS) = grid % cells_n(1,-NbC)
            grid % faces_n(2,NS) = grid % cells_n(2,-NbC)
            grid % faces_n(3,NS) = grid % cells_n(3,-NbC)
            grid % faces_n(4,NS) = grid % cells_n(4,-NbC)
          else if(j <= 5) then
            grid % cells_n_nodes(-NbC) = 3
            grid % cells_n(1,-NbC) = grid % cells_n(f5n(j,1), i)
            grid % cells_n(2,-NbC) = grid % cells_n(f5n(j,2), i)
            grid % cells_n(3,-NbC) = grid % cells_n(f5n(j,3), i)
 
            grid % faces_n_nodes(NS) = 3
            grid % faces_n(1,NS) = grid % cells_n(1,-NbC)
            grid % faces_n(2,NS) = grid % cells_n(2,-NbC)
            grid % faces_n(3,NS) = grid % cells_n(3,-NbC)
          end if
        end if

      end if
    end do 
  end do
  write(*,*) 'Boundary cells: ', NBc

  end subroutine Grid_Topology
