!==============================================================================!
  subroutine Convective_Outflow(grid)
!------------------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use Grid_Mod
  use Work_Mod, only: t_x => r_cell_01,  &
                      t_y => r_cell_02,  &
                      t_z => r_cell_03           
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
!==============================================================================!

  call Compute_Fluxes(grid)

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! On the boundary perform the extrapolation
    if(c2  < 0) then
      if( (TypeBC(c2) == CONVECT) ) then
        U % n(c2) = U % n(c2) - ( Ubulk(material(c1)) * u % x(c1)         & 
                                + Vbulk(material(c1)) * u % y(c1)         &
                                + Wbulk(material(c1)) * u % z(c1) ) * dt
        V % n(c2) = V % n(c2) - ( Ubulk(material(c1)) * v % x(c1)         & 
                                + Vbulk(material(c1)) * v % y(c1)         &
                                + Wbulk(material(c1)) * v % z(c1) ) * dt
        W % n(c2) = W % n(c2) - ( Ubulk(material(c1)) * w % x(c1)         & 
                                + Vbulk(material(c1)) * w % y(c1)         &
                                + Wbulk(material(c1)) * w % z(c1) ) * dt
      end if
    end if
  end do

  if(HOT==YES) then
    call GraPhi(grid, t % n, 1, t_x, .TRUE.)     ! dT/dx
    call GraPhi(grid, t % n, 2, t_y, .TRUE.)     ! dT/dy
    call GraPhi(grid, t % n, 3, t_z, .TRUE.)     ! dT/dz
    call GraCorNew(grid, T % n, t_x, t_y, t_z) ! needed ?
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2  < 0) then
        if( (TypeBC(c2) == CONVECT) ) then
          T % n(c2) = T % n(c2) - ( Ubulk(material(c1)) * t_x(c1)        & 
                                  + Vbulk(material(c1)) * t_y(c1)        &
                                  + Wbulk(material(c1)) * t_z(c1) ) * dt
        end if
      end if
    end do
  end if

  end subroutine
