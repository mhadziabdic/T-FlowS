!==============================================================================!
  subroutine Convective_Outflow(grid, dt)
!------------------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod
  use all_mod
  use Flow_Mod
  use rans_mod
  use Grid_Mod
  use Bulk_Mod
  use Control_Mod
  use Work_Mod, only: t_x => r_cell_01,  &
                      t_y => r_cell_02,  &
                      t_z => r_cell_03           
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: dt
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
!==============================================================================!

  call Bulk_Mod_Compute_Fluxes(grid, bulk, flux)

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! On the boundary perform the extrapolation
    if(c2  < 0) then
      if( (Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT) ) then
        U % n(c2) = U % n(c2) - ( bulk(material(c1)) % u * u % x(c1)         & 
                                + bulk(material(c1)) % v * u % y(c1)         &
                                + bulk(material(c1)) % w * u % z(c1) ) * dt
        V % n(c2) = V % n(c2) - ( bulk(material(c1)) % u * v % x(c1)         & 
                                + bulk(material(c1)) % v * v % y(c1)         &
                                + bulk(material(c1)) % w * v % z(c1) ) * dt
        W % n(c2) = W % n(c2) - ( bulk(material(c1)) % u * w % x(c1)         & 
                                + bulk(material(c1)) % v * w % y(c1)         &
                                + bulk(material(c1)) % w * w % z(c1) ) * dt
      end if
    end if
  end do

  if(heat_transfer == YES) then
    call GraPhi(grid, t % n, 1, t_x, .true.)     ! dT/dx
    call GraPhi(grid, t % n, 2, t_y, .true.)     ! dT/dy
    call GraPhi(grid, t % n, 3, t_z, .true.)     ! dT/dz
    call GraCorNew(grid, T % n, t_x, t_y, t_z) ! needed ?
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2  < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT) ) then
          T % n(c2) = T % n(c2) - ( bulk(material(c1)) % u * t_x(c1)        & 
                                  + bulk(material(c1)) % v * t_y(c1)        &
                                  + bulk(material(c1)) % w * t_z(c1) ) * dt
        end if
      end if
    end do
  end if

  end subroutine
