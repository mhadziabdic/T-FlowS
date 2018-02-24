!==============================================================================!
  subroutine Calculate_Vis_T_Spalart_Allmaras(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Flow_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c
  real              :: x_rat, f_v1, lf
!==============================================================================!

  if(turbulence_model == DES_SPALART) then
    do c = 1, grid % n_cells
      x_rat    = VIS % n(c)/viscosity
      f_v1     = x_rat**3/(x_rat**3 + Cvis1**3)
      vis_t(c) = density * f_v1 * VIS % n(c)
    end do
  end if

  if(turbulence_model == SPALART_ALLMARAS) then
    do c = 1, grid % n_cells
      x_rat    = VIS % n(c)/viscosity
      f_v1     = x_rat**3/(x_rat**3 + Cvis1**3)
      vis_t(c) = density * f_v1 * VIS % n(c)
    end do
  end if

  call Exchange(grid, vis_t)  

  end subroutine