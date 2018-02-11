!==============================================================================!
  subroutine Calcvist_SPA_ALL(grid, n) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Constants_Pro_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n 
  real    :: Xrat, Fv1, lf, Cs 
!==============================================================================!

  if(SIMULA == DES_SPA) then
    do c = 1, grid % n_cells
      Xrat     = VIS % n(c)/VISc
      Fv1      = Xrat**3/(Xrat**3 + Cvis1**3)
      vis_t(c) = DENc(material(c)) * Fv1 * VIS % n(c)
    end do
  end if

  if(SIMULA == SPA_ALL) then
    do c = 1, grid % n_cells
      Xrat     = VIS % n(c)/VISc
      Fv1      = Xrat**3/(Xrat**3 + Cvis1**3)
      vis_t(c) = DENc(material(c)) * Fv1 * VIS % n(c)
    end do
  end if

  call Exchange(grid, vis_t)  

  end subroutine
