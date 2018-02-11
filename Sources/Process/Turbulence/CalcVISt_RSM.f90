!==============================================================================!
  subroutine Calcvist_RSM(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RSM models (EBM and HJ).              !
!   If hybrid option is used turbulent diffusivity is modeled by vis_t.         !
!   Otherwise, vis_t is used as false diffusion in order to increase            !
!   stability of computation.                                                  !
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
  integer :: c
  real    :: Cmu_mod                                        
!==============================================================================!

  call Compute_Shear_And_Vorticity(grid)

  if(SIMULA == HJ) then
    do c=1, grid % n_cells
      kin % n(c) = 0.5*max(uu%n(c)+vv%n(c)+ww%n(c),1.0e-12)

      Cmu_mod = max(-(  uu % n(c) * u % x(c)  &
                      + vv % n(c) * v % y(c)  &
                      + ww % n(c) * w % z(c)  &
                      + uv % n(c) * (v % x(c) + u % y(c))  &
                      + uw % n(c) * (u % z(c) + w % x(c))  &
                      + vw % n(c) * (v % z(c) + w % y(c))) &
               / max(kin % n(c)**2 / eps_tot(c) * Shear(c)**2, 1.0e-12), 0.0)

      Cmu_mod = min(0.12,Cmu_mod) 
      vis_t(c) = Cmu_mod * DENc(material(c)) * kin % n(c)**2 / eps_tot(c)
    end do 
  else if(SIMULA==EBM) then
    do c=1, grid % n_cells
      kin % n(c) = 0.5*max(uu%n(c)+vv%n(c)+ww%n(c),1.0e-12)

      Cmu_mod = max(-(  uu % n(c) * u % x(c)  &
                      + vv % n(c) * v % y(c)  &
                      + ww % n(c) * w % z(c)  &
                      + uv % n(c) * (v % x(c) + u % y(c))  &
                      + uw % n(c) * (u % z(c) + w % x(c))  &
                      + vw % n(c) * (v % z(c) + w % y(c))) &
               / max(kin % n(c)**2 / eps % n(c) * Shear(c)**2, 1.0e-12), 0.0)

      Cmu_mod = min(0.12,Cmu_mod)
      vis_t(c) = Cmu_mod*DENc(material(c)) * kin%n(c) * kin%n(c) / eps % n(c)
    end do
  end if

  call Exchange(grid, vis_t)  

  end subroutine
