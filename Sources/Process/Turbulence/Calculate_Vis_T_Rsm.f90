!==============================================================================!
  subroutine Calculate_Vis_T_Rsm(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RSM models ('EBM' and 'HJ').          !
!   If hybrid option is used turbulent diffusivity is modeled by vis_t.        !
!   Otherwise, vis_t is used as false diffusion in order to increase           !
!   stability of computation.                                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
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
  real              :: cmu_mod                                        
!==============================================================================!

  call Calculate_shear_And_Vorticity(grid)

  if(turbulence_model .eq. HANJALIC_JAKIRLIC) then
    do c = 1, grid % n_cells
      kin % n(c) = 0.5*max(uu % n(c)+vv % n(c)+ww % n(c), 1.0e-12)

      cmu_mod = max(-(  uu % n(c) * u % x(c)               &
                      + vv % n(c) * v % y(c)               &
                      + ww % n(c) * w % z(c)               &
                      + uv % n(c) * (v % x(c) + u % y(c))  &
                      + uw % n(c) * (u % z(c) + w % x(c))  &
                      + vw % n(c) * (v % z(c) + w % y(c))) &
               / max(kin % n(c)**2 / ( eps_tot(c) + TINY) * shear(c)**2, 1.0e-12), 0.0)

      cmu_mod = min(0.12, cmu_mod) 
      vis_t(c) = cmu_mod * density * kin % n(c)**2 / ( eps_tot(c) + TINY)
    end do 
  else if(turbulence_model .eq. REYNOLDS_STRESS_MODEL) then
    do c=1, grid % n_cells
      kin % n(c) = 0.5*max(uu % n(c)+vv % n(c)+ww % n(c), 1.0e-12)

      cmu_mod = max(-(  uu % n(c) * u % x(c)  &
                      + vv % n(c) * v % y(c)  &
                      + ww % n(c) * w % z(c)  &
                      + uv % n(c) * (v % x(c) + u % y(c))  &
                      + uw % n(c) * (u % z(c) + w % x(c))  &
                      + vw % n(c) * (v % z(c) + w % y(c))) &
               / max(kin % n(c)**2 / eps % n(c) * shear(c)**2, 1.0e-12), 0.0)

      cmu_mod = min(0.12,cmu_mod)
      vis_t(c) = cmu_mod*density * kin%n(c) * kin%n(c) / eps % n(c)
    end do
  end if

  call Comm_Mod_Exchange(grid, vis_t)  

  end subroutine
