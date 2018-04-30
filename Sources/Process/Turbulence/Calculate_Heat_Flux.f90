!==============================================================================!
  subroutine Calculate_Heat_Flux(grid)
!------------------------------------------------------------------------------!
!   Computes turbulent heat fluxes                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Grid_Mod
  use Grad_Mod
  use Flow_Mod 
  use rans_mod
  use Work_Mod, only: t_x => r_cell_01,  &
                      t_y => r_cell_02,  &
                      t_z => r_cell_03    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c 
  real    :: beta, pr, pr_t
!==============================================================================!

  call Grad_Mod_For_Phi(grid, t % n, 1, t_x, .true.)
  call Grad_Mod_For_Phi(grid, t % n, 2, t_y, .true.)
  call Grad_Mod_For_Phi(grid, t % n, 3, t_z, .true.)

!-------------------------------------------!
!    Compute the sources in the interior    !
!-------------------------------------------!
  pr   = 0.71  ! bad, hard coded
  beta = 1.0
  pr_t = 0.9

  if(turbulence_model == K_EPS               .or.  &
     turbulence_model == K_EPS_ZETA_F        .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F .or.  &
     turbulence_model == DES_SPALART) then

    do c = 1, grid % n_cells
      pr_t = 1.0/(   0.5882                                                &
                   + 0.228  * (vis_t(c)/(viscosity+1.0e-12))               &
                   - 0.0441 * (vis_t(c)/(viscosity+1.0e-12))**2.0          & 
                   * (1.0 - exp(-5.165*( viscosity/(vis_t(c)+1.0e-12) )))  &
                 )

      ut % n(c) = -vis_t(c) / pr_t * t_x(c)
      vt % n(c) = -vis_t(c) / pr_t * t_y(c)
      wt % n(c) = -vis_t(c) / pr_t * t_z(c)
    end do

  else if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
          turbulence_model == HANJALIC_JAKIRLIC) then
    do c = 1, grid % n_cells

      ! Pr_t computed but not used?
      pr_t = 1.0/(   0.5882                                                &
                   + 0.228  * (vis_t(c)/(viscosity+1.0e-12))               &
                   - 0.0441 * (vis_t(c)/(viscosity+1.0e-12))**2.0          & 
                   * (1.0 - exp(-5.165*( viscosity/(vis_t(c)+1.0e-12) )))  &
                 )

      ut % n(c) =  -0.22*Tsc(c) * (uu % n(c) * t_x(c) +  &
                                   uv % n(c) * t_y(c) +  &
                                   uw % n(c) * t_z(c))
      vt % n(c) =  -0.22*Tsc(c) * (uv % n(c) * t_x(c) +  &
                                   vv % n(c) * t_y(c) +  &
                                   vw % n(c) * t_z(c))
      wt % n(c) =  -0.22*Tsc(c) * (uw % n(c) * t_x(c) +  &
                                   vw % n(c) * t_y(c) +  &
                                   ww % n(c) * t_z(c))

    end do
  end if

  if(buoyancy == YES) then 
    ut % n(c) = min( 0.01 * t_ref, ut % n(c))
    ut % n(c) = max(-0.01 * t_ref, ut % n(c))
    vt % n(c) = min( 0.01 * t_ref, vt % n(c))
    vt % n(c) = max(-0.01 * t_ref, vt % n(c))
    wt % n(c) = min( 0.01 * t_ref, wt % n(c))
    wt % n(c) = max(-0.01 * t_ref, wt % n(c))
    Pbuoy(c) = -beta*(  grav_x * ut % n(c)  &
                      + grav_y * vt % n(c)  &
                      + grav_z * wt % n(c))
    Pbuoy(c) = max(Pbuoy(c),0.0)
  end if

  end subroutine
