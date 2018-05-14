!==============================================================================!
  subroutine Time_And_Length_Scale(grid)
!------------------------------------------------------------------------------!
!   Calculates time scale and leght scale in manner to avoid singularity       !
!   in eps equation.                                                           !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Flow_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
  use Work_Mod, only: t1 => r_cell_01,  &  ! [s]
                      t2 => r_cell_02,  &  ! [s]
                      t3 => r_cell_03,  &  ! [s]
                      l1 => r_cell_04,  &  ! [m]
                      l2 => r_cell_05,  &  ! [m]
                      l3 => r_cell_06      ! [m]
!------------------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!----------------------------------[Locals]------------------------------------!
  real    :: kin_visc   ! kinematic viscosity [m^2/s]
!==============================================================================!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!------------------------------------------------------------------------------!
  kin_visc = viscosity/density

  t1(1:) = kin % n(1:)/(eps % n(1:) + TINY)
  t2(1:) = c_t*sqrt(kin_visc/(eps % n(1:) + TINY))

  l1(1:) = kin % n(1:)**1.5/(eps % n(1:) + TINY)
  l2(1:) = Cni*(kin_visc**3/(eps % n(1:) + TINY))**0.25
  

  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_K_EPS_ZETA_F) then
    if(ROUGH .eq. YES) then
      t_scale(1:) =     max(t1(1:),t2(1:))
      l_scale(1:) = c_l*max(l1(1:),l2(1:))
    else
      t3(1:) = 0.6/(sqrt(3.0)*c_mu_d * zeta % n(1:) * shear(1:) + TINY)
      l3(1:) = sqrt(kin % n(1:)/3.0)/(c_mu_d * zeta % n(1:) * shear(1:) + TINY)
      t_scale(1:) =     max(min(t1(1:),t3(1:)),t2(1:))
      l_scale(1:) = c_l*max(min(l1(1:),l3(1:)),l2(1:))
    end if
  else if(turbulence_model .eq. REYNOLDS_STRESS_MODEL) then
    kin % n(1:) = max(0.5*(uu % n(1:) + vv % n(1:) + ww % n(1:)), TINY)
    t_scale(1:) =     max(t1(1:),t2(1:))
    l_scale(1:) = c_l*max(l1(1:),l2(1:))
  end if

  end subroutine
