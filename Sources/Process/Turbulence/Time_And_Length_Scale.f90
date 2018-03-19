!==============================================================================!
  subroutine Time_And_Length_Scale(grid)
!------------------------------------------------------------------------------!
!   Calculates time scale and leght scale in manner to avoid singularity       !
!   in eps equation.                                                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c 
  real    :: t1, t2, t3, l1, l2, l3
!==============================================================================!

  if(turbulence_model == K_EPS_ZETA_F .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F) then
    if(ROUGH == YES) then
      do c = 1, grid % n_cells 
        t1 = kin % n(c)/(eps % n(c))
        t2 = Ct*sqrt((viscosity/density)/eps % n(c))

        l1 = kin % n(c)**1.5/eps % n(c)
        l2 = Cni*((viscosity/density)**3/eps % n(c))**0.25

        Tsc(c) = max(t1,t2)
        Lsc(c) = Cl*max(l1,l2)
      end do
    else
      do c = 1, grid % n_cells
        t1 = kin % n(c)/(eps % n(c) + tiny)
        t2 = Ct*sqrt((viscosity/density)/eps % n(c))
        t3 = 0.6/(sqrt(3.0)*CmuD * zeta % n(c) * Shear(c))

        l1 = kin % n(c)**1.5/eps % n(c)
        l2 = Cni*((viscosity/density)**3/eps % n(c))**0.25
        l3 = sqrt(kin % n(c)/3.0)/(CmuD * zeta % n(c) * Shear(c)) 

        Tsc(c) = max(min(t1,t3),t2)
        Lsc(c) = Cl*max(min(l1,l3),l2)
      end do
    end if 
!  Purpose:                                                            !
  else if(turbulence_model == REYNOLDS_STRESS_MODEL) then
    do c = 1, grid % n_cells
      kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)),1.0e-12)
      t1 = kin % n(c)/(eps % n(c))
      t2 = Ct*sqrt(viscosity/abs(eps % n(c)))

      l1 = kin % n(c)**1.5/eps % n(c)
      l2 = Cni*(viscosity**3/abs(eps % n(c)))**0.25

      Tsc(c) = max(t1,t2)
      Lsc(c) = Cl*max(l1,l2)
    end do
  end if

  end subroutine
