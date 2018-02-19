!======================================================================!
  subroutine Time_And_Length_Scale(grid)
!----------------------------------------------------------------------!
!  Purpose:                                                            !
!  Calculates time scale and leght scale in manner to avoid singularity!
!  in eps equation.                                                     !
!                                                                      !
!  Authors: Muhamed Hadziabdic 
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use Flow_Mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!----------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!-------------------------------[Locals]-------------------------------!
  integer :: c 
  real    :: T1, T2, L1, L2, L3, T3
!======================================================================!

  if(turbulence_model == K_EPS_V2) then
    do c = 1, grid % n_cells 
      T1 = kin%n(c)/(eps%n(c) )
      T2 = Ct*(abs(viscosity/(eps%n(c) )))**0.5
      T3 = 0.6*kin % n(c) / ( CmuD * v2 % n(c) * Shear(c) *&
           sqrt(3.0))
      L1 = kin%n(c)**1.5/(eps%n(c))
      L2 = Cni*(viscosity**3/(eps%n(c)))**0.25
      L3 = kin % n(c)**1.5 / ( CmuD * v2 % n(c) *&
           Shear(c) * sqrt(3.0))
      Tsc(c) = max(min(T1,T3),T2)
      Lsc(c) = Cl*max(min(L1,L3),L2)
    end do
  else if(turbulence_model == K_EPS_ZETA_F .or.  &
          turbulence_model == HYBRID_K_EPS_ZETA_F) then
    if(ROUGH == YES) then
      do c = 1, grid % n_cells 
        T1 = kin%n(c)/(eps%n(c))
        T2 = Ct*sqrt(viscosity/eps%n(c))

        L1 = kin%n(c)**1.5/eps%n(c)
        L2 = Cni*(viscosity**3/eps%n(c))**0.25

        Tsc(c) = max(T1,T2)
        Lsc(c) = Cl*max(L1,L2)
      end do
    else
      do c = 1, grid % n_cells
        T1 = kin%n(c)/(eps%n(c) + tiny)
        T2 = Ct*sqrt(viscosity/eps%n(c))
        T3 = 0.6/(sqrt(3.0)*CmuD * v2 % n(c) * Shear(c))

        L1 = kin%n(c)**1.5/eps%n(c)
        L2 = Cni*(viscosity**3/eps%n(c))**0.25
        L3 = sqrt(kin % n(c)/3.0)/(CmuD * v2 % n(c) * Shear(c)) 

        Tsc(c) = max(min(T1,T3),T2)
        Lsc(c) = Cl*max(min(L1,L3),L2)
      end do
    end if 
  else if(turbulence_model == REYNOLDS_STRESS_MODEL) then
    do c = 1, grid % n_cells
      kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)),1.0e-12)
      T1 = kin%n(c)/(eps%n(c))
      T2 = Ct*sqrt(viscosity/abs(eps%n(c)))

      L1 = kin%n(c)**1.5/eps%n(c)
      L2 = Cni*(viscosity**3/abs(eps%n(c)))**0.25

      Tsc(c) = max(T1,T2)
      Lsc(c) = Cl*max(L1,L2)
    end do
  end if

  end subroutine
