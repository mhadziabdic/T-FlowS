!======================================================================!
  subroutine Time_And_Length_Scale(grid)
!----------------------------------------------------------------------!
!  Purpose:                                                            !
!  Calculates time scale and leght scale in manner to avoid singularity!
!  in Eps equation.                                                     !
!                                                                      !
!  Authors: Muhamed Hadziabdic 
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use Grid_Mod
  use Constants_Pro_Mod
!----------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!-------------------------------[Locals]-------------------------------!
  integer c 
  real T1, T2, L1, L2, L3, T3
!======================================================================!

  if(SIMULA==K_EPS_VV) then
    do c = 1, grid % n_cells 
      T1 = Kin%n(c)/(Eps%n(c) )
      T2 = Ct*(abs(VISc/(Eps%n(c) )))**0.5
      T3 = 0.6*Kin % n(c) / ( CmuD * v_2 % n(c) * Shear(c) *&
           sqrt(3.0))
      L1 = Kin%n(c)**1.5/(Eps%n(c))
      L2 = Cni*(VISc**3/(Eps%n(c)))**0.25
      L3 = Kin % n(c)**1.5 / ( CmuD * v_2 % n(c) *&
           Shear(c) * sqrt(3.0))
      Tsc(c) = max(min(T1,T3),T2)
      Lsc(c) = Cl*max(min(L1,L3),L2)
    end do
  else if(SIMULA == ZETA.or.SIMULA==HYB_ZETA) then
    if(ROUGH == YES) then
      do c = 1, grid % n_cells 
        T1 = Kin%n(c)/(Eps%n(c))
        T2 = Ct*sqrt(VISc/Eps%n(c))

        L1 = Kin%n(c)**1.5/Eps%n(c)
        L2 = Cni*(VISc**3/Eps%n(c))**0.25

        Tsc(c) = max(T1,T2)
        Lsc(c) = Cl*max(L1,L2)
      end do
    else
      do c = 1, grid % n_cells
        T1 = Kin%n(c)/(Eps%n(c) + tiny)
        T2 = Ct*sqrt(VISc/Eps%n(c))
        T3 = 0.6/(sqrt(3.0)*CmuD * v_2 % n(c) * Shear(c))

        L1 = Kin%n(c)**1.5/Eps%n(c)
        L2 = Cni*(VISc**3/Eps%n(c))**0.25
        L3 = sqrt(Kin % n(c)/3.0)/(CmuD * v_2 % n(c) * Shear(c)) 

        Tsc(c) = max(min(T1,T3),T2)
        Lsc(c) = Cl*max(min(L1,L3),L2)
      end do
    end if 
  else if(SIMULA == EBM) then
    do c = 1, grid % n_cells
      Kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)),1.0e-12)
      T1 = Kin%n(c)/(Eps%n(c))
      T2 = Ct*sqrt(VISc/abs(Eps%n(c)))

      L1 = Kin%n(c)**1.5/Eps%n(c)
      L2 = Cni*(VISc**3/abs(Eps%n(c)))**0.25

      Tsc(c) = max(T1,T2)
      Lsc(c) = Cl*max(L1,L2)
    end do
  end if

  end subroutine
