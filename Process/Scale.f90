!======================================================================!
  SUBROUTINE Scale()
!----------------------------------------------------------------------!
!  Purpose:                                                            !
!  Calculates time scale and leght scale in manner to avoid singularity!
!  in Eps equation.                                                     !
!                                                                      !
!  Authors: Muhamed Hadziabdic 
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER c 
  REAL T1, T2, L1, L2, L3, T3
!--------------------------------[CVS]---------------------------------!
!  $Id: Scale.f90,v 1.2 2017/08/31 22:09:09 mhadziabdic Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Process/Scale.f90,v $              
!======================================================================!

  if(SIMULA==K_EPS_VV) then
    do c = 1, NC 
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
      do c = 1, NC 
        T1 = Kin%n(c)/(Eps%n(c))
        T2 = Ct*sqrt(VISc/Eps%n(c))

        L1 = Kin%n(c)**1.5/Eps%n(c)
        L2 = Cni*(VISc**3/Eps%n(c))**0.25

        Tsc(c) = max(T1,T2)
        Lsc(c) = Cl*max(L1,L2)
      end do
    else
      do c = 1, NC
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
    do c = 1, NC
      Kin % n(c) = max(0.5*(uu % n(c) + vv % n(c) + ww % n(c)),1.0e-12)
      T1 = Kin%n(c)/(Eps%n(c))
      T2 = Ct*sqrt(VISc/abs(Eps%n(c)))

      L1 = Kin%n(c)**1.5/Eps%n(c)
      L2 = Cni*(VISc**3/abs(Eps%n(c)))**0.25

      Tsc(c) = max(T1,T2)
      Lsc(c) = Cl*max(L1,L2)
    end do
  end if

   RETURN
   END SUBROUTINE Scale
