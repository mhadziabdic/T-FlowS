!======================================================================!
  SUBROUTINE SourceEpsKEPSV2F()
!-----------------------------------------------------------------------!
! Purpose:                                                              !
! Calculates source terms in equation of dissipation of turbulent energy! 
! and imposes  boundary condition                                       ! 
! model : k-eps-v2f and zeta-f                                          !  
! Authors: Muhamed Hadziabdic 
!-----------------------------------------------------------------------!
!------------------------------[Modules]--------------------------------!
  USE all_mod
  USE pro_mod
  USE rans_mod
  USE les_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, s, c1, c2,j	 
  REAL    :: Esor, Ce_11, Gblend, fp, fa, Rey, Ret, Ck, EBF
  REAL    :: Utan, UnorSq, Unor, UtotSq, dely, Stot, EpsWall, EpsHom
  REAL    :: BL_EPS, Pro, Pk_turb, Pk_vis, Yplus
!======================================================================!


!----------------------------------------------------------------------!
!    In dissipation of turbulent kinetic energy equation exist two     !
!    source terms which have form:                                     !
!                                                                      !
!    /                                                                 !
!   |                                                                  !
!   |((Cv_e1*PkE - Cv_11 DENc*Eps)/Tsc)*dV                             !
!   |                                                                  !
!  /                                                                   !
!                                                                      !
!  First, positive , source term is solved and added to source         !
!  coefficient b(c) on right hand side.                                !
!  Second, negative, source term is added to main diagonal left hand   !
!  side coefficient matrix in order to increase stability of solver    !
!  It is nessesary to calculate coefficient Cv_11 using Kin, Cv_e2, vi2!
!  and coefficient A1                                                  !
!----------------------------------------------------------------------!      

  call Scale()

  if(SIMULA == ZETA.or.SIMULA==HYB_ZETA) then
    do c = 1,NC 
      Esor = volume(c)/(Tsc(c)+tiny)
      Ce_11 = Ce1*(1.0 + alpha*(1.0/(v_2%n(c)+tiny) ))    
      b(c) = b(c) + Ce_11*Pk(c)*Esor
 
!----- Fill in a diagonal of coefficient matrix 
      Aval(Adia(c)) =  Aval(Adia(c)) + Ce2*Esor*DENc(material(c))
    end do                   
  end if
!----  Imposing a boundary condition on wall for Eps 

    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
       if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
          UtotSq = U % n(c1) * U % n(c1) &
                 + V % n(c1) * V % n(c1) &
                 + W % n(c1) * W % n(c1) 
          Unor = ( U % n(c1) * Sx(s)     &
                 + V % n(c1) * Sy(s)     &
                 + W % n(c1) * Sz(s) )   &
                 /sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))  
          UnorSq = Unor*Unor

          if( UtotSq   >  UnorSq) then
            Utan = sqrt(UtotSq - UnorSq)
          else
            Utan = TINY
          end if

          EpsWall = 2.0 * VISc * Kin%n(c1) / WallDs(c1)**2.0
          EpsHom = Cmu75 * Kin%n(c1)**1.5 / (WallDs(c1) * kappa)
          Uf(c1) = Cmu25 * Kin%n(c1)**0.5

          Pk_turb = Cmu75 * Kin%n(c1)**1.5 / (WallDs(c1) * kappa)
          Pk_vis = VISt(c1)*Shear(c1)*Shear(c1)                                   !standard
!--Kader
          Yplus = Cmu25 * sqrt(Kin%n(c1)) * WallDs(c1) / VISc + TINY     !standard
          EBF = 0.01 * Yplus**4.0 / (1.0 + 5.0*Yplus) + TINY           !original
          Pro = Pk_vis * exp(-1.0 * EBF) + Pk_turb * exp(-1.0 / EBF) 
          BL_EPS = min((Pk_turb * exp(-1.0 / EBF))/Pro,1.0)
          Ynd(c1) = Uf(c1)*WallDs(c1)/VISc 

          Ck = Cmu**0.25*Kin%n(c1)**0.5
          Ynd(c1) = Ck*WallDs(c1)/VISc 
          EBF  = 0.001*Ynd(c1)**4.0/(1.0+Ynd(c1))
          Eps%n(c1) = EpsWall * exp(-1.0 * EBF) + EpsHom * exp(-1.0 / EBF) 
          
          if(ROUGH == YES) then
            Eps%n(c1) = Cmu75 * Kin%n(c1)**1.5 / ((WallDs(c1)) * kappa)
          end if
!-----Adjusting a coefficient to fix Eps value in near wall calls
          do j=Acol(c1), Acol(c1+1)-1
            Aval(j) = 0.0
          end do
          b(c1) = Eps % n(c1)
          Aval(Adia(c1)) = 1.0
       end if
      end if
    end do  

  if(SIMULA == K_EPS_VV) then
    do c = 1,NC
      Esor = volume(c)/Tsc(c)
      Ce_11 = Ce1*(1.0 + alpha*(Kin%n(c)/(v_2%n(c)) + tiny)**0.5)
      b(c) = b(c) + Ce_11*Pk(c)*Esor

!----- Fill in a diagonal of coefficient matrix
      Aval(Adia(c)) =  Aval(Adia(c)) + Ce2*Esor*DENc(material(c))
    end do
  end if 

  RETURN
  END SUBROUTINE SourceEpsKEPSV2F 
