!======================================================================!
  SUBROUTINE CalcVISt_KepsV2F_rough() 
!----------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                  !
!                                                                      !
!   In the domain:                                                     !
!   ~~~~~~~~~~~~~~                                                     !
!   For k-eps model :                                                  !
!                       2                                              !
!   VISt = Cmu * rho * K  * Eps                                        ! 
!                                                                      !
!   On the boundary (wall viscosity):                                  !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                   !
!            +          kappa                                          !
!   VIStw = y  * VISt ----------                                       ! 
!                     E * ln(y+)                                       !
!                                                                      !
!    For k-eps-v2f model :                                             !
!                                                                      !
!    VISt = CmuD * rho * Tsc  * vv                                     !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!------------------------------[Calling]-------------------------------!
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, c1, c2, s
  REAL    :: Utan, UnorSq, Unor, UtotSq 
  REAL    :: lf, Gblend, Ustar, Ck, yPlus, Uplus, EBF
  REAL    :: Cmu1, beta, Prmol, Prturb

!--------------------------------[CVS]---------------------------------!
!  $Id: CalcVISt_KEPSV2F_rough.f90,v 1.1 2017/08/31 22:35:25 mhadziabdic Exp $
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Process/CalcVISt_KEPSV2F_rough.f90,v $
!======================================================================!

!======================================================================!
!  K-EPS-V2F model
!======================================================================! 
  call Scale()

  if(SIMULA == K_EPS_VV) then
    do c = 1,NC
      VISt(c) = CmuD*v_2%n(c)*Tsc(c)
    end do
  else if(SIMULA == ZETA) then
    do c = 1,NC
      VISt(c) = min(CmuD*v_2%n(c)*Kin % n(c)*Tsc(c),0.09*Kin % n(c)*Tsc(c))
    end do
  else if(SIMULA == HYB_ZETA) then
    do c = 1,NC
      VISt(c) = CmuD*v_2%n(c)*Kin % n(c)*Tsc(c)
      VISt_eff(c) = max(VISt(c),VISt_sgs(c))
    end do
    call Exchng(VISt_eff)  
  end if

    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
    
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
!----- Compute tangential velocity component
          UtotSq = U % n(c1) * U % n(c1) &
                 + V % n(c1) * V % n(c1) &
                 + W % n(c1) * W % n(c1)
          Unor = ( U % n(c1) * Sx(s)     &
                 + V % n(c1) * Sy(s)     &
                 + W % n(c1) * Sz(s) )   &
                 /(Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))**0.5
          UnorSq = Unor*Unor

          if( UtotSq   >  UnorSq) then
            Utan = sqrt(UtotSq - UnorSq)
          else
            Utan = TINY
          end if


          Uf(c1)  = Cmu**0.25*Kin%n(c1)**0.5
          Ynd(c1) = (WallDs(c1)+Zo)*Uf(c1)/VISc 
          Gblend  = 0.01*Ynd(c1)**4.0/(1.0+5.0*Ynd(c1))

          Yplus = Ynd(c1)
          Uplus = log((WallDs(c1)+Zo)/Zo)/(kappa + TINY) + TINY
          EBF = 0.01*Yplus**4.0/(1.0 + Yplus) + TINY
          VISwall(c1) = Ynd(c1)*VISc/(Yplus*exp(-1.0*EBF)+Uplus*exp(-1.0/EBF) + TINY)

          if(HOT==YES) then
            Prturb = 1.0 / ( 0.5882 + 0.228*VISt(c1)/VISc   &
                   - 0.0441 * (VISt(c1)/VISc)**2.0  &
                   * (1.0 - exp(-5.165*VISc/(VISt(c1)+tiny))) )
            Prmol = VISc * CAPc(material(c1)) / CONc(material(c1))
            beta = 9.24 * ((Prmol/Prturb)**0.75 - 1.0) * (1.0 + 0.28 * exp(-0.007*Prmol/Prturb))
            EBF = 0.01 * (Prmol*Yplus)**4.0 / (1.0 + 5.0 * Prmol**3.0 * Yplus) + TINY
            CONwall(c1) = Yplus*VISc*CAPc(material(c1))/(Yplus*Prmol*exp(-1.0 * EBF) &
                        + (Uplus + beta)*Prturb*exp(-1.0/EBF) + TINY)
          end if
        end if  ! TypeBC(c2)==WALL or WALLFL
      end if    ! c2 < 0
    end do
  RETURN
 
  call Exchng(VISt)  
  call Exchng(VISwall)  
  RETURN

  END SUBROUTINE CalcVISt_KepsV2F_rough  
