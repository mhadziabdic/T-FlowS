!======================================================================!
  SUBROUTINE SourceKinKEPSV2F_rough()
!----------------------------------------------------------------------!

!------------------------------[Modules]-------------------------------!
!----------------------------------------------------------------------!
!   Computes the source terms in Kin transport equation,               !
!   wall shear stress (wall function approuch)                         !
!                                                                      !
!   model : k-eps-v2f                                                  !
!                                                                      !
!   Authors: Muhamed Hadziabdic and Bojan Niceno                       !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: c, c1, c2, s
  REAL    :: Utan, UnorSq, Unor, UtotSq, dely, Stot
  REAL    :: lf, Gblend, Ustar, Ck, yPlus, Uplus
!--------------------------------[CVS]---------------------------------!
!  $Id: SourceKinKEPSV2F_rough.f90,v 1.1 2017/08/31 22:35:25 mhadziabdic Exp $  
!  $UserKDSource: /home/muhamed/.CVSROOT/T-Rex/Pro/KSource.f90,v $  
!======================================================================! 

!===================!
!                   !
!     Pk terms      !
!     in domain     !
!                   !
!===================!

!----------------------------------------------------------------------!
!                                                                      !
!   Pk  is a production of turbulent kinematic energy.                 !
!                                                                      !
!   The form of Pk which is solving in this subroutine is :            !
!                                                                      !
!                                                                      !
!    Pk = 2*VISk{ [ (dU/dx)**2 + (dV/dy)**2 + (dW/dz)**2 ] +           ! 
!        0.5( dV/dx + dU/dy )**2 + 0.5( dW/dy + dV/dz )**2 +           !
!        0.5(dU/dz + dW/dx )**2 }                                      !
!                                                                      !
!   Dimension of the Pk is [ kg m /s**3 ]                              !
!   In kinetic energy equation exist two source terms which have form: !
!                                                                      !
!    /             /                                                   !
!   |             |                                                    !
!   |Pk dV  -     |DENc Eps dV                                         !
!   |             |                                                    !
!  /             /                                                     !
!                                                                      !
!----------------------------------------------------------------------!

!  Zo = 0.000025

  do c=1,NC
!----- Production:
    b(c) = b(c) + VISt(c) * Shear(c) * Shear(c) * volume(c)
 
!----------------------------------------!
!----- Dissipation:
    Aval(Adia(c)) = Aval(Adia(c)) +                                    &
                    DENc(material(c))*Eps%n(c)/(Kin%n(c)+tiny)*volume(c)
    Pk(c) =  VISt(c) * Shear(c) * Shear(c) 
     
  end do

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
                 /sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))
          UnorSq = Unor*Unor

          if( UtotSq   >  UnorSq) then
            Utan = sqrt(UtotSq - UnorSq)
          else
            Utan = TINY
          end if
         
          Uf(c1)  = Cmu**0.25*Kin%n(c1)**0.5
          Ynd(c1) = (WallDs(c1)+Zo)*Uf(c1)/VISc 
          Gblend  = 0.01*Ynd(c1)**4.0/(1.0+5.0*Ynd(c1))

          if(Ynd(c1) > 3.0) then 
            TauWall(c1) = DENc(material(c1))*kappa*Uf(c1)*Utan    &
                         /(log((WallDs(c1)+Zo)/Zo))    
            Pk(c1) = TauWall(c1)*Uf(c1)/(kappa*(WallDs(c1)+Zo))

            b(c1) = b(c1) + Pk(c1) * volume(c1)
            b(c1) = b(c1) - VISt(c1) * Shear(c1) * Shear(c1) * volume(c1)
            Kin%n(c2) = TauWall(c1)/0.09**0.5
          end if  
        end if  ! TypeBC(c2)==WALL or WALLFL
      end if    ! c2 < 0
    end do

  RETURN
 
  END SUBROUTINE SourceKinKEPSV2F_rough 
