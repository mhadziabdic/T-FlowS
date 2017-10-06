!======================================================================!
  subroutine SourceKinKEPSV2F()
!----------------------------------------------------------------------!

!------------------------------[Modules]-------------------------------!
!----------------------------------------------------------------------!
!   Computes the source terms in Kin transport equation.               !
!                                                                      !
!   model : k-eps-v2f and zeta-f                                       !
!                                                                      !
!   Authors: Muhamed Hadziabdic
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c, c1, c2, s
  real    :: Utan, UnorSq, Unor, UtotSq, dely, Stot
  real    :: lf, Gblend, Ustar, Ck, yPlus, Uplus, Pk_turb, Pk_vis
  real    :: EBF, EBF1
  real    :: ALPHA1, Lrans, Lsgs
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

  if(SIMULA == HYB_ZETA) then
    do c=1,NC
      lf = volume(c)**0.3333
      Lsgs  = 0.8*lf
      Lrans = 0.41*WallDs(c)
      ALPHA1 = max(1.0,Lrans/Lsgs)
!----- Production:
      b(c) = b(c) + VISt(c) * Shear(c) * Shear(c) * volume(c)
!----------------------------------------!
!----- Dissipation:
      if(ALPHA1 < 1.05) then
        Aval(Adia(c)) = Aval(Adia(c)) +                                    &
                        DENc(material(c))*Eps%n(c)/(Kin%n(c) + tiny)*volume(c)
      else
        Aval(Adia(c)) = Aval(Adia(c)) +                                    &
                        DENc(material(c))*min(ALPHA1**1.45*Eps%n(c),Kin%n(c)**1.5  &
                       /(lf*0.01))/(Kin%n(c) + tiny)*volume(c)
      end if
      Pk(c) =  VISt(c) * Shear(c) * Shear(c)
    end do
  else
    do c=1,NC
!----- Production:
      b(c) = b(c) + VISt(c) * Shear(c) * Shear(c) * volume(c)
 
!----------------------------------------!
!----- Dissipation:
      Aval(Adia(c)) = Aval(Adia(c)) +                                    &
                      DENc(material(c))*Eps%n(c)/(Kin%n(c)+tiny)*volume(c)
      Pk(c) =  VISt(c) * Shear(c) * Shear(c) 
      if (BUOY == YES) then ! XXXXX 6 Jun 2014
        buoyBeta(c) = 1.0
        Gbuoy(c) = -buoyBeta(c)*(grav_x*ut%n(c) + grav_y*vt%n(c) + grav_z*wt%n(c))   ! XXXXX 5 Jul 2014
        b(c) = b(c) + max(0.0,Gbuoy(c)*volume(c))
        Aval(Adia(c))= Aval(Adia(c))                                       &
                     + max(0.0,-Gbuoy(c)*volume(c)/(Kin%n(c) + tiny))
      end if
    end do
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
                 /sqrt(Sx(s)*Sx(s) + Sy(s)*Sy(s) + Sz(s)*Sz(s))
          UnorSq = Unor*Unor

          if( UtotSq   >  UnorSq) then
            Utan = sqrt(UtotSq - UnorSq)
          else
            Utan = TINY
          end if
    
          if(Ynd(c1) > 3.0) then 
            if(ROUGH==NO) then
              TauWall(c1) = DENc(material(c1))*kappa*Uf(c1)*Utan    &   
                           /(log(Elog*Ynd(c1)))    
              Pk(c1) = TauWall(c1)*Uf(c1)/(kappa*WallDs(c1))
            else if(ROUGH==YES) then
              TauWall(c1) = DENc(material(c1))*kappa*Uf(c1)*Utan    &   
                             /(log((WallDs(c1)+Zo)/Zo))    
              Pk(c1) = TauWall(c1)*Uf(c1)/(kappa*(WallDs(c1)+Zo))
              Kin%n(c2) = TauWall(c1)/0.09**0.5
            end if
            b(c1) = b(c1) + Pk(c1) * volume(c1)
            b(c1) = b(c1) - VISt(c1) * Shear(c1) * Shear(c1) * volume(c1)
          end if  
        end if  ! TypeBC(c2)==WALL or WALLFL
      end if    ! c2 < 0 
    end do

  RETURN
 
  end subroutine SourceKinKEPSV2F 
