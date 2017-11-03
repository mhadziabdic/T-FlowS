!======================================================================!
  subroutine SourceVisSpalart(dPHIdx,dPHIdy,dPHIdz)
!----------------------------------------------------------------------!
!   Computes the source terms in VIS transport equation.               !
!   model : spa-all                                                    !
!                                                                      !
!   Authors: Muhamed Hadziabdic and Bojan Niceno                       ! 
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c 
  real    :: Xrat, Fv1, Fv2, Fw, SS, DistV, ProdV, R, GG, Dif
  real    :: dist
  real    :: dPHIdx(-NbC:NC), dPHIdy(-NbC:NC), dPHIdz(-NbC:NC)
!======================================================================!

  if(SIMULA == SPA_ALL) then
    do c = 1, NC
!--------------------------------------------------!
!    Compute the production term                   !
!--------------------------------------------------!
      Xrat  = VIS % n(c)/VISc
      Fv1   = Xrat**3.0/(Xrat**3.0 + Cvis1**3.0)
      Fv2   = 1.0 - Xrat/(1.0 + Xrat*Fv1)
      SS    = Vort(c) + VIS % n(c)*Fv2/(kappa**2.0*WallDs(c)**2.0)
      ProdV = Cb1 * DENc(material(c)) * SS * VIS % n(c)
      b(c)  = b(c) + ProdV * volume(c)
!--------------------------------------------------!
!    Compute the destruction  term                 !
!--------------------------------------------------!
      R     = VIS % n(c)/(SS * kappa**2.0 * WallDs(c)**2.0)
      GG    = R + Cw2*(R**6.0 - R)
      Fw    = GG*((1.0 + Cw3**6.0)/(GG**6.0 + Cw3**6.0))**(1.0/6.0)
      DistV = Cw1* DENc(material(c)) * Fw * (VIS % n(c)/WallDs(c)**2.0)
      A % val(A % dia(c)) = A % val(A % dia(c)) + DistV * volume(c)
 
!--------------------------------------------------!
!    Compute the first-order diffusion term        !
!--------------------------------------------------!
      Dif   = Cb2 * DENc(material(c)) * (dPHIdx(c)**2.0 + dPHIdy(c)**2.0           &
              + dPHIdz(c)**2.0)/SIGMAv
      b(c)  = b(c) + Dif * volume(c)
    end do
  else if(SIMULA == DES_SPA) then
    do c = 1, NC

      dist = min(WallDs(c),0.65*delta(c))
!--->>>>      dist = WallDs(c)
!--------------------------------------------------!
!    Compute the production term                   !
!--------------------------------------------------!
      Xrat  = VIS % n(c)/VISc
      Fv1   = Xrat**3.0/(Xrat**3.0 + Cvis1**3.0)
      Fv2   = 1.0 - Xrat/(1.0 + Xrat*Fv1)
      SS    = Vort(c) + VIS % n(c)*Fv2/(kappa**2.0*dist**2.0)
      ProdV = Cb1 * DENc(material(c)) * SS * VIS % n(c)
      b(c)  = b(c) + ProdV * volume(c)
      
!--------------------------------------------------!
!    Compute the destruction  term                 !
!--------------------------------------------------!
      R     = VIS % n(c)/(SS * kappa**2.0 * dist**2.0)
      GG    = R + Cw2*(R**6.0 - R)
      Fw    = GG*((1.0 + Cw3**6.0)/(GG**6.0 + Cw3**6.0))**(1.0/6.0)
      DistV = Cw1* DENc(material(c)) * Fw * (VIS % n(c)/dist**2.0)
      A % val(A % dia(c)) = A % val(A % dia(c)) + DistV * volume(c)

!--------------------------------------------------!
!    Compute the first-order diffusion term        !
!--------------------------------------------------!
      Dif   = Cb2 * DENc(material(c)) * (PHIx(c) + PHIy(c) + PHIz(c))**2.0/SIGMAv
      b(c)  = b(c) + Dif * volume(c)
    end do 
  end if

  RETURN
  end subroutine SourceVisSpalart
