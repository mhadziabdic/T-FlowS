!======================================================================!
  SUBROUTINE SourceEpsKeps
!----------------------------------------------------------------------!
!   Computes the source terms in the Eps transport equation,           !
!   wall shear stress (wall function approuch)                         !
!                                                                      !
!   model : k-eps                                                      !
!                                                                      !
!   Authors: Muhamed Hadziabdic and Bojan Niceno                       !
!----------------------------------------------------------------------! 
!                                                                      !
!                           2                                          !
!       Eps              Eps                                           !
!   Ce1 --- Gk - Ce2 rho ---                                           !
!        K                K                                            !
!                                                                      !
!   assigns epsilon from the wall function:                            !
!                                                                      !
!               3/4      3/2                                           !
!   Eps_w = Cmu^    * K^     / (Kappa/y)                               !
!                                                                      !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE pro_mod
  USE les_mod
  USE rans_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  INTEGER :: s, c, c1, c2, j, i
  REAL    :: Ret, Fmu, L1, L2, YAP, T1, yStar, Ce2star, nu_rng, Lf
!--------------------------------[CVS]---------------------------------!
! '$Id: SourceEpsKEps.f90,v 1.2 2017/08/31 22:32:06 mhadziabdic Exp $' 
! '$Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Process/SourceEpsKEps.f90,v $' 
!======================================================================!

  if(MODE == HRe) then
    do c=1,NC
!----- Positive contribution:
      b(c) = b(c) & 
           + Ce1 * Eps % n(c) / Kin % n(c) * Pk(c) * volume(c)
    
!----- Negative contribution:
      Aval(Adia(c)) = Aval(Adia(c)) &
                    + Ce2 * DENc(material(c)) * Eps % n(c) / Kin % n(c) * volume(c)
    end do 
!------------------------------------!
!   Cut-off the wall influence in order to
!   impose the boundary condition for EPS
!------------------------------------!
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)

      if(c2 < 0.and.TypeBC(c2) /= BUFFER ) then  
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then
!----- This will fix the value of Eps in the first cell

!          if(ROUGH==NO) then
            Eps % n(c1) = Cmu75 * (Kin % n(c1))**1.5 / (kappa*WallDs(c1))
!          else if(ROUGH==YES) then
!            Eps % n(c1) = Cmu75*(Kin%n(c1))**1.5/(kappa*(WallDs(c1)+Zo))
!          end if
          do j=Acol(c1), Acol(c1+1) -1
            Aval(j) = 0.0
          end do   
          Aval(Adia(c1)) = 1.0
          b(c1) = Eps % n(c1)
        end if  ! TypeBC(c2)==WALL or WALLFL
      end if    ! c2 < 0
    end do
  end if        ! end if mode = wf   


!=====================================================================
!   Jones-Launder model 
!=====================================================================
  if(MODE == LRe) then

   call CalcShear(U % n, V % n, W % n, Shear)
   call GraPhi(Shear,1,PHIx, .TRUE.)  ! dU/dx
   call GraPhi(Shear,2,PHIy, .TRUE.)  ! dW/dy
   call GraPhi(Shear,3,PHIz, .TRUE.)  ! dV/dz

    do c=1,NC
!----- Positive contribution:
      b(c) = b(c) & 
           + Ce1 * Eps % n(c) / Kin % n(c) * Pk(c) * volume(c)        & 
           + 2.0 * VISc * VISt(c) * &
           (PHIx(c)*PHIx(c)+PHIy(c)*PHIy(c)+PHIz(c)*PHIz(c)) * volume(c)
    
!----- Negative contribution:
      Ret = Kin % n(c)*Kin % n(c)/(VISc*Eps % n(c))
      Fmu = 1.0 - 0.3*exp(-(Ret*Ret))
      Aval(Adia(c)) = Aval(Adia(c))                                   &
      + Fmu * Ce2 * DENc(material(c)) * Eps % n(c) / Kin % n(c) * volume(c)        

!----- Yap correction
       L1 = Kin % n(c)**1.5/Eps % n(c)
       L2 = 2.55 * WallDs(c)
       YAP = 0.83 * Eps % n(c) * Eps % n(c)/Kin % n(c) * max((L1/L2 -1.0)*(L1/L2)**2.0,0.0)
       b(c) = b(c) + YAP * volume(c) 
    end do 
!------------------------------------!
!   Boundary condition fo EPS        !
!------------------------------------!
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)

      if(c2 < 0.and.TypeBC(c2) /= BUFFER ) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then

          Eps % n(c2) = 0.0

        end if  ! TypeBC(c2)==WALL or WALLFL
      end if    ! c2 < 0
    end do
  end if

  if(SIMULA == HYB_PITM) then
    do c=1,NC
!----- Positive contribution:
      b(c) = b(c) &
           + Ce1 * Eps % n(c) / Kin % n(c) * Pk(c) * volume(c)

      Lf = volume(c)**0.33333

!----- Negative contribution:
      Ret = Kin % n(c)*Kin % n(c)/(VISc*Eps % n(c))
      yStar = (VISc * Eps % n(c))**0.25 * WallDs(c) / VISc
      Fmu = (1.0 - exp(-yStar/3.1))**2.0*(1.0                             &
            - 0.25*exp(-(Ret/6.)*(Ret/6.)))

      Fmu = min(Fmu,1.0)

      Ce2 =  1.5 + 0.4/(1.0 + 2.4*(0.41*WallDs(c)/Lf)**0.66666)

      Aval(Adia(c)) = Aval(Adia(c)) &
      + (Ce1 + (Ce2 - Ce1) * Fmu ) * DENc(material(c)) * Eps % n(c) / Kin % n(c) * volume(c)

    end do

!------------------------------------!
!   Boundary condition fo EPS        !
!------------------------------------!
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)

      if(c2 < 0.and.TypeBC(c2) /= BUFFER ) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then

          Eps % n(c2) = 2.0 * VISc * Kin % n(c1)/(WallDs(c1)*WallDs(c1))

        end if  ! TypeBC(c2)==WALL or WALLFL
      end if    ! c2 < 0
    end do
  end if        ! end if mode = stan

  return
  END SUBROUTINE SourceEpsKeps   

