!======================================================================!
  subroutine TanConv(Var, Cell, Convection)
!----------------------------------------------------------------------!
!   Calculates tangential momentum convective term for AnalyticWF      !
!----------------------------------------------------------------------!
!     
!  The convective term is calculated as:
!
!                                 ____
!                       dPhi      \     PhiTanS MFluxS
!    PhiConv =  rho Ui ------  =   \   ----------------
!                        dxi       /       volume(c)
!                                 /___
!
!  where:
!   
!    Phi = Utan ; for var = 1
!    Phi = T    ; for var = 2
!
!
!-----------------------------[Modules]--------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!----------------------------[Arguments]-------------------------------!
  integer :: Var, Cell
  real    :: Convection
!------------------------------[Locals]--------------------------------!
  integer :: j, s, c1, c2, c
  real    :: DENs, Us, Vs, Ws, MFluxS, PhiTanS, Stot
!--------------------------------[CVS]---------------------------------!
!  $Id: TanConv.f90,v 1.1 2017/09/06 06:17:30 mhadziabdic Exp $
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Process/TanConv.f90,v $
!======================================================================!


!-----Initialize
  Convection = 0.0
  PhiTanS = 0.0
  MFluxS = 0.0
  Us = 0.0
  Vs = 0.0
  Ws = 0.0

!-----Browse through cell faces
  do j = 1,8  !8 defined in UnkAlloc.f90
    if(CellFace(Cell,j) > 0) then !others are -1
      s  = CellFace(Cell,j)
      c1 = SideC(1,s)
      c2 = SideC(2,s)

!-----Interpolate density
      if( StateMat(material(c1))==FLUID .and.      &
          StateMat(material(c2))==FLUID) then
        DENs =      f(s)  * DENc(material(c1))     &
             + (1.0-f(s)) * DENc(material(c2))
      else if( StateMat(material(c1))==FLUID .and. &
               StateMat(material(c2))==SOLID) then
        DENs = DENc(material(c1))
      else if( StateMat(material(c1))==SOLID .and. &
               StateMat(material(c2))==FLUID) then
        DENs = DENc(material(c2))
      else
        DENs =      f(s)  * DENc(material(c1))     &
             + (1.0-f(s)) * DENc(material(c2))
      end if

!-----Interpolate velocities
      if( c2 > 0  .or.  (c2 < 0 .and. TypeBC(c2)==BUFFER) )then
        Us = f(s) * U % n(c1) + (1.0 - f(s)) * U % n(c2)
        Vs = f(s) * V % n(c1) + (1.0 - f(s)) * V % n(c2)
        Ws = f(s) * W % n(c1) + (1.0 - f(s)) * W % n(c2)
      else          !if c2 < 0
        if( TypeBC(c2) == INFLOW  .or. &
            TypeBC(c2) == OUTFLOW .or. &
            TypeBC(c2) == CONVECT ) then
          Us = U % n(c2)
          Vs = V % n(c2)
          Ws = W % n(c2)
!        else if(TypeBC(c2)==SYMMETRY) then !iskljucio simetriju
!          Us = U % n(c1)
!          Vs = V % n(c1)
!          Ws = W % n(c1)
        else           !it's WALL or WALLFL
          Us = 0.0
          Vs = 0.0
          Ws = 0.0
        end if
      end if

!-----Flux through cell face
      MFluxS = DENs * ( Us*Sx(s) + Vs*Sy(s) + Ws*Sz(s) )

!-----Calculate tangential velocity on the face
      if(Var==1) then
        Stot = sqrt( Sx(WallFace(Cell))*Sx(WallFace(Cell))  &
                   + Sy(WallFace(Cell))*Sy(WallFace(Cell))  &
                   + Sz(WallFace(Cell))*Sz(WallFace(Cell))  )
        if( c2 > 0  .or.  (c2 < 0 .and. TypeBC(c2)==BUFFER) ) then
          PhiTanS =         f(s) * ( U%n(c1)*sqrt(1.0-(Sx(WallFace(Cell))/Stot)**2.0)   &
                                   + V%n(c1)*sqrt(1.0-(Sy(WallFace(Cell))/Stot)**2.0)   &
                                   + W%n(c1)*sqrt(1.0-(Sz(WallFace(Cell))/Stot)**2.0) ) &
                  + (1.0 - f(s)) * ( U%n(c2)*sqrt(1.0-(Sx(WallFace(Cell))/Stot)**2.0)   &
                                   + V%n(c2)*sqrt(1.0-(Sy(WallFace(Cell))/Stot)**2.0)   &
                                   + W%n(c2)*sqrt(1.0-(Sz(WallFace(Cell))/Stot)**2.0) )
!  PhiTanS = f(s) * U%n(c1) + (1.0 - f(s)) * U%n(c2)
!  PhiTanS = (f(s)*U%n(c1) + (1.0-f(s))*U%n(c2)) * xsp(s)/sqrt(xsp(s)**2 + ysp(s)**2 + TINY) &
!          + (f(s)*V%n(c1) + (1.0-f(s))*V%n(c2)) * ysp(s)/sqrt(xsp(s)**2 + ysp(s)**2 + TINY)
        else             !if c2 < 0
          if( TypeBC(c2) == INFLOW  .or. &
              TypeBC(c2) == OUTFLOW .or. &
              TypeBC(c2) == CONVECT ) then
            PhiTanS = U%n(c2) * sqrt(1.0 - (Sx(WallFace(Cell))/Stot)**2.0)  &
                    + V%n(c2) * sqrt(1.0 - (Sy(WallFace(Cell))/Stot)**2.0)  &
                    + W%n(c2) * sqrt(1.0 - (Sz(WallFace(Cell))/Stot)**2.0)
!  PhiTanS = U%n(c2)
!  PhiTanS = U%n(c2) * xc(c)/sqrt(xc(c)**2 + yc(c)**2 + TINY) &
!          + V%n(c2) * yc(c)/sqrt(xc(c)**2 + yc(c)**2 + TINY)
          else if(TypeBC(c2) == SYMMETRY) then
            PhiTanS = U%n(c1) * sqrt(1.0 - (Sx(WallFace(Cell))/Stot)**2.0)  &
                    + V%n(c1) * sqrt(1.0 - (Sy(WallFace(Cell))/Stot)**2.0)  &
                    + W%n(c1) * sqrt(1.0 - (Sz(WallFace(Cell))/Stot)**2.0)
!  PhiTanS = U%n(c1)
!  PhiTanS = U%n(c1) * xc(c)/sqrt(xc(c)**2 + yc(c)**2 + TINY) &
!          + V%n(c1) * yc(c)/sqrt(xc(c)**2 + yc(c)**2 + TINY)
          end if
        end if

!-----Calculate temperature on the face
      else if(Var==2) then
        if(c2 > 0  .or.  (c2 < 0 .and. TypeBC(c2)==BUFFER) )then
          PhiTanS = f(s) * T % n(c1) + (1.0 - f(s)) * T % n(c2)
        else             !if c2 < 0
          if( TypeBC(c2) == INFLOW  .or. &
              TypeBC(c2) == OUTFLOW .or. &
              TypeBC(c2) == CONVECT ) then
            PhiTanS = T % n(c2)
          else if(TypeBC(c2) == SYMMETRY) then
            PhiTanS = T % n(c1)
          end if
        end if
      end if        !end if Var==1,2,3

!-----Calculate convection
      if(Cell==c1) Convection = Convection + MFluxS * PhiTanS / volume(Cell)
      if(Cell==c2) Convection = Convection - MFluxS * PhiTanS / volume(Cell)
!  Convection = Convection - MFluxS * PhiTanS / volume(Cell)

    end if  !end if CellFace > 0
  end do    !end do j=1,8


  RETURN 
  end subroutine TanConv
