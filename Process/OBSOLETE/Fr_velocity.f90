!======================================================================!
  subroutine Fr_velocity()
!----------------------------------------------------------------------!
!  Purpose:                                                            !
!  Calculate Yplus in the near wall cells in order to perform swiching !
!  from  wall function approach to up to the wall approach             !
!                                                                      !
!  Authors: Muhamed Hadziabdic and Bojan Niceno                        !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer c1,c2, s 
  real :: UtotSq, Unor, UnorSq, Utan
!--------------------------------[CVS]---------------------------------!
!  $Id: Fr_velocity.f90,v 1.2 2017/08/31 21:53:38 mhadziabdic Exp $
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Process/Fr_velocity.f90,v $
!======================================================================!

  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)

    if(c2 < 0 .and. TypeBC(c2) /= BUFFER) then
      if(TypeBC(c2) == WALL .or. TypeBC(c2) == WALLFL) then
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

        Uf(c1)  = (CmuD*Kin%n(c1)*v_2%n(c1))**0.25  
        Ynd(c1) = Uf(c1)*WallDs(c1)/VISc
      end if
    end if
  end do 

  call Exchng(VISt)

  RETURN 
  end subroutine Fr_velocity
