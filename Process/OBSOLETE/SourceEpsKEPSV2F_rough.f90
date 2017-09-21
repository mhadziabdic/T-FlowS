!======================================================================!
  SUBROUTINE SourceEpsKEPSV2F_rough()
!-----------------------------------------------------------------------!
! Purpose:                                                              !
! Calculates source terms in equation of dissipation of turbulent energy! 
!  and imposes  a boundary condition                                    ! 
! model : k-eps-v2f                                                     !  
! Authors: Muhamed Hadziabdic and Bojan Niceno                          !
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
  REAL    :: Esor, Ce_11, Gblend, fp, fa, Rey, Ret
!--------------------------------[CVS]---------------------------------!
!  $Id: SourceEpsKEPSV2F_rough.f90,v 1.1 2017/08/31 22:35:25 mhadziabdic Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Process/SourceEpsKEPSV2F_rough.f90,v $  
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

    do c = 1,NC 
      Esor = volume(c)/(Tsc(c)+tiny)
      Ce_11 = Ce1*(1.0 + alpha*(1.0/(v_2%n(c)+tiny) ))    
      b(c) = b(c) + Ce_11*Pk(c)*Esor
 
!----- Fill in a diagonal of coefficient matrix 
      Aval(Adia(c)) =  Aval(Adia(c)) + Ce2*Esor*DENc(material(c))
    end do                   

!----  Imposing a boundary condition on wall for Eps 

    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
        if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then

          Eps%n(c1) = Cmu75 * Kin%n(c1)**1.5 / (WallDs(c1) * kappa)

!----- Adjusting a coefficient in order to get a fixed
!----- value of Eps in near wall calls

          do j=Acol(c1), Acol(c1+1) -1
            Aval(j) = 0.0
          end do
          Aval(Adia(c1)) = 1.0 
          b(c1)       = Eps % n(c1)
        end if
      end if
    end do  

  RETURN
  END SUBROUTINE SourceEpsKEPSV2F_rough  
