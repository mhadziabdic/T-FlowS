!======================================================================!
  subroutine GraPhiCor(phi, phi_x, phi_y, phi_z)
!----------------------------------------------------------------------!
! Calculates gradient in the cells adjencent to material interface     !
! boundaries. Assumes that "tentative" gradients are just calculated   !
! and stored in "phi_x", "phi_y" and "phi_z" arrays.                   !
!                                                                      !
! It also assumes that the gradients "phi_x", "phi_y" and "phi_z"      !
! are fresh in buffers.                                                !
!                                                                      !
! This entire procedure is for two materials.                          !
!                                                                      !
! It is not desiged, and probably won't work in the following          !
! situations:                                                          !
!                                                                      !
!    +---+---+---+                                                     !
!    | F | F | F |                                                     !
!    +---+---+---+                                                     !
!    | F | S | S |  ->  The cell in the middle will probably not work  !
!    +---+---+---+                                                     !
!    | F | S | S |                                                     !
!    +---+---+---+                                                     !
!                                                                      !
! Further, it will probably not work on periodic boundaries.           !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use Solver_Mod, only: p1, p2  ! bad practice
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  real :: phi(-NbC:NC),    &
          phi_x(-NbC:NC),  &
          phi_y(-NbC:NC),  &
          phi_z(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  integer :: s, c, c1, c2
  real    :: Dphi1, Dphi2, Dxc1, Dyc1, Dzc1, Dxc2, Dyc2, Dzc2 
  real    :: f1, f2, phi_f
!======================================================================!

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!---- Take care of material interfaces          
    if( StateMat(material(c1))==FLUID .and. &  
        StateMat(material(c2))==SOLID       &  
        .or.                                &  
        StateMat(material(c1))==SOLID .and. &  
        StateMat(material(c2))==FLUID ) then   

      Dxc1 = xsp(s)-xc(c1)                     
      Dyc1 = ysp(s)-yc(c1)                     
      Dzc1 = zsp(s)-zc(c1)                     
      Dxc2 = xsp(s)-xc(c2)                     
      Dyc2 = ysp(s)-yc(c2)                     
      Dzc2 = zsp(s)-zc(c2)                     

!---- Missing parts of the gradient vector
      p1(c1) = CONc(material(c1)) *  &
           ( (G(1,c1)*Dxc1+G(4,c1)*Dyc1+G(5,c1)*Dzc1) * Sx(s) + &
             (G(4,c1)*Dxc1+G(2,c1)*Dyc1+G(6,c1)*Dzc1) * Sy(s) + & 
             (G(5,c1)*Dxc1+G(6,c1)*Dyc1+G(3,c1)*Dzc1) * Sz(s) )
      if(c2 > 0) then               
        p2(c2) = CONc(material(c2)) *  &
              ( (G(1,c2)*Dxc2+G(4,c2)*Dyc2+G(5,c2)*Dzc2) * Sx(s) + &
                (G(4,c2)*Dxc2+G(2,c2)*Dyc2+G(6,c2)*Dzc2) * Sy(s) + &
                (G(5,c2)*Dxc2+G(6,c2)*Dyc2+G(3,c2)*Dzc2) * Sz(s) )
      else if(TypeBC(c2) == BUFFER) then ! prepare to exchange
        p2(c1) = -p1(c1)
      end if
    end if    
  end do
 
  call Exchng(p2)

  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)

!---- Take care of material interfaces          
    if( StateMat(material(c1))==FLUID .and. &  
        StateMat(material(c2))==SOLID       &  
        .or.                                &  
        StateMat(material(c1))==SOLID .and. &  
        StateMat(material(c2))==FLUID ) then   

!---- Flux from cell 1 towards the material interface
      f1 = CONc(material(c1)) *  &
           (phi_x(c1)*Sx(s) + phi_y(c1)*Sy(s) + phi_z(c1)*Sz(s))

!---- Flux from cell 2 towards the material interface
      f2 = CONc(material(c2)) *  &
           (phi_x(c2)*Sx(s) + phi_y(c2)*Sy(s) + phi_z(c2)*Sz(s))

!---- The two fluxes (q1 and q2) should be the same
      phi_f = (f2 - f1) / (p1(c1) - p2(c2) + TINY)
      phiside(s) = phi_f

      Dxc1 = xsp(s)-xc(c1)                     
      Dyc1 = ysp(s)-yc(c1)                     
      Dzc1 = zsp(s)-zc(c1)                     
      Dxc2 = xsp(s)-xc(c2)                     
      Dyc2 = ysp(s)-yc(c2)                     
      Dzc2 = zsp(s)-zc(c2)                     

!---- Now update the gradients
      phi_x(c1)=phi_x(c1)+phi_f*(G(1,c1)*Dxc1+G(4,c1)*Dyc1+G(5,c1)*Dzc1)
      phi_y(c1)=phi_y(c1)+phi_f*(G(4,c1)*Dxc1+G(2,c1)*Dyc1+G(6,c1)*Dzc1) 
      phi_z(c1)=phi_z(c1)+phi_f*(G(5,c1)*Dxc1+G(6,c1)*Dyc1+G(3,c1)*Dzc1)

      if(c2 > 0) then
       phi_x(c2)=phi_x(c2)+phi_f*(G(1,c2)*Dxc2+G(4,c2)*Dyc2+G(5,c2)*Dzc2)
       phi_y(c2)=phi_y(c2)+phi_f*(G(4,c2)*Dxc2+G(2,c2)*Dyc2+G(6,c2)*Dzc2)
       phi_z(c2)=phi_z(c2)+phi_f*(G(5,c2)*Dxc2+G(6,c2)*Dyc2+G(3,c2)*Dzc2)
      end if

    end if    
  end do

  call Exchng(phi_x)
  call Exchng(phi_y)
  call Exchng(phi_z)

  end subroutine GraPhiCor
