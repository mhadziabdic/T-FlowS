!======================================================================!
  subroutine GraPhi(phi, i, phii, Boundary)
!----------------------------------------------------------------------!
! Calculates gradient of generic variable phi by a least squares       !
! method.                                                              !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  integer :: i
  real    :: phi(-NbC:NC), phii(-NbC:NC) 
  logical :: Boundary
!-------------------------------[Locals]-------------------------------!
  integer :: s, c, c1, c2
  real    :: Dphi1, Dphi2, Dxc1, Dyc1, Dzc1, Dxc2, Dyc2, Dzc2 
!======================================================================!

  call Exchng(phi)

  do c=1,NC
    phii(c)=0.0
  end do

  if(i == 1) then
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      Dxc1 = Dx(s)
      Dyc1 = Dy(s)
      Dzc1 = Dz(s)
      Dxc2 = Dx(s)
      Dyc2 = Dy(s)
      Dzc2 = Dz(s)
      Dphi1 = phi(c2)-phi(c1) 
      Dphi2 = phi(c2)-phi(c1) 
      if(c2 < 0 .and. TypeBC(c2) == SYMMETRY) then
        Dphi1 = 0.0  
        Dphi2 = 0.0  
      end if

!---- Take care of material interfaces           ! 2mat
!      if( StateMat(material(c1))==FLUID .and. &  ! 2mat
!          StateMat(material(c2))==SOLID       &  ! 2mat
!          .or.                                &  ! 2mat
!          StateMat(material(c1))==SOLID .and. &  ! 2mat
!          StateMat(material(c2))==FLUID ) then   ! 2mat
!        Dxc1 = xsp(s)-xc(c1)                     ! 2mat
!        Dyc1 = ysp(s)-yc(c1)                     ! 2mat
!        Dzc1 = zsp(s)-zc(c1)                     ! 2mat
!        Dxc2 = xsp(s)-xc(c2)                     ! 2mat
!        Dyc2 = ysp(s)-yc(c2)                     ! 2mat
!        Dzc2 = zsp(s)-zc(c2)                     ! 2mat
!        Dphi1 = 0.0-phi(c1)                      ! 2mat
!        Dphi2 = 0.0-phi(c2)                      ! 2mat 
!      end if                                     ! 2mat
              
      if(Boundary) then
        phii(c1)=phii(c1)+Dphi1*(G(1,c1)*Dxc1+G(4,c1)*Dyc1+G(5,c1)*Dzc1) 
        if(c2 > 0) then
	  phii(c2)=phii(c2)+Dphi2*(G(1,c2)*Dxc2+G(4,c2)*Dyc2+G(5,c2)*Dzc2)
        end if
      else
        if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER ) & 
                    phii(c1)=phii(c1)+Dphi1*(G(1,c1)*Dxc1+G(4,c1)*Dyc1+G(5,c1)*Dzc1) 
	if(c2 > 0)  phii(c2)=phii(c2)+Dphi2*(G(1,c2)*Dxc2+G(4,c2)*Dyc2+G(5,c2)*Dzc2)
      end if ! Boundary
    end do
  end if

  if(i == 2) then
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      Dxc1 = Dx(s)
      Dyc1 = Dy(s)
      Dzc1 = Dz(s)
      Dxc2 = Dx(s)
      Dyc2 = Dy(s)
      Dzc2 = Dz(s)
      Dphi1 = phi(c2)-phi(c1) 
      Dphi2 = phi(c2)-phi(c1) 
      if(c2 < 0 .and. TypeBC(c2) == SYMMETRY) then
        Dphi1 = 0.0  
        Dphi2 = 0.0  
      end if

!---- Take care of material interfaces           ! 2mat
!      if( StateMat(material(c1))==FLUID .and. &  ! 2mat
!          StateMat(material(c2))==SOLID       &  ! 2mat
!          .or.                                &  ! 2mat
!          StateMat(material(c1))==SOLID .and. &  ! 2mat
!          StateMat(material(c2))==FLUID ) then   ! 2mat
!        Dxc1 = xsp(s)-xc(c1)                     ! 2mat
!        Dyc1 = ysp(s)-yc(c1)                     ! 2mat
!        Dzc1 = zsp(s)-zc(c1)                     ! 2mat
!        Dxc2 = xsp(s)-xc(c2)                     ! 2mat
!        Dyc2 = ysp(s)-yc(c2)                     ! 2mat
!        Dzc2 = zsp(s)-zc(c2)                     ! 2mat
!        Dphi1 = 0.0-phi(c1)                      ! 2mat
!        Dphi2 = 0.0-phi(c2)                      ! 2mat 
!      end if                                     ! 2mat

      if(Boundary) then
        phii(c1)=phii(c1)+Dphi1*(G(4,c1)*Dxc1+G(2,c1)*Dyc1+G(6,c1)*Dzc1) 
        if(c2  > 0) then
	  phii(c2)=phii(c2)+Dphi2*(G(4,c2)*Dxc2+G(2,c2)*Dyc2+G(6,c2)*Dzc2)
        end if
      else
        if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER ) & 
                    phii(c1)=phii(c1)+Dphi1*(G(4,c1)*Dxc1+G(2,c1)*Dyc1+G(6,c1)*Dzc1) 
        if(c2  > 0) phii(c2)=phii(c2)+Dphi2*(G(4,c2)*Dxc2+G(2,c2)*Dyc2+G(6,c2)*Dzc2)
      end if ! Boundary
    end do
  end if

  if(i == 3) then
    do s=1,NS
      c1=SideC(1,s)
      c2=SideC(2,s)
      Dxc1 = Dx(s)
      Dyc1 = Dy(s)
      Dzc1 = Dz(s)
      Dxc2 = Dx(s)
      Dyc2 = Dy(s)
      Dzc2 = Dz(s)
      Dphi1 = phi(c2)-phi(c1) 
      Dphi2 = phi(c2)-phi(c1) 
      if(c2 < 0 .and. TypeBC(c2) == SYMMETRY) then
        Dphi1 = 0.0  
        Dphi2 = 0.0  
      end if

!---- Take care of material interfaces           ! 2mat
!      if( StateMat(material(c1))==FLUID .and. &  ! 2mat
!          StateMat(material(c2))==SOLID       &  ! 2mat
!          .or.                                &  ! 2mat
!          StateMat(material(c1))==SOLID .and. &  ! 2mat
!          StateMat(material(c2))==FLUID ) then   ! 2mat
!        Dxc1 = xsp(s)-xc(c1)                     ! 2mat
!        Dyc1 = ysp(s)-yc(c1)                     ! 2mat
!        Dzc1 = zsp(s)-zc(c1)                     ! 2mat
!        Dxc2 = xsp(s)-xc(c2)                     ! 2mat
!        Dyc2 = ysp(s)-yc(c2)                     ! 2mat
!        Dzc2 = zsp(s)-zc(c2)                     ! 2mat
!        Dphi1 = 0.0-phi(c1)                      ! 2mat
!        Dphi2 = 0.0-phi(c2)                      ! 2mat 
!      end if                                     ! 2mat 

      if(Boundary) then
        phii(c1)=phii(c1)+Dphi1*(G(5,c1)*Dxc1+G(6,c1)*Dyc1+G(3,c1)*Dzc1) 
        if(c2 > 0) then
	  phii(c2)=phii(c2)+Dphi2*(G(5,c2)*Dxc2+G(6,c2)*Dyc2+G(3,c2)*Dzc2)
        end if
      else
        if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER ) & 
                   phii(c1)=phii(c1)+Dphi1*(G(5,c1)*Dxc1+G(6,c1)*Dyc1+G(3,c1)*Dzc1) 
        if(c2 > 0) phii(c2)=phii(c2)+Dphi2*(G(5,c2)*Dxc2+G(6,c2)*Dyc2+G(3,c2)*Dzc2) 
      end if ! Boundary
    end do
  end if

  call Exchng(phii)

  if(.not. Boundary) call CorBad(phii)

  end subroutine GraPhi
