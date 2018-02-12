!==============================================================================!
  subroutine GraPhi(grid, phi, i, phii, Boundary)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by a least squares method.     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: i
  real            :: phi (-grid % n_bnd_cells:grid % n_cells),  &
                     phii(-grid % n_bnd_cells:grid % n_cells) 
  logical         :: Boundary
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2
  real    :: Dphi1, Dphi2, dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2 
!==============================================================================!

  call Exchange(grid, phi)

  do c=1,grid % n_cells
    phii(c)=0.0
  end do

  if(i == 1) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      dx_c1 = grid % dx(s)
      dy_c1 = grid % dy(s)
      dz_c1 = grid % dz(s)
      dx_c2 = grid % dx(s)
      dy_c2 = grid % dy(s)
      dz_c2 = grid % dz(s)
      Dphi1 = phi(c2)-phi(c1) 
      Dphi2 = phi(c2)-phi(c1) 
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY) then
        Dphi1 = 0.0  
        Dphi2 = 0.0  
      end if

!---- Take care of material interfaces           ! 2mat
!      if( StateMat(material(c1))==FLUID .and. &  ! 2mat
!          StateMat(material(c2))==SOLID       &  ! 2mat
!          .or.                                &  ! 2mat
!          StateMat(material(c1))==SOLID .and. &  ! 2mat
!          StateMat(material(c2))==FLUID ) then   ! 2mat
!        dx_c1 = xsp(s)-xc(c1)                     ! 2mat
!        dy_c1 = ysp(s)-yc(c1)                     ! 2mat
!        dz_c1 = zsp(s)-zc(c1)                     ! 2mat
!        dx_c2 = xsp(s)-xc(c2)                     ! 2mat
!        dy_c2 = ysp(s)-yc(c2)                     ! 2mat
!        dz_c2 = zsp(s)-zc(c2)                     ! 2mat
!        Dphi1 = 0.0-phi(c1)                      ! 2mat
!        Dphi2 = 0.0-phi(c2)                      ! 2mat 
!      end if                                     ! 2mat
              
      if(Boundary) then
        phii(c1)=phii(c1)+Dphi1*(g(1,c1)*dx_c1+g(4,c1)*dy_c1+g(5,c1)*dz_c1) 
        if(c2 > 0) then
          phii(c2)=phii(c2)+Dphi2*(g(1,c2)*dx_c2+g(4,c2)*dy_c2+g(5,c2)*dz_c2)
        end if
      else
        if(c2 > 0 .or.  &
           c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER ) & 
             phii(c1)=phii(c1)+Dphi1*(g(1,c1)*dx_c1+g(4,c1)*dy_c1+g(5,c1)*dz_c1) 
        if(c2 > 0)  &
             phii(c2)=phii(c2)+Dphi2*(g(1,c2)*dx_c2+g(4,c2)*dy_c2+g(5,c2)*dz_c2)
      end if ! Boundary
    end do
  end if

  if(i == 2) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      dx_c1 = grid % dx(s)
      dy_c1 = grid % dy(s)
      dz_c1 = grid % dz(s)
      dx_c2 = grid % dx(s)
      dy_c2 = grid % dy(s)
      dz_c2 = grid % dz(s)
      Dphi1 = phi(c2)-phi(c1) 
      Dphi2 = phi(c2)-phi(c1) 
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY) then
        Dphi1 = 0.0  
        Dphi2 = 0.0  
      end if

      ! Take care of material interfaces            ! 2mat
      !  if( StateMat(material(c1))==FLUID .and. &  ! 2mat
      !      StateMat(material(c2))==SOLID       &  ! 2mat
      !      .or.                                &  ! 2mat
      !      StateMat(material(c1))==SOLID .and. &  ! 2mat
      !      StateMat(material(c2))==FLUID ) then   ! 2mat
      !    dx_c1 = xsp(s)-xc(c1)                     ! 2mat
      !    dy_c1 = ysp(s)-yc(c1)                     ! 2mat
      !    dz_c1 = zsp(s)-zc(c1)                     ! 2mat
      !    dx_c2 = xsp(s)-xc(c2)                     ! 2mat
      !    dy_c2 = ysp(s)-yc(c2)                     ! 2mat
      !    dz_c2 = zsp(s)-zc(c2)                     ! 2mat
      !    Dphi1 = 0.0-phi(c1)                      ! 2mat
      !    Dphi2 = 0.0-phi(c2)                      ! 2mat 
      !  end if                                     ! 2mat

      if(Boundary) then
        phii(c1)=phii(c1)+Dphi1*(g(4,c1)*dx_c1+g(2,c1)*dy_c1+g(6,c1)*dz_c1) 
        if(c2  > 0) then
          phii(c2)=phii(c2)+Dphi2*(g(4,c2)*dx_c2+g(2,c2)*dy_c2+g(6,c2)*dz_c2)
        end if
      else
        if(c2 > 0 .or. c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER ) & 
                    phii(c1) = phii(c1) + Dphi1*(  g(4,c1)*dx_c1  &
                                                 + g(2,c1)*dy_c1  &
                                                 + g(6,c1)*dz_c1) 
        if(c2  > 0) phii(c2) = phii(c2) + Dphi2*(  g(4,c2)*dx_c2  &
                                                 + g(2,c2)*dy_c2  &
                                                 + g(6,c2)*dz_c2)
      end if ! Boundary
    end do
  end if

  if(i == 3) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      dx_c1 = grid % dx(s)
      dy_c1 = grid % dy(s)
      dz_c1 = grid % dz(s)
      dx_c2 = grid % dx(s)
      dy_c2 = grid % dy(s)
      dz_c2 = grid % dz(s)
      Dphi1 = phi(c2)-phi(c1) 
      Dphi2 = phi(c2)-phi(c1) 
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY) then
        Dphi1 = 0.0  
        Dphi2 = 0.0  
      end if

      ! Take care of material interfaces            ! 2mat
      !  if( StateMat(material(c1))==FLUID .and. &  ! 2mat
      !      StateMat(material(c2))==SOLID       &  ! 2mat
      !      .or.                                &  ! 2mat
      !      StateMat(material(c1))==SOLID .and. &  ! 2mat
      !      StateMat(material(c2))==FLUID ) then   ! 2mat
      !    dx_c1 = xsp(s)-xc(c1)                     ! 2mat
      !    dy_c1 = ysp(s)-yc(c1)                     ! 2mat
      !    dz_c1 = zsp(s)-zc(c1)                     ! 2mat
      !    dx_c2 = xsp(s)-xc(c2)                     ! 2mat
      !    dy_c2 = ysp(s)-yc(c2)                     ! 2mat
      !    dz_c2 = zsp(s)-zc(c2)                     ! 2mat
      !    Dphi1 = 0.0-phi(c1)                      ! 2mat
      !    Dphi2 = 0.0-phi(c2)                      ! 2mat 
      !  end if                                     ! 2mat 

      if(Boundary) then
        phii(c1)=phii(c1)+Dphi1*(g(5,c1)*dx_c1+g(6,c1)*dy_c1+g(3,c1)*dz_c1) 
        if(c2 > 0) then
          phii(c2)=phii(c2)+Dphi2*(g(5,c2)*dx_c2+g(6,c2)*dy_c2+g(3,c2)*dz_c2)
        end if
      else
        if(c2 > 0 .or.  &
           c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER ) & 
                   phii(c1) = phii(c1) + Dphi1 * (  g(5,c1)*dx_c1  &
                                                  + g(6,c1)*dy_c1  &
                                                  + g(3,c1)*dz_c1 ) 
        if(c2 > 0) phii(c2) = phii(c2) + Dphi2 * (  g(5,c2)*dx_c2  &
                                                  + g(6,c2)*dy_c2  &
                                                  + g(3,c2)*dz_c2) 
      end if ! Boundary
    end do
  end if

  call Exchange(grid, phii)

  if(.not. Boundary) call Correct_Bad(grid, phii)

  end subroutine
