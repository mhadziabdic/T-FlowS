!==============================================================================!
  subroutine GradP3(grid, phi, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi. phi may stand either          !
!   for pressure (P) or pressure corrections (PP). This procedure also         !
!   handles different materials.                                               ! 
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phi(-NbC:NC),                              &
                     phi_x(-NbC:NC),                            &
                     phi_y(-NbC:NC),                            &
                     phi_z(-NbC:NC)
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2
  real    :: phi_s, xs, ys, zs 
!==============================================================================!
 
  call Exchng(phi)

  Ps = 0.0

  do c=1,NC
    phi_x(c)=0.0
    phi_y(c)=0.0
    phi_z(c)=0.0
  end do

  !------------------------------------------------------------!
  !   First step: without any wall influence, except outflow   !
  !------------------------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 > 0                            .or. &
       c2 < 0 .and. TypeBC(c2) == BUFFER .or. &
       c2 < 0 .and. TypeBC(c2) == OUTFLOW) then  
      if( StateMat(material(c1))==FLUID .and. &
          StateMat(material(c2))==FLUID ) then  
        phi_s = f(s)*phi(c1)+(1.0-f(s))*phi(c2)  
        phi_x(c1) = phi_x(c1) + phi_s * grid % sx(s)
        phi_y(c1) = phi_y(c1) + phi_s * grid % sy(s)
        phi_z(c1) = phi_z(c1) + phi_s * grid % sz(s)
        phi_x(c2) = phi_x(c2) - phi_s * grid % sx(s)
        phi_y(c2) = phi_y(c2) - phi_s * grid % sy(s)
        phi_z(c2) = phi_z(c2) - phi_s * grid % sz(s)
      end if
    end if
  end do

  do c=1,NC
    if(StateMat(material(c))==FLUID) then
      phi_x(c)=phi_x(c)/volume(c)
      phi_y(c)=phi_y(c)/volume(c)
      phi_z(c)=phi_z(c)/volume(c)
    end if
  end do

  !------------------------------------------------------------!
  !   Second step: extrapolate to boundaries, except outflow   !
  !------------------------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0               .and. &
       TypeBC(c2) /= BUFFER .and. &
       TypeBC(c2) /= OUTFLOW) then  
      if(StateMat(material(c1))==FLUID) then
        Ps(s) = (   phi(c1)                     +     &
                    phi_x(c1) * (grid % xc(c2)-grid % xc(c1)) +     &
                    phi_y(c1) * (grid % yc(c2)-grid % yc(c1)) +     &
                    phi_z(c1) * (grid % zc(c2)-grid % zc(c1))   ) / &
              ( 1.0 - (  grid % sx(s) * (grid % xc(c2)-grid % xc(c1))      & 
                       + grid % sy(s) * (grid % yc(c2)-grid % yc(c1))      &
                       + grid % sz(s) * (grid % zc(c2)-grid % zc(c1))  ) / volume(c1)  )
      end if
    end if

    ! Handle two materials
    if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER) then  
      if( StateMat(material(c1))==FLUID .and. &
          StateMat(material(c2))==SOLID ) then  
        xs = grid % xf(s) 
        ys = grid % yf(s) 
        zs = grid % zf(s) 
        Ps(s) = (   phi(c1)                 +     &
                    phi_x(c1) * (xs-grid % xc(c1)) +     &
                    phi_y(c1) * (ys-grid % yc(c1)) +     &
                    phi_z(c1) * (zs-grid % zc(c1))   ) / &
              ( 1.0 - (  grid % sx(s) * (xs-grid % xc(c1))      & 
                       + grid % sy(s) * (ys-grid % yc(c1))      &
                       + grid % sz(s) * (zs-grid % zc(c1))  ) / volume(c1)  )
      end if
      if( StateMat(material(c1))==SOLID .and. &
          StateMat(material(c2))==FLUID ) then  
        xs = grid % xf(s) 
        ys = grid % yf(s) 
        zs = grid % zf(s) 
        Ps(s) = (   phi(c2)                 +     &
                    phi_x(c2) * (xs-grid % xc(c2)) +     &
                    phi_y(c2) * (ys-grid % yc(c2)) +     &
                    phi_z(c2) * (zs-grid % zc(c2))   ) / &
              ( 1.0 + (  grid % sx(s) * (xs-grid % xc(c2))      & 
                       + grid % sy(s) * (ys-grid % yc(c2))      &
                       + grid % sz(s) * (zs-grid % zc(c2))  ) / volume(c2)  )
      end if
    end if ! c2 < 0
  end do

  !---------------------------------------------!
  !   Third step: compute the final gradients   !
  !---------------------------------------------!
  do s=1,NS
    c1=SideC(1,s)
    c2=SideC(2,s)
    if(c2 < 0               .and. &
       TypeBC(c2) /= BUFFER .and. &
       TypeBC(c2) /= OUTFLOW) then  
      phi_x(c1) = phi_x(c1) + Ps(s) * grid % sx(s)/volume(c1)
      phi_y(c1) = phi_y(c1) + Ps(s) * grid % sy(s)/volume(c1)
      phi_z(c1) = phi_z(c1) + Ps(s) * grid % sz(s)/volume(c1)
    end if

    ! Handle two materials
    if(c2 > 0 .or. c2 < 0 .and. TypeBC(c2) == BUFFER) then  
      if( StateMat(material(c1))==FLUID .and. &
          StateMat(material(c2))==SOLID ) then  
        phi_x(c1) = phi_x(c1) + Ps(s) * grid % sx(s)/volume(c1)
        phi_y(c1) = phi_y(c1) + Ps(s) * grid % sy(s)/volume(c1)
        phi_z(c1) = phi_z(c1) + Ps(s) * grid % sz(s)/volume(c1)
      end if 
      if( StateMat(material(c1))==SOLID .and. &
          StateMat(material(c2))==FLUID ) then  
        phi_x(c2) = phi_x(c2) - Ps(s) * grid % sx(s)/volume(c2)
        phi_y(c2) = phi_y(c2) - Ps(s) * grid % sy(s)/volume(c2)
        phi_z(c2) = phi_z(c2) - Ps(s) * grid % sz(s)/volume(c2)
      end if 
    end if  ! c2 < 0
  end do

  call Exchng(phi_x)
  call Exchng(phi_y)
  call Exchng(phi_z)

  end subroutine
