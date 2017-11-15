!==============================================================================!
  subroutine Compute_Minimum_Maximum(grid, phi)
!------------------------------------------------------------------------------!
!                                                                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phi(-grid % n_bnd_cells:grid % n_cells) 
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
!==============================================================================!

  phi_max = phi 
  phi_min = phi 

  do s = 1, grid % n_faces
    c1=SideC(1,s)
    c2=SideC(2,s)

    if( (c2>0) .or. (c2<0 .and. TypeBC(c2)==BUFFER) ) then
      phi_max(c1) = max(phi_max(c1), phi(c2))
      phi_min(c1) = min(phi_min(c1), phi(c2))
      phi_max(c2) = max(phi_max(c2), phi(c1))
      phi_min(c2) = min(phi_min(c2), phi(c1))
    end if

  end do

  call Exchange(grid, phi_max)
  call Exchange(grid, phi_min)

  end subroutine
