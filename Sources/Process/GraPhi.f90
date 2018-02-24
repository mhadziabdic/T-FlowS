!==============================================================================!
  subroutine GraPhi(grid, phi, i, phii, boundary)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by a least squares method.     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Flow_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: i
  real            :: phi (-grid % n_bnd_cells:grid % n_cells),  &
                     phii(-grid % n_bnd_cells:grid % n_cells) 
  logical         :: boundary
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2
  real    :: dphi1, dphi2, dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2 
!==============================================================================!

  call Comm_Mod_Exchange(grid, phi)

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
      dphi1 = phi(c2)-phi(c1) 
      dphi2 = phi(c2)-phi(c1) 
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY) then
        dphi1 = 0.0  
        dphi2 = 0.0  
      end if

      if(boundary) then
        phii(c1)=phii(c1)+dphi1*(g(1,c1)*dx_c1+g(4,c1)*dy_c1+g(5,c1)*dz_c1) 
        if(c2 > 0) then
          phii(c2)=phii(c2)+dphi2*(g(1,c2)*dx_c2+g(4,c2)*dy_c2+g(5,c2)*dz_c2)
        end if
      else
        if(c2 > 0 .or.  &
           c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER ) & 
             phii(c1)=phii(c1)+dphi1*(g(1,c1)*dx_c1+g(4,c1)*dy_c1+g(5,c1)*dz_c1) 
        if(c2 > 0)  &
             phii(c2)=phii(c2)+dphi2*(g(1,c2)*dx_c2+g(4,c2)*dy_c2+g(5,c2)*dz_c2)
      end if ! boundary
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
      dphi1 = phi(c2)-phi(c1) 
      dphi2 = phi(c2)-phi(c1) 
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY) then
        dphi1 = 0.0  
        dphi2 = 0.0  
      end if

      if(boundary) then
        phii(c1)=phii(c1)+dphi1*(g(4,c1)*dx_c1+g(2,c1)*dy_c1+g(6,c1)*dz_c1) 
        if(c2  > 0) then
          phii(c2)=phii(c2)+dphi2*(g(4,c2)*dx_c2+g(2,c2)*dy_c2+g(6,c2)*dz_c2)
        end if
      else
        if(c2 > 0 .or. c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER ) & 
                    phii(c1) = phii(c1) + dphi1*(  g(4,c1)*dx_c1  &
                                                 + g(2,c1)*dy_c1  &
                                                 + g(6,c1)*dz_c1) 
        if(c2  > 0) phii(c2) = phii(c2) + dphi2*(  g(4,c2)*dx_c2  &
                                                 + g(2,c2)*dy_c2  &
                                                 + g(6,c2)*dz_c2)
      end if ! boundary
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
      dphi1 = phi(c2)-phi(c1) 
      dphi2 = phi(c2)-phi(c1) 
      if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY) then
        dphi1 = 0.0  
        dphi2 = 0.0  
      end if

      if(boundary) then
        phii(c1)=phii(c1)+dphi1*(g(5,c1)*dx_c1+g(6,c1)*dy_c1+g(3,c1)*dz_c1) 
        if(c2 > 0) then
          phii(c2)=phii(c2)+dphi2*(g(5,c2)*dx_c2+g(6,c2)*dy_c2+g(3,c2)*dz_c2)
        end if
      else
        if(c2 > 0 .or.  &
           c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER ) & 
                   phii(c1) = phii(c1) + dphi1 * (  g(5,c1)*dx_c1  &
                                                  + g(6,c1)*dy_c1  &
                                                  + g(3,c1)*dz_c1 ) 
        if(c2 > 0) phii(c2) = phii(c2) + dphi2 * (  g(5,c2)*dx_c2  &
                                                  + g(6,c2)*dy_c2  &
                                                  + g(3,c2)*dz_c2) 
      end if ! boundary
    end do
  end if

  call Comm_Mod_Exchange(grid, phii)

  if(.not. boundary) call Correct_Bad(grid, phii)

  end subroutine
