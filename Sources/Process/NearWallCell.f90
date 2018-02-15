!==============================================================================!
  subroutine NearWallCell(grid)
!------------------------------------------------------------------------------!
! The subroutine links interior cells to the closes wall cell. This is         !
! needed for Standard Smagorinsky SGS model used in 'LES'.                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer          :: k,  c, nearest_cell  
  real             :: new_distance, old_distance
  real             :: Distance
!==============================================================================!  

  if(this_proc  < 2)                                                     &
  print *, '# Searching for corresponding wall cells!'

  nearest_cell = 0
  near = 0
  old_distance = HUGE
  do c = 1, grid % n_cells
    old_distance = HUGE
    do k = 1, grid % n_cells
      if(IsNearWall(k)) then
        new_distance = Distance(grid % xc(k),grid % yc(k),grid % zc(k),grid % xc(c),grid % yc(c),grid % zc(c))
        if(new_distance <= old_distance) then
          nearest_cell =  k
          old_distance = new_distance
        end if 
      end if
    end do
    near(c) = nearest_cell 
  end do

  if(this_proc < 2) print *, '# Searching finished'

  end subroutine

