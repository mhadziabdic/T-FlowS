!======================================================================!
  subroutine NearWallCell(grid)
!----------------------------------------------------------------------!
! The subroutine links interior cells to the closes wall cell. This is
! needed for Standard Smagorinsky SGS model used in LES.  
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use rans_mod
  use Grid_Mod
!----------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!-------------------------------[Locals]-------------------------------!
  integer          :: k,  c, nearest_cell  
  real             :: new_distance, old_distance
  real             :: Distance
!======================================================================!  

!--------------------------------------------------------------------------------!
! Purpose: This program is finding for every inside domain cell corresponding cell
! ~~~~~~~~ on the nearest_cell wall. This is needed for calculation y+  
!--------------------------------------------------------------------------------!

  if(this_proc  < 2)                                                     &
  write(*,*) '# Now searching for corresponding wall cells!'

  nearest_cell = 0
  near = 0
  old_distance = HUGE
  do c = 1, NC
    old_distance = HUGE
    do k = 1, NC
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

  if(this_proc < 2) write(*,*) '# Searching finished'

  end subroutine NearWallCell

