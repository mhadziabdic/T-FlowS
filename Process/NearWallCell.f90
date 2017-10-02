!======================================================================!
  subroutine NearWallCell()
!----------------------------------------------------------------------!
! The subroutine links interior cells to the closes wall cell. This is
! needed for Standard Smagorinsky SGS model used in LES.  
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer          :: k,  c, nearest  
  real             :: DISTnew, DISTold
  real             :: Dist
!======================================================================!  

!--------------------------------------------------------------------------------!
! Purpose: This program is finding for every inside domain cell corresponding cell
! ~~~~~~~~ on the nearest wall. This is needed for calculation y+  
!--------------------------------------------------------------------------------!

  if(this  < 2)                                                     &
  write(*,*) '# Now searching for corresponding wall cells!'

  nearest = 0
  near = 0
  DISTold = HUGE
  do c = 1, NC
    DISTold = HUGE
    do k = 1, NC
      if(IsNearWall(k)) then
        DISTnew = Dist(xc(k),yc(k),zc(k),xc(c),yc(c),zc(c))
        if(DISTnew <= DISTold) then
          nearest =  k
          DISTold = DISTnew
        end if 
      end if
    end do
    near(c) = nearest 
  end do

  if(this < 2) write(*,*) '# Searching finished'

  end subroutine NearWallCell

