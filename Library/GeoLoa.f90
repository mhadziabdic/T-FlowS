!======================================================================!
  subroutine GeoLoa
!----------------------------------------------------------------------!
! Reads:  NAME.geo                                                     !
! ~~~~~~                                                               !
!------------------------------[Modules]-------------------------------!
  use all_mod
  use par_mod  ! needs to know THIS
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer      :: c, s
  character*80 :: nameIn  
!======================================================================!

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!       Read the file with       !
!     geometrical quantities     !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call NamFil(THIS, nameIn, '.geo', len_trim('.geo')) 
  open(9, FILE=nameIn, FORM='UNFORMATTED')
  if(this < 2) write(*,*) '# Now reading the file:', nameIn

  read(9) (xc(c), c=1,NC)
  read(9) (yc(c), c=1,NC) 
  read(9) (zc(c), c=1,NC)

  read(9) (xc(c), c=-1,-NBC,-1)  
  read(9) (yc(c), c=-1,-NBC,-1)
  read(9) (zc(c), c=-1,-NBC,-1) 

  read(9) (volume(c), c=1,NC)
  read(9) (delta(c),  c=1,NC)

  read(9) (WallDs(c), c=1,NC)

  read(9) (Sx(s), s=1,NS)
  read(9) (Sy(s), s=1,NS)
  read(9) (Sz(s), s=1,NS)

  read(9) (Dx(s), s=1,NS)
  read(9) (Dy(s), s=1,NS)
  read(9) (Dz(s), s=1,NS)

  read(9) (f(s), s=1,NS)

  read(9) (xsp(s), s=1,NS)
  read(9) (ysp(s), s=1,NS)
  read(9) (zsp(s), s=1,NS)

  close(9) 

!->>> do s=1,NS
!->>>   write(*,*) '================'
!->>>	write(*,*) s
!->>>	c1 = sidec(1,s)
!->>>	c2 = sidec(2,s)
!->>>	write(*,'(3F10.5)') xc(c1),yc(c1),zc(c1)  
!->>>	write(*,'(3F10.5)') xc(c2),yc(c2),zc(c2)  
!->>>	write(*,'(3F10.5)') Sx(s),Sy(s),Sz(s)  
!->>> end do

  end subroutine GeoLoa
