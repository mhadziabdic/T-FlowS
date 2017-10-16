!==============================================================================!
  subroutine Load_Geo
!------------------------------------------------------------------------------!
!   Reads:  name.geo                                                           !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use par_mod, only: this_proc
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, s
  character(len=80) :: name_in  
!==============================================================================!

  !----------------------------!
  !     Read the file with     !
  !   geometrical quantities   !
  !----------------------------!
  call Name_File(this_proc, name_in, '.geo', len_trim('.geo')) 
  open(9, file=name_in, FORM='unformatted')
  if(this_proc < 2) write(*,*) '# Now reading the file:', name_in

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

  end subroutine Load_Geo
