!==============================================================================!
  subroutine Save_Com()
!------------------------------------------------------------------------------!
! Writes: *.com                                                                !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use div_mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_out
  integer           :: sub
!==============================================================================!

  do sub=1,n_sub 
    call Name_File(sub, name_out, '.com')
    open(9, file=name_out)
    write(*,*) '# Creating the file:', trim(name_out)
    write(9,*) name_out(1:len_trim(name_out)-4) 
    close(9)
  end do

 end subroutine
