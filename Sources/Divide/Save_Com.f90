!==============================================================================!
  subroutine Save_Com()
!------------------------------------------------------------------------------!
! Writes: *.com                                                                !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use par_mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_out
  integer           :: sub
!==============================================================================!

  do sub=1,n_sub 
    call NamFil(sub, name_out, '.com', len_trim('.com'))
    open(9, file=name_out)
    write(6, *) 'Now creating the file:', name_out
    call NamFil(sub, name_out, '.com', 0)
    write(9,*) name_out 
    close(9)
  end do

 end subroutine Save_Com
