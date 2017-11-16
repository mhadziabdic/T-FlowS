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
    call Name_File(sub, name_out, '.com', len_trim('.com'))
    open(9, file=name_out)
    write(*,*) '# Now creating the file:', trim(name_out)
    call Name_File(sub, name_out, '.com', 0)
    write(9,*) name_out 
    close(9)
  end do

 end subroutine Save_Com
