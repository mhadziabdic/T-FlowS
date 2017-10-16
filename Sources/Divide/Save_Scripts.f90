!==============================================================================!
  subroutine Save_Scripts
!------------------------------------------------------------------------------!
!   Writes scripts for converting ASCII to binary and vice versa.              !
!------------------------------------------------------------------------------!
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

  open(9, file = 'convert.scr')
  do sub=1,n_sub
    call Name_File(sub, name_out, '.com', len_trim('.com'))
    write(9,'(A8,A80)') './B2A < ', name_out
  end do
  close(9)

  open(9, file = 'CONVERT.SCR')
  do sub=1,n_sub
    call Name_File(sub, name_out, '.com', len_trim('.com'))
    write(9,'(A8,A80)') './A2B < ', name_out
  end do
  close(9)

  end subroutine Save_Scripts
