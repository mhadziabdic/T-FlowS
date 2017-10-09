!======================================================================!
  subroutine ComSav()
!----------------------------------------------------------------------!
! Writes: *.com, convert.scr and CONVERT.scr                           !
! ~~~~~~~                                                              ! 
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  character(len=80) :: name_out
  integer           :: sub
!======================================================================!

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!                         !
!     create com file     !
!                         !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  do sub=1,n_sub 
    call NamFil(sub, name_out, '.com', len_trim('.com'))
    open(9, FILE=name_out)
    write(6, *) 'Now creating the file:', name_out
    call NamFil(sub, name_out, '.com', 0)
    write(9,*) name_out 
    close(9)
  end do

  open(9, FILE = 'convert.scr')
  do sub=1,n_sub
    call NamFil(sub, name_out, '.com', len_trim('.com'))
    write(9,'(A8,A80)') './B2A < ', name_out
  end do
  close(9)

  open(9, FILE = 'CONVERT.SCR')
  do sub=1,n_sub
    call NamFil(sub, name_out, '.com', len_trim('.com'))
    write(9,'(A8,A80)') './A2B < ', name_out
  end do
  close(9)

  end subroutine ComSav
