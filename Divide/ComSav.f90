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
  character :: namOut*80
  integer   :: sub
!======================================================================!

!<<<<<<<<<<<<<<<<<<<<<<<<<!
!                         !
!     create com file     !
!                         !
!<<<<<<<<<<<<<<<<<<<<<<<<<!
  do sub=1,Nsub 
    call NamFil(sub, namOut, '.com', len_trim('.com'))
    open(9, FILE=namOut)
    write(6, *) 'Now creating the file:', namOut
    call NamFil(sub, namOut, '.com', 0)
    write(9,*) namOut 
    close(9)
  end do

  open(9, FILE = 'convert.scr')
  do sub=1,Nsub
    call NamFil(sub, namOut, '.com', len_trim('.com'))
    write(9,'(A8,A80)') './B2A < ', namOut
  end do
  close(9)

  open(9, FILE = 'CONVERT.SCR')
  do sub=1,Nsub
    call NamFil(sub, namOut, '.com', len_trim('.com'))
    write(9,'(A8,A80)') './A2B < ', namOut
  end do
  close(9)

  end subroutine ComSav
