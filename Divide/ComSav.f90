!======================================================================!
  SUBROUTINE ComSav()
!----------------------------------------------------------------------!
! Writes: *.com, convert.scr and CONVERT.scr                           !
! ~~~~~~~                                                              ! 
!------------------------------[Modules]-------------------------------!
  USE all_mod
  USE gen_mod
  USE par_mod
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------[Locals]-------------------------------!
  CHARACTER :: namOut*80
  INTEGER   :: sub
!--------------------------------[CVS]---------------------------------!
!  $Id: ComSav.f90,v 1.1 2014/11/24 11:32:32 muhamed Exp $  
!  $Source: /home/mhadziabdic/Dropbox/cvsroot/T-FlowS-CVS/Divide/ComSav.f90,v $             
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

  END SUBROUTINE ComSav
