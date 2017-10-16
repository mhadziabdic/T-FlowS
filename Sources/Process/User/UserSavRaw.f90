!======================================================================!
  subroutine UserSavRaw
!----------------------------------------------------------------------!
! Writes: NAME.rawdata                                                 !
! ~~~~~~~                                                              !
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer   :: c, NCtot
  character :: namOut*80
!======================================================================!

!------------------------------!
!     Create raw data file     !
!------------------------------!
  namOut = name
  call NamFil(this, namOut, '.rawdata', len_trim('.rawdata') )
  open(9, file=namOut)
  write(6, *) 'Now creating the file:', namOut

!---- Total number of cells
  NCtot = NC
  call IGlSum(NCtot)
  if(this < 2) then
    write(9,'(I9)') NCtot
  end if 

!---- Write raw data file 
  do c=1,NC
    write(9,'(6E16.6)') xc(c),yc(c),zc(c),U % n(c),V % n(c),W % n(c)
  end do

  close(9)

  end subroutine UserSavRaw
