!======================================================================!
  subroutine BCelLoa
!----------------------------------------------------------------------!
! Reads:  NAME.faces.gmv  NAME.shadow.gmv                              !
! ~~~~~~                                                               !                                                                           
!------------------------------[Modules]-------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer           :: c1, c2, n, s, dum_i
  character(len=80) :: dum_s, name_in
!======================================================================!

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!     Read the GMV file with      !
!        cell face nodes          !
!                                 !
!     (it has been done by        !
!      copying and modifying      ! 
!          GenSav.f90)            !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  name_in = name
  name_in(len_trim(name)+1:len_trim(name)+10) = '.faces.gmv'
  write(*,*) '# Now reading the file:', name_in
  open(9, FILE=name_in)

!---------------!
!     start     !
!---------------!
  read(9,'(A80)') dum_s 

!---------------!
!     nodes     !
!---------------!
  read(9,'(A80)') dum_s 
  do n=1,3*NN
    read(9,'(A80)') dum_s 
  end do  

!----------------------!
!     cell section     !
!----------------------!
  read(9,'(A80)') dum_s 
  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    read(9,*) dum_s, dum_i
    if(dum_s == 'tri') then 
      SideN(s,0) = 3
      read(9,*) &
        SideN(s,1), SideN(s,2), SideN(s,3)
    else if(dum_s == 'quad') then
      SideN(s,0) = 4
      read(9,*) &
        SideN(s,1), SideN(s,2), SideN(s,3), SideN(s,4)  
    else
      write(*,*) 'Unsupported cell-face type:', dum_s
      write(*,*) 'Exiting'
      stop
    end if
  end do  

  close(9)

!>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                          !
!     read shadow file     !
!                          !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  call NamFil(0, name_in, '.shadow.gmv', len_trim('.shadow.gmv'))
  open(9, FILE=name_in)
  write(6, *) 'Now reading the file:', name_in

  do s=NS+1,NS+NSsh
    read(9,*) SideN(s,0)
    if(SideN(s,0)==3) then
      read(9,*) SideN(s,1), SideN(s,2), SideN(s,3),             &
                SideC(1,s), SideC(2,s)
    else if(SideN(s,0)==4) then
      read(9,*) SideN(s,1), SideN(s,2), SideN(s,3), SideN(s,4), &
                SideC(1,s), SideC(2,s)
    end if
  end do

  close(9)

  end subroutine BCelLoa
