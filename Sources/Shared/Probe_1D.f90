!==============================================================================!
  subroutine Probe_1D
!------------------------------------------------------------------------------!
!   This subroutine finds the coordinate of cell-centers in non-homogeneous    !
!   direction and write them in file called "name.1D"                          !
!------------------------------------------------------------------------------!
  use all_mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical :: isit
!----------------------------------[Calling]-----------------------------------! 
  include "Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n_prob, p, c
  real              :: zp(1000)
  character(len=80) :: name_prob
  character(len=80) :: answer
!==============================================================================!

  write(*,*) '========================================'
  write(*,*) ' Looking for non-homogeneous directions '
  write(*,*) '----------------------------------------'
  write(*,*) 'Insert non-homogeneous direction (x,y,z or skip)'
  read(*,*) answer
  call touppr(answer)
  if(answer=='SKIP') return

  n_prob = 0
  zp=0.0

  !-----------------------------!
  !   Browse through all cells  !
  !-----------------------------!
  do c=1,NC

    ! Try to find the cell among the probes
    do p=1,n_prob
      if(answer == 'X') then
        if( Approx(xc(c), zp(p), 1.0e-9)) go to 1
      else if(answer == 'Y') then
        if( Approx(yc(c), zp(p), 1.0e-9)) go to 1
      else if(answer == 'Z') then
        if( Approx(zc(c), zp(p), 1.0e-9)) go to 1
      end if
    end do 

    ! Couldn't find a cell among the probes, add a new one
    n_prob = n_prob + 1
    if(answer=='X') zp(n_prob)=xc(c)
    if(answer=='Y') zp(n_prob)=yc(c)
    if(answer=='Z') zp(n_prob)=zc(c)

    if(n_prob == 1000) then
      write(*,*) 'Probe 1D: Not a 1D (channel flow) problem.'
      isit = .false.
      return
    end if
1 end do

  isit = .true.

  !--------------------!
  !   Create 1D file   !
  !--------------------!
  name_prob = name
  name_prob(len_trim(name)+1:len_trim(name)+4) = '.1Dc'
  write(6, *) 'Now creating the file:', name_prob
  open(9, file=name_prob)

  ! Write the number of probes 
  write(9,'(I8)') n_prob

  ! Write the probe coordinates out
  do p=1,n_prob
    write(9,'(I8,1PE17.8)') p, zp(p)
  end do

  close(9)

  end subroutine Probe_1D
