!======================================================================!
  subroutine Probe_2D 
!----------------------------------------------------------------------!
! Finds coordinates of all the planes for the channel flow.            !
! It assumes that homogeneous directions of the flow are x and y.      !
!----------------------------------------------------------------------!
  use all_mod
!----------------------------------------------------------------------!
  implicit none
!------------------------------[Calling]-------------------------------! 
  include "Approx.int"
!-------------------------------[Locals]-------------------------------!
  integer           :: n_prob, p, c
  real              :: yp(20000), zp(20000)
  character(len=80) :: name_prob
  character(len=80) :: answer
!======================================================================!

  write(*,*) '==============================='
  write(*,*) ' Looking for homogeneous plane '
  write(*,*) '-------------------------------'
  write(*,*) 'Insert homogeneous direction (xy,yz,zx or skip)'
  read(*,*) answer
  call To_Upper_Case(answer)
  if(answer=='SKIP') return

  n_prob = 0
  zp=0.0
  yp=0.0

  !-----------------------------!
  !   Browse through all cells  !
  !-----------------------------!
  do c=1,NC

    ! Try to find the cell among the probes
    do p=1,n_prob
      if(answer=='YZ') then
        if( Approx(yc(c), yp(p)) .and. &  
            Approx(zc(c), zp(p)) ) go to 1
      else if(answer=='ZX') then
        if( Approx(xc(c), yp(p)) .and. &  
            Approx(zc(c), zp(p)) ) go to 1
      else if(answer=='XY') then
        if( Approx(xc(c), yp(p)) .and. &  
            Approx(yc(c), zp(p)) ) go to 1
      end if
    end do 

    ! Couldn't find a cell among the probes, add a new one
    n_prob = n_prob+1
    if(answer=='YZ') then
      yp(n_prob)=yc(c)
      zp(n_prob)=zc(c)
    else if(answer=='ZX') then
      yp(n_prob)=xc(c)
      zp(n_prob)=zc(c)
    else if(answer=='XY') then
      yp(n_prob)=xc(c)
      zp(n_prob)=yc(c)
    end if 

    if(n_prob == 20000) then
      write(*,*) 'Probe 2D: Not a 2D problem.'
      return
    end if
1 end do

  !--------------------!
  !   Create 2D file   !
  !--------------------!
  name_prob = name
  name_prob(len_trim(name)+1:len_trim(name)+3) = '.2D'
  write(6, *) 'Now creating the file:', name_prob
  open(9, file=name_prob)

  ! Write the number of probes 
  write(9,'(I8)') n_prob

  ! Write the probe coordinates out
  do p=1,n_prob
    write(9,'(I8,1PE17.8,1PE17.8)') p, yp(p), zp(p)
  end do

  close(9)

  end subroutine Probe_2D
