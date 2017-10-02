!======================================================================!
  subroutine Probe2D 
!----------------------------------------------------------------------!
! Finds coordinates of all the planes for the channel flow.            !
! It assumes that homogeneous directions of the flow are x and y.      !
!----------------------------------------------------------------------!
  use all_mod
!----------------------------------------------------------------------!
  implicit none
!------------------------------[Calling]-------------------------------! 
  interface
    logical function Approx(A,B,tol)
      real           :: A,B
      real, optional :: tol
    end function Approx
  end interface
!-------------------------------[Locals]-------------------------------!
  integer   :: Nprob, p, c
  real      :: yp(20000), zp(20000)
  character :: namPro*80
  character :: answer*80
!======================================================================!

  write(*,*) '==============================='
  write(*,*) ' Looking for homogeneous plane '
  write(*,*) '-------------------------------'
  write(*,*) 'Insert homogeneous direction (xy,yz,zx or skip)'
  read(*,*) answer
  call touppr(answer)
  if(answer=='SKIP') return

  NProb = 0
  zp=0.0
  yp=0.0

  do c=1,NC

!---- try to find the cell among the probes
    do p=1,Nprob
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

!---- couldn't find a cell among the probes, add a new one
    Nprob = Nprob+1
    if(answer=='YZ') then
      yp(Nprob)=yc(c)
      zp(Nprob)=zc(c)
    else if(answer=='ZX') then
      yp(Nprob)=xc(c)
      zp(Nprob)=zc(c)
    else if(answer=='XY') then
      yp(Nprob)=xc(c)
      zp(Nprob)=yc(c)
    end if 

    if(Nprob == 20000) then
      write(*,*) 'Probe 2D: Not a 2D problem.'
      return
    end if
1 end do

!<<<<<<<<<<<<<<<<<<<<<<<<!
!     create 2D file     !
!<<<<<<<<<<<<<<<<<<<<<<<<!
  namPro = name
  namPro(len_trim(name)+1:len_trim(name)+3) = '.2D'
  write(6, *) 'Now creating the file:', namPro
  open(9, FILE=namPro)

!---- write the number of probes 
  write(9,'(I8)') Nprob

!---- write the probe coordinates out
  do p=1,Nprob
    write(9,'(I8,1PE17.8,1PE17.8)') p, yp(p), zp(p)
  end do

  close(9)

  end subroutine Probe2D
