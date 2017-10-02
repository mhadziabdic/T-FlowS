!======================================================================!
  subroutine Probe1D
!----------------------------------------------------------------------!
! This program finds the coordinate of cell-centers in non-homogeneous
! direction and write them in file name.1D
!----------------------------------------------------------------------!
  use all_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  logical :: isit
!------------------------------[Calling]-------------------------------! 
  interface
    logical function Approx(A,B,tol)
      real           :: A,B
      real, OPTIONAL :: tol
    end function Approx
  end interface
!-------------------------------[Locals]-------------------------------!
  integer   :: Nprob, p, c
  real      :: zp(1000)
  character :: namPro*80
  character :: answer*80
!======================================================================!

  write(*,*) '========================================'
  write(*,*) ' Looking for non-homogeneous directions '
  write(*,*) '----------------------------------------'
  write(*,*) 'Insert non-homogeneous direction (x,y,z or skip)'
  read(*,*) answer
  call touppr(answer)
  if(answer=='SKIP') return

  NProb = 0
  zp=0.0

  do c=1,NC

!---- try to find the cell among the probes
    do p=1,Nprob
      if(answer == 'X') then
        if( Approx(xc(c), zp(p),1.0e-9)) go to 1
      else if(answer == 'Y') then
        if( Approx(yc(c), zp(p),1.0e-9)) go to 1
      else if(answer == 'Z') then
        if( Approx(zc(c), zp(p),1.0e-9)) go to 1
      end if
    end do 

!---- couldn't find a cell among the probes, add a new one
    Nprob = Nprob+1
    if(answer=='X') zp(Nprob)=xc(c)
    if(answer=='Y') zp(Nprob)=yc(c)
    if(answer=='Z') zp(Nprob)=zc(c)

    if(Nprob == 1000) then
      write(*,*) 'Probe 1D: Not a 1D (channel flow) problem.'
      isit = .false.
      return
    end if
1 end do

  isit = .true.

!<<<<<<<<<<<<<<<<<<<<<<<<!
!     create 1D file     !
!<<<<<<<<<<<<<<<<<<<<<<<<!
  namPro = name
  namPro(len_trim(name)+1:len_trim(name)+4) = '.1Dc'
  write(6, *) 'Now creating the file:', namPro
  open(9, FILE=namPro)

!---- write the number of probes 
  write(9,'(I8)') Nprob

!---- write the probe coordinates out
  do p=1,Nprob
    write(9,'(I8,1PE17.8)') p, zp(p)
  end do

  close(9)

  end subroutine Probe1D
