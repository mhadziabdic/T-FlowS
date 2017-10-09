!======================================================================!
  subroutine Probe1D_nodes_gen
!----------------------------------------------------------------------!
! This program finds the coordinate of nodes in non-homogeneous  
! direction and write them in file name.1D
!----------------------------------------------------------------------!
  use all_mod
  use gen_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  logical :: isit
!------------------------------[Calling]-------------------------------! 
  interface
    logical function Approx(A,B,tol)
      real           :: A,B
      real, optional :: tol
    end function Approx
  end interface
!-------------------------------[Locals]-------------------------------!
  integer   :: Nprob, p, c, n
  real      :: N_p(10000)
  character :: namPro*80
  character :: answer*80
!======================================================================!

  write(*,*) '==========================================='
  write(*,*) ' Creating 1D file with the node '
  write(*,*) ' coordinates in non-homogeneous directions '
  write(*,*) '-------------------------------------------'
  write(*,*) 'Insert non-homogeneous direction (x_node, y_node, z_node, Rx, Ry, Rz or skip)'
  read(*,*) answer
  call touppr(answer)
  if(answer=='SKIP') return
 
  NProb = 0
  N_p   = 0.0

  do n =  1, NN 
!---- try to find the cell among the probes
    do p=1,Nprob
        if(answer == 'X') then
          if( Approx(x_node(n), N_p(p)) ) go to 1
        else if(answer == 'Y') then
          if( Approx(y_node(n), N_p(p)) ) go to 1
        else if(answer == 'Z') then
          if( Approx(z_node(n), N_p(p)) ) go to 1
        else if(answer == 'RX') then
          if( Approx( (z_node(n)**2.0 + y_node(n)**2.0)**0.5, N_p(p)) ) go to 1
        else if(answer == 'RY') then
          if( Approx( (x_node(n)**2.0 + z_node(n)**2.0)**0.5, N_p(p)) ) go to 1
        else if(answer == 'RZ') then
          if( Approx( (x_node(n)**2.0 + y_node(n)**2.0)**0.5, N_p(p)) ) go to 1
        end if
    end do 
  
!---- couldn't find a cell among the probes, add a new one
    Nprob = Nprob+1
    if(answer=='X') N_p(Nprob)=x_node(n)
    if(answer=='Y') N_p(Nprob)=y_node(n)
    if(answer=='Z') N_p(Nprob)=z_node(n)

    if(answer=='RX') N_p(Nprob)=(z_node(n)**2.0 + y_node(n)**2.0)**0.5
    if(answer=='RY') N_p(Nprob)=(x_node(n)**2.0 + z_node(n)**2.0)**0.5
    if(answer=='RZ') N_p(Nprob)=(x_node(n)*x_node(n) + y_node(n)*y_node(n))**0.5

    if(Nprob == 10000) then
      write(*,*) 'Probe 1D: Not a 1D (channel flow) problem.'
      isit = .false.
      return
    end if
1  end do

  isit = .true.

!<<<<<<<<<<<<<<<<<<<<<<<<!
!     create 1D file     !
!<<<<<<<<<<<<<<<<<<<<<<<<!
  namPro = name
  namPro(len_trim(name)+1:len_trim(name)+3) = '.1D'
  write(6, *) 'Now creating the file:', namPro
  open(9, FILE=namPro)
!---- write the number of probes 
  write(9,'(I8)') Nprob

!  call SSORT (N_p, 1, Nprob*3, 1)
!  call  RISort(N_p, Nprob, Nprob*6,2)

  call SORT2(N_p, Nprob*2, Nprob)

!---- write the probe coordinates out
  do p=1, Nprob
    write(9,'(I8,1E17.8)') p, N_p(p)
  end do

  close(9)

  end subroutine Probe1D_nodes_gen
