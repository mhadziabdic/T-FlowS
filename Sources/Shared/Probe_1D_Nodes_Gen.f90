!==============================================================================!
  subroutine Probe_1D_Nodes_Gen(grid)
!------------------------------------------------------------------------------!
!   This subroutine finds the coordinate of nodes in non-homogeneous           !
!   direction and write them in file name.1D                                   !
!------------------------------------------------------------------------------!
  use all_mod
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------! 
  include "Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n_prob, p, n
  real              :: n_p(10000)
  character(len=80) :: name_prob
  character(len=80) :: answer
  logical           :: isit
!==============================================================================!

  write(*,*) '#==========================================='
  write(*,*) '# Creating 1D file with the node '
  write(*,*) '# coordinates in non-homogeneous directions '
  write(*,*) '#-------------------------------------------'
  write(*,*) '# Insert non-homogeneous direction '
  write(*,*) '# (x, y, z, rx, ry, rz or skip)'
  write(*,*) '# -------------------------------------------'
  read(*,*) answer
  call To_Upper_Case(answer)
  if(answer=='SKIP') return
 
  n_prob = 0
  n_p    = 0.0

  !-----------------------------!
  !   Browse through all nodes  !
  !-----------------------------!
  do n = 1, grid % n_nodes 

    ! Try to find the cell among the probes
    do p=1,n_prob
      if(answer == 'X') then
        if( Approx(grid % xn(n), n_p(p)) ) go to 1
      else if(answer == 'Y') then
        if( Approx(grid % yn(n), n_p(p)) ) go to 1
      else if(answer == 'Z') then
        if( Approx(grid % zn(n), n_p(p)) ) go to 1
      else if(answer == 'RX') then
        if( Approx( (grid % zn(n)**2 +   &
                     grid % yn(n)**2)**.5, n_p(p)) ) go to 1
      else if(answer == 'RY') then
        if( Approx( (grid % xn(n)**2 +   &
                     grid % zn(n)**2)**.5, n_p(p)) ) go to 1
      else if(answer == 'RZ') then
        if( Approx( (grid % xn(n)**2 +   &
                     grid % yn(n)**2)**.5, n_p(p)) ) go to 1
      end if
    end do 
  
    ! Couldn't find a cell among the probes, add a new one
    n_prob = n_prob+1
    if(answer=='X') n_p(n_prob)= grid % xn(n)
    if(answer=='Y') n_p(n_prob)= grid % yn(n)
    if(answer=='Z') n_p(n_prob)= grid % zn(n)

    if(answer=='RX') n_p(n_prob)= (grid % zn(n)**2 +      &
                                   grid % yn(n)**2)**0.5
    if(answer=='RY') n_p(n_prob)= (grid % xn(n)**2 +      &
                                   grid % zn(n)**2)**0.5
    if(answer=='RZ') n_p(n_prob)= (grid % xn(n)**2 +      &
                                   grid % yn(n)**2)**0.5

    if(n_prob == 10000) then
      write(*,*) '# Probe 1D: Not a 1D (channel flow) problem.'
      isit = .false.
      return
    end if
1  end do

  isit = .true.

  !--------------------!
  !   Create 1D file   !
  !--------------------!
  name_prob = name
  name_prob(len_trim(name)+1:len_trim(name)+3) = '.1D'
  write(6, *) 'Now creating the file:', name_prob
  open(9, file=name_prob)
  ! Write the number of probes 
  write(9,'(I8)') n_prob

!  call SSORT (n_p, 1, n_prob*3, 1)
!  call  RISort(n_p, n_prob, n_prob*6,2)
  call SORT2(n_p, n_prob*2, n_prob)

  ! Write the probe coordinates out
  do p=1, n_prob
    write(9,'(I8,1E17.8)') p, n_p(p)
  end do

  close(9)

  end subroutine
