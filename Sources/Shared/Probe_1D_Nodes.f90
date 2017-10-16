!==============================================================================!
  subroutine Probe_1D_Nodes
!------------------------------------------------------------------------------!
!   This subroutine finds the coordinate of cell-centers in non-homogeneous    !
!   direction and write them in file called "name.1D"                          !
!------------------------------------------------------------------------------!
  use all_mod
  use gen_mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical :: isit
!----------------------------------[Calling]-----------------------------------! 
  include "Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n_prob, p, c, n
  real              :: n_p(10000)
  character(len=80) :: name_prob
  character(len=80) :: answer
!==============================================================================!

  write(*,*) '==========================================='
  write(*,*) ' Creating 1D file with the node '
  write(*,*) ' coordinates in non-homogeneous directions '
  write(*,*) '-------------------------------------------'
  write(*,*) 'Insert non-homogeneous direction '// &
             '(x_node, y_node, z_node, Rx, Ry, Rz or skip)'
  read(*,*) answer
  call To_Upper_Case(answer)
  if(answer=='SKIP') return
 
  n_prob = 0
  n_p   = 0.0

  !-----------------------------!
  !   Browse through all cells  !
  !-----------------------------!
  do c=-NbC, NC
    do n=1,CellN(c,0)

      ! Try to find the cell among the probes
      do p=1,n_prob
          if(answer == 'X') then
            if( Approx(x_node(CellN(c,n)), n_p(p)) ) go to 1
          else if(answer == 'Y') then
            if( Approx(y_node(CellN(c,n)), n_p(p)) ) go to 1
          else if(answer == 'Z') then
            if( Approx(z_node(CellN(c,n)), n_p(p)) ) go to 1
          else if(answer == 'RX') then
            if( Approx( (z_node(CellN(c,n))**2.0 +   &
                         y_node(CellN(c,n))**2.0)**0.5, n_p(p)) ) go to 1
          else if(answer == 'RY') then
            if( Approx( (x_node(CellN(c,n))**2.0 +   &
                         z_node(CellN(c,n))**2.0)**0.5, n_p(p)) ) go to 1
          else if(answer == 'RZ') then
            if( Approx( (x_node(CellN(c,n))**2.0 +   &
                         y_node(CellN(c,n))**2.0)**0.5, n_p(p)) ) go to 1
          end if
      end do 
  
      ! Couldn't find a cell among the probes, add a new one
      n_prob = n_prob+1
      if(answer=='X') n_p(n_prob)= x_node(CellN(c,n))
      if(answer=='Y') n_p(n_prob)= y_node(CellN(c,n))
      if(answer=='Z') n_p(n_prob)= z_node(CellN(c,n))

      if(answer=='RX') n_p(n_prob)= (z_node(CellN(c,n))**2.0 +  &
                                     y_node(CellN(c,n))**2.0)**0.5
      if(answer=='RY') n_p(n_prob)= (x_node(CellN(c,n))**2.0 +  &
                                     z_node(CellN(c,n))**2.0)**0.5
      if(answer=='RZ') n_p(n_prob)= (x_node(CellN(c,n))**2.0 +  &
                                     y_node(CellN(c,n))**2.0)**0.5

      if(n_prob == 10000) then
        write(*,*) 'Probe 1D: Not a 1D (channel flow) problem.'
        isit = .false.
        return
      end if
    end do
1 end do

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

  end subroutine Probe_1D_Nodes
