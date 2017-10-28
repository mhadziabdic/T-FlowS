!==============================================================================!
  subroutine Print_Grid_Statistics 
!------------------------------------------------------------------------------!
!   Prints some statistical data about the grid on the standard output.        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k, numb, nonz, stencw
!==============================================================================!

  write(*,*) '#==============================================================='
  write(*,*) '# Grid statistics'
  write(*,*) '#---------------------------------------------------------------'
  write(*,*) '# Number of nodes         :', NN
  write(*,*) '# Number of cells         :', NC
  write(*,*) '# Number of sides         :', NS
  write(*,*) '# Number of boundary cells:', NbC
  write(*,*) '#---------------------------------------------------------------'

  ! Find the number of non zero entries
  nonz=0
  do i = 1,NC
    stencw=1            ! it used to be zero
    do j=1,24
      if( grid % cells(i) % c(j) > 0 ) stencw=stencw + 1
    end do
    nonz = nonz + stencw
  end do

  write(*,*) '# Number of non zero matrix entries:', nonz
  write(*,*) '# Average stencil size:', real(nonz)/real(NC)
  write(*,*) '# Max number of nodes and cells:',   grid % max_n_nodes
  write(*,*) '# Max number of boundary cells:',    grid % max_n_boundary_cells
  write(*,*) '#---------------------------------------------------------------'

  ! Neighbours
  do j=1,24
    numb=0
    do i=1,NC
      stencw=0
      do k=1,24
        if( grid % cells(i) % c(k)  > 0 ) stencw=stencw+1
      end do
      if(stencw  ==  j) numb=numb+1
    end do
    if(numb /= 0) then
      write(*,*) '# Number of cells with ',j, ' neighbours: ',numb
    endif
  end do

  ! Twins
  do j=1,8
    numb=0
    do i=1,NN
      if(TwinN(i,0)  ==  j) numb=numb+1
    end do
    if(numb /= 0) then
      write(*,*) '# Number of nodes with ',j, ' twins: ',numb
    endif 
  end do

  write(*,*) '#---------------------------------------------------------------'

  end subroutine Print_Grid_Statistics
