!==============================================================================!
  subroutine Print_Grid_Statistics 
!------------------------------------------------------------------------------!
!   Prints some statistical data about the grid on the standard output.        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
!------------------------------------------------------------------------------! 
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k, numb, nonz, stencw
!==============================================================================!

  write(*,*) '==============='
  write(*,*) 'Grid statistics'
  write(*,*) '==============='
  write(*,*) '  number of nodes         :', NN
  write(*,*) '  number of cells         :', NC
  write(*,*) '  number of sides         :', NS
  write(*,*) '  number of boundary cells:', NbC
  write(*,*) '----------------------------------'

  ! Find the number of non zero entries
  nonz=0
  do i = 1,NC
    stencw=1            ! it used to be zero
    do j=1,24
      if( CellC(i,j) > 0 ) stencw=stencw + 1
    end do
    nonz = nonz + stencw
  end do

  write(*,*) '  number of non zero matrix entries:', nonz
  write(*,*) '  average stencil size:', real(nonz)/real(NC)
  write(*,*) '  max number of nodes and cells:',   MAXN
  write(*,*) '  max number of boundary cells:',    MAXB
  write(*,*) '----------------------------------'

  ! Neighbours
  do j=1,24
    numb=0
    do i=1,NC
      stencw=0
      do k=1,24
        if( CellC(i,k)  > 0 ) stencw=stencw+1
      end do
      if(stencw  ==  j) numb=numb+1
    end do
    if(numb /= 0) then
    write(*,*) '  number of cells with ',j, ' neighbours: ',numb
    endif
  end do

  ! Twins
  do j=1,8
    numb=0
    do i=1,NN
      if(TwinN(i,0)  ==  j) numb=numb+1
    end do
    if(numb /= 0) then
    write(*,*) '  number of nodes with ',j, ' twins: ',numb
    endif 
  end do

  end subroutine Print_Grid_Statistics
