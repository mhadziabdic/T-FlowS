!==============================================================================!
  program Divisor
!------------------------------------------------------------------------------!
!   Divides the domain in equaly balanced subdomains.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use Tokenizer_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type)  :: grid         ! grid to be divided
  integer          :: n_sub_tot    ! total number of subdomains
  integer          :: n_div        ! total number of divisions
  integer          :: chunks(128)
  integer          :: i, j, c
  character(len=8) :: answer
  real             :: start, finish
!==============================================================================!

  call cpu_time(start)

  call Logo_Div

  write(*,'(A41)') '# Input problem name: (without extension)'
  call Tokenizer_Mod_Read_Line(5)  
  read(line % tokens(1), *)  name

  ! Load the finite volume grid
  call Load_Cns           (grid, 0)
  call Allocate_Additional(grid)
  call Load_Geo           (grid, 0)
  call Load_Gmv_Faces     (grid)

  ! Initialize processor numbers
  do c = 1, grid % n_cells
    proces(c)=1
  end do

  ! Sort the cells
  print *, '# Sorting the cells'
  do i = 1, grid % n_cells
    ix(i) = i
    criter(i) = grid % xc(i) + 0.01 * grid % yc(i) + 0.0001 * grid % zc(i)
  end do
  call Sort_Real_By_Index(criter(1),ix(1),grid % n_cells,2)
  do i = 1, grid % n_cells
    iy(i) = i
    criter(i) = grid % yc(i) + 0.01 * grid % zc(i) + 0.0001 * grid % xc(i)
  end do
  call Sort_Real_By_Index(criter(1),iy(1),grid % n_cells,2)
  do i = 1, grid % n_cells
    iz(i) = i
    criter(i) = grid % zc(i) + 0.01 * grid % xc(i) + 0.0001 * grid % yc(i)
  end do
  call Sort_Real_By_Index(criter(1),iz(1),grid % n_cells,2)
  print *, '# Finished sorting'

  call Load_Geo(grid, 0)

  print *, '# Number of subdomains:'
  read(*,*)  n_sub_tot
  n_sub = n_sub_tot      ! Needed for EpsPar

  allocate (subNC(n_sub))

  print *, '#==============================='
  print *, '# Algorythm for decomposition:'
  print *, '# COO -> Coordinate multisection' 
  print *, '# INE -> Inertial multisection' 
  print *, '#-------------------------------'
  read(*,*) answer
  call To_Upper_Case(answer)
  if(answer == 'COO') then
    division_algorithm = COORDINATE
  else if(answer == 'INE') then
    division_algorithm = INERTIAL
  else
    print *, '# Error in input. Exiting!' 
    stop
  end if

  n_sub=1
  subNC(n_sub) = grid % n_cells

  !----------------------------------!
  !   Find the number of divisions   !
  !----------------------------------!
  call Factorize_Number(n_div, n_sub_tot, chunks)
  print *, '# Number of divisions:', n_div

  !-------------------------------! 
  !   Through all the divisions   !
  !-------------------------------!
  do i = 1, n_div

    do j = 1, n_sub
      print *, '# Dividing', j, ' into', chunks(i), ' chunks'
      call Split_Subdomain(grid, j, chunks(i))
    end do

  end do     

  do j = 1, n_sub
    print *, '# Processor:', j, ' cells:', subNC(j)
  end do

  call Create_Buffers_And_Save(grid)

  call Save_Com
  call Save_Scripts
  
  call cpu_time(finish)
  print '("# Time = ",f14.3," seconds.")',finish-start

  end program
