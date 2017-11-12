!==============================================================================!
  program Divisor
!------------------------------------------------------------------------------!
!   Divides the domain in equaly balanced subdomains.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use par_mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer          :: n_sub_tot            ! total number of subdomains
  integer          :: n_div                ! total number of divisions
  integer          :: chunks(128)
  integer          :: i, j, c
  character(len=8) :: answer
  real             :: start, finish
!==============================================================================!

  call cpu_time(start)

  ! Test the precision
  open(90,FORM='unformatted',file='Divisor.real');
  write(90) 3.1451592
  close(90)
               
  call Logo

  write(*,'(A41)') '# Input problem name: (without extension)'
  call ReadC(5,inp,tn,ts,te)  
  read(inp, '(A80)')  name

  ! Load the finite volume grid
  call Load_Cns
  call Load_Gmv
  call Load_Geo
  call BCelLoa

  ! Initialize processor numbers
  do c=1,NC
    proces(c)=1
  end do

  ! Sort the cells
  write(*,*) '# Sorting the cells'
  do i=1,NC
    ix(i) = i
    criter(i) = xc(i) + 0.01 * yc(i) + 0.0001 * zc(i)
  end do
  call Sort_Real_By_Index(criter(1),ix(1),NC,2)
  do i=1,NC
    iy(i) = i
    criter(i) = yc(i) + 0.01 * zc(i) + 0.0001 * xc(i)
  end do
  call Sort_Real_By_Index(criter(1),iy(1),NC,2)
  do i=1,NC
    iz(i) = i
    criter(i) = zc(i) + 0.01 * xc(i) + 0.0001 * yc(i)
  end do
  call Sort_Real_By_Index(criter(1),iz(1),NC,2)
  write(*,*) '# Finished sorting'

  call Load_Geo

  write(*,*) '# Number of subdomains:'
  read(*,*)  n_sub_tot
  n_proc = n_sub_tot      ! Needed for EpsPar

  allocate (subNC(n_proc))
  allocate (NBBs(0:n_proc))
  allocate (NBBe(0:n_proc))

  write(*,*) '#==============================='
  write(*,*) '# Algorythm for decomposition:'
  write(*,*) '# COO -> Coordinate multisection' 
  write(*,*) '# INE -> Inertial multisection' 
  write(*,*) '#-------------------------------'
  read(*,*) answer
  call To_Upper_Case(answer)
  if(answer == 'COO') then
    division_algorithm = COORDINATE
  else if(answer == 'INE') then
    division_algorithm = INERTIAL
  else
    write(*,*) '# Error in input. Exiting!' 
    stop
  end if

  n_sub=1
  subNC(n_sub) = NC

  !----------------------------------!
  !   Find the number of divisions   !
  !----------------------------------!
  call Factorize_Number(n_div, n_sub_tot, chunks)
  write(*,*) '# Number of divisions:', n_div

  !-------------------------------! 
  !   Through all the divisions   !
  !-------------------------------!
  do i = 1, n_div

    do j=1,n_sub
      write(*,*) '# Dividing', j, ' into', chunks(i), ' chunks'
      call Split_Subdomain(j, chunks(i))
    end do

  end do     

  do j=1,n_sub
    write(*,*) '# Processor:', j, ' cells:', subNC(j)
  end do

  call Number 

  call Save_Com
  call Save_Scripts
  
  call cpu_time(finish)
  print '("# Time = ",f14.3," seconds.")',finish-start

  end program Divisor