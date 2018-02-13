!==============================================================================!
  subroutine Cgns_Mod_Write_Base_Info(base)
!------------------------------------------------------------------------------!
!   Reads main info from base node base_id
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: base
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id    ! base index number
  character(len=80) :: base_name  ! name of the base
  integer           :: cell_dim   ! cell dimensions (3->volume cell)
  integer           :: phys_dim   ! number of coordinates to create vector
  integer           :: error
!==============================================================================!

  ! Set input parameters
  base_id   = base
  base_name = cgns_base(base) % name    
  cell_dim  = cgns_base(base) % cell_dim
  phys_dim  = cgns_base(base) % phys_dim

  ! Read CGNS base information
  call Cg_Base_Write_F(file_id,    &
                       base_name,  &
                       cell_dim,   &
                       phys_dim,   &
                       base_id,    &
                       error)

  if (error .ne. 0) then
    print *, "# Failed to get base info"
    call Cg_Error_Exit_F()
  endif

  ! Print some info
  if(verbose) then
    print *, '#   ============================'
    print *, '#   '
    print *, '#   Base name: ',      base_name
    print *, '#   '
    print *, '#   ============================'
    print *, '#   Base name: ',      base_name
    print *, '#   Cell dimension: ', cell_dim
    print *, '#   Phys dimension: ', phys_dim
  end if

  end subroutine
