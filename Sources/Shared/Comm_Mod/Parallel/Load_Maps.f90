!==============================================================================!
  subroutine Comm_Mod_Load_Maps(grid)
!------------------------------------------------------------------------------!
!   Reads: name.map file                                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Tokenizer_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c
  character(len=80) :: name_in
!==============================================================================!

  !------------------------------------------------------------------------!
  !   For run with one processor, no needd to read the map, just form it   !
  !------------------------------------------------------------------------!
  if(n_proc < 2) then

    nc_s = grid % n_cells
    nb_s = grid % n_bnd_cells
    nc_t = nc_s
    nb_t = nb_s

    allocate(cell_map    (nc_s))
    allocate(bnd_cell_map(nb_s))

    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, nc_t
      cell_map(c) = c - 1
    end do
  
    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, nb_t
      bnd_cell_map(c) = c - 1
    end do

  !-------------------------------------------------!
  !   For parallel runs, you need to read the map   !
  !-------------------------------------------------!
  else

    call Name_File(this_proc, name_in, '.map')
    open(9, file=name_in)
    if(this_proc < 2) print *, '# Now reading the file:', name_in

    !-----------------------!
    !   Load cell mapping   !
    !-----------------------!
    read(9, '(2i9)') nc_s, nb_s

    nc_t = nc_s
    nb_t = nb_s
    call Comm_Mod_Global_Sum_Int(nc_t)
    call Comm_Mod_Global_Sum_Int(nb_t)

    allocate(bnd_cell_map(nb_s))
    allocate(cell_map    (nc_s))

    do c = 1, nc_s
      read(9, '(i9)') cell_map(c)
    end do

    ! Correct cell mapping to start from zero
    cell_map = cell_map - 1

    ! Read boundary cell map
    do c = 1, nb_s
      read(9, '(i9)') bnd_cell_map(c)
    end do
    
    ! Correct boundary cell mapping to be positive and start from zero
    bnd_cell_map = bnd_cell_map + nb_t

    close(9)

  end if

  end subroutine
