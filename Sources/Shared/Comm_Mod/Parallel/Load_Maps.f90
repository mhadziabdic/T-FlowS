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
  integer           :: c, s         
  character(len=80) :: name_in
!==============================================================================!

  !------------------------------------------------------------------------!
  !   For run with one processor, no needd to read the map, just form it   !
  !------------------------------------------------------------------------!
  if(n_proc < 2) then

    nc_s  = grid % n_cells
    nb_s  = grid % n_bnd_cells
    nf_s  = grid % n_faces
    nbf_s = 0
    nc_t  = nc_s
    nb_t  = nb_s
    nf_t  = nf_s

    allocate(cell_map    (nc_s))
    allocate(bnd_cell_map(nb_s))
    allocate(face_map    (nf_s))

    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, nc_t
      cell_map(c) = c - 1
    end do
  
    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, nb_t
      bnd_cell_map(c) = c - 1
    end do

    ! -1 is to start from zero, as needed by MPI functions
    do s = 1, nf_t
      face_map(s) = s - 1
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
    read(9, '(4i9)') nc_s, nb_s, nf_s, nbf_s

    nc_t  = nc_s
    nb_t  = nb_s
    nf_t  = nf_s
    call Comm_Mod_Global_Sum_Int(nc_t)
    call Comm_Mod_Global_Sum_Int(nb_t)
    call Comm_Mod_Global_Sum_Int(nf_t)

    allocate(cell_map    (nc_s));   cell_map     = 1       
    allocate(bnd_cell_map(nb_s));   bnd_cell_map = 1
    allocate(face_map    (nf_s));   face_map     = 1
    allocate(buf_face_map(nbf_s));  buf_face_map = 1 

    ! Read cell map
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

    ! Read face map
    do s = 1, nf_s
      read(9, '(i9)') face_map(s)
    end do

    ! Correct cell mapping to start from zero
    face_map = face_map - 1

    ! Read face map
    do s = 1, nbf_s
      read(9, '(i9)') buf_face_map(s)
    end do

    ! Correct cell mapping to start from zero
    buf_face_map = buf_face_map - 1

    close(9)

  end if

  end subroutine
