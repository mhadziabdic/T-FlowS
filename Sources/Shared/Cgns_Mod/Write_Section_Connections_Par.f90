!==============================================================================!
  subroutine Cgns_Mod_Write_Section_Connections_Par(base, block, sect, grid)
!------------------------------------------------------------------------------!
!   Gets n_sects from block
!------------------------------------------------------------------------------!
!   Array structures in current function are strictly followings:              !
!   Processor:    |        P_1        |               P_2               | ...  !
!   x,y,z:        |      (1 : NN_1)   |       NN_1 + 1 : NN_1 + NN_2    | ...  ! 
!   Connections:  |   (8, 1 : NC_1)   |   (8, NC_1 + 1 : NC_1 + NC_2)   | ...  ! 
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use par_mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer         :: base, block, sect
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base_id       ! base index number
  integer              :: block_id      ! block index number
  integer              :: sect_id       ! element section index
  character(len=80)    :: sect_name     ! name of the Elements_t node
  integer              :: cell_type     ! types of elements in the section
  integer              :: n_bnd         ! index of last boundary element
  integer              :: error
  integer              :: n_nodes, c, cnt, i, j
  integer, allocatable :: cell_n(:,:)
  integer              :: first_cell    ! look at array structure at the header
  integer              :: last_cell     ! look at array structure at the header
  integer              :: first_node    ! look at array structure at the header
!==============================================================================!

  ! Set input parameters
  base_id    = base
  block_id   = block
  sect_id    = sect
  sect_name  = trim(cgns_base(base_id)%block(block_id)%section(sect_id)%name)
  cell_type  = cgns_base(base_id)%block(block_id)%section(sect_id)%cell_type 
  i          = cgns_base(base_id)%block(block_id)%section(sect_id)%first_cell
  j          = cgns_base(base_id)%block(block_id)%section(sect_id)%last_cell
  ! cells of cell_type on this_proc
  cnt        = j - i + 1

  if ( cnt .ne. 0 ) then
  !----------------------------------------!
  !   Create empty elem node in DB block   !
  !----------------------------------------!

    ! total cells of cell_type
    i = 1
    c = cnt
    call wait
    call IglSum(c)
    n_bnd = 0 ! unsorted boundary elements

    !if (this_proc == 1 ) then ! this if is needed since sect_name ='' for n_nroc>1
      !---------- Create empty elem node in DB
      call Cgp_Section_Write_F( &
        file_id,                &
        base_id,                &
        block_id,               &
        sect_name,              &
        cell_type,              &
        i,                      &
        c,                      &
        n_bnd,                  &
        sect_id,                &
        error)
      if (error .ne. 0) then
         print*, '*FAILED* to create empty ', trim(sect_name)
         call Cgp_Error_Exit_F()
      endif
    !end if
    call wait

  !--------------------------------------!
  !   Fill empty elem node in DB block   !
  !--------------------------------------!
    i = cnt ! cnt_hex/pyr/wed/tet
    call Cgns_Mod_Get_Arrays_Dimensions_Par(first_cell, i)

    last_cell = first_cell + cnt - 1

    ! Allocate memory
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) n_nodes = 8
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) n_nodes = 5
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') n_nodes = 6
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') n_nodes = 4
    allocate(cell_n(1:n_nodes, first_cell:last_cell), stat = error)

    if (error .ne. 0) then
       print*, '*FAILED* to allocate ', trim(sect_name)
       call Cgp_Error_Exit_F()
    endif

    ! convert T-FlowS -> CGNS [same as VTK]
    i = first_cell
    do c = 1, grid % n_cells
      if (grid % cells_n_nodes(c) .eq. 8) then ! hex
        cell_n (1, i) = grid % cells_n(1, c)
        cell_n (2, i) = grid % cells_n(2, c)
        cell_n (3, i) = grid % cells_n(4, c)
        cell_n (4, i) = grid % cells_n(3, c)
        cell_n (5, i) = grid % cells_n(5, c)
        cell_n (6, i) = grid % cells_n(6, c)
        cell_n (7, i) = grid % cells_n(8, c)
        cell_n (8, i) = grid % cells_n(7, c)
        i = i + 1
      elseif (grid % cells_n_nodes(c) .eq. 6) then ! wedge
        cell_n (1:6, i) = grid % cells_n(1:6, c)
        i = i + 1
      elseif (grid % cells_n_nodes(c) .eq. 5) then ! pyramid
        cell_n (1, i) = grid % cells_n(5, c)
        cell_n (2, i) = grid % cells_n(1, c)
        cell_n (3, i) = grid % cells_n(2, c)
        cell_n (4, i) = grid % cells_n(4, c)
        cell_n (5, i) = grid % cells_n(3, c)
        i = i + 1
      elseif (grid % cells_n_nodes(c) .eq. 4) then ! tetra
        cell_n (1:4, i) = grid % cells_n(1:4, c)
        i = i + 1
      end if
    end do

    ! shift cell_n according to a number of nodes in coord array
    i = grid % n_nodes
    call Cgns_Mod_Get_Arrays_Dimensions_Par(first_node, i)
    cell_n = cell_n + first_node - 1

    !---------- Fill that node with grid coordinates 
    call Cgp_Elements_Write_Data_F( &
      file_id,                      &
      base_id,                      &
      block_id,                     &
      sect_id,                      &
      first_cell,                   &
      last_cell,                    &
      cell_n,                       &
      error)

    if (error .ne. 0) then
       print*, '*FAILED* to fill ', sect_id
       call Cgp_Error_Exit_F()
    endif

    deallocate(cell_n)

    ! Print some info
    if(verbose .and. this_proc.eq.1) then
      print *, '#         ---------------------------------'
      print *, '#         Cell section name: ', sect_name
      print *, '#         ---------------------------------'
      print *, '#         Cell section idx:    ', sect_id
      print *, '#         Cell section type:   ', ElementTypeName(cell_type)
    end if
    if(verbose) then
      print *, '#         Number of cells: ', cnt, " (P:",this_proc,")"
      print *, '#         First cell:', first_cell, " (P:",this_proc,")"
      print *, '#         Last cell: ', last_cell, " (P:",this_proc,")"
    end if

  end if

  end subroutine