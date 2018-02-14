!==============================================================================!
  subroutine Cgns_Mod_Write_Section_Connections_Seq(base, block, sect, grid)
!------------------------------------------------------------------------------!
!   Gets n_sects from block
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer              :: base, block, sect
  type(Grid_Type)      :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: base_id       ! base index number
  integer              :: block_id      ! block index number
  integer              :: sect_id       ! element section index
  character(len=80)    :: sect_name     ! name of the Elements_t node
  integer              :: cell_type     ! types of elements in the section
  integer              :: first_cell    ! index of first element
  integer              :: last_cell     ! index of last element
  integer              :: n_bnd         ! index of last boundary element
  integer              :: error
  integer              :: n_nodes, c, i
  integer, allocatable :: cell_n(:,:)
!==============================================================================!



  ! Set input parameters
  base_id    = base
  block_id   = block
  sect_id    = sect

  sect_name  = trim(cgns_base(base_id)%block(block_id)%section(sect_id)%name)
  cell_type  = cgns_base(base_id)%block(block_id)%section(sect_id)%cell_type 
  first_cell = cgns_base(base_id)%block(block_id)%section(sect_id)%first_cell
  last_cell  = cgns_base(base_id)%block(block_id)%section(sect_id)%last_cell 


      ! Print some info
    if(verbose ) then
      print *, '#         ---------------------------------'
      print *, '#         Cell section name: ', sect_name
      print *, '#         ---------------------------------'
      print *, '#         Cell section idx:    ', sect_id
      print *, '#         Cell section type:   ', ElementTypeName(cell_type)
    end if
    if(verbose) then
      print *, '#         Number of cells: ', last_cell - first_cell + 1
      print *, '#         First cell:', first_cell
      print *, '#         Last cell: ', last_cell
    end if

  !---------- create and fill "Hexagons" node in DB
  if ( last_cell .ne. 0 ) then

    n_bnd = 0 ! unsorted boundary elements

    ! Allocate memory
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) n_nodes = 8
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) n_nodes = 5
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') n_nodes = 6
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') n_nodes = 4
    allocate(cell_n(1:n_nodes, first_cell:last_cell), stat = error)

    if (error .ne. 0) then
       print*, '*FAILED* to allocate ', trim(sect_name)
       call cg_error_exit_f()
    endif

    ! convert T-FlowS -> CGNS [same as VTK]
    i = 1
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
      else
        print *, '# Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), ' nodes.'
        print *, '# Exiting'
        stop
      end if
    end do

    ! Write element data 
    call Cg_Section_Write_F(file_id,    &
                            base_id,    &
                            block_id,   &
                            sect_name,  &
                            cell_type,  &
                            first_cell, &
                            last_cell,  &
                            n_bnd,      &
                            cell_n,     &
                            sect_id,    &
                            error)

    if (error .ne. 0) then
       print*, '*FAILED* to write ', trim(sect_name), ' connections'
       call cg_error_exit_f()
    endif

    deallocate(cell_n)
  end if

  end subroutine