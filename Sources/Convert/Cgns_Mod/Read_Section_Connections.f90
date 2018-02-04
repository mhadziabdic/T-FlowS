!==============================================================================!
  subroutine Cgns_Mod_Read_Section_Connections(base, block, sect, grid)
!------------------------------------------------------------------------------!
!   Read elements connection info for current sect
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer*8       :: base, block, sect
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer*8         :: base_id       ! base index number
  integer*8         :: block_id      ! block index number
  integer*8         :: sect_id       ! element section index
  character(len=80) :: sect_name     ! name of the Elements_t node
  integer*8         :: cell_type     ! types of elements in the section
  integer*8         :: first_cell    ! index of first element
  integer*8         :: last_cell     ! index of last element
  integer*8         :: n_bnd         ! index of last boundary element
  integer*8         :: parent_flag
  integer*8         :: error
  integer*8              :: n_nodes, c, n, cell, dir, cnt, bc
  integer*8, allocatable :: cell_n(:,:)
  integer*8, allocatable :: face_n(:,:)
  integer*8, allocatable :: parent_data(:,:)
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  sect_id  = sect

  ! Introduce some abbraviations 
  sect_name   = cgns_base(base) % block(block) % section(sect) % name
  cell_type   = cgns_base(base) % block(block) % section(sect) % cell_type
  first_cell  = cgns_base(base) % block(block) % section(sect) % first_cell
  last_cell   = cgns_base(base) % block(block) % section(sect) % last_cell
  parent_flag = cgns_base(base) % block(block) % section(sect) % parent_flag

  ! Number of cells in this section
  cnt = last_cell - first_cell + 1 ! cells in this sections

  if(parent_flag .eq. 1) then
    allocate(parent_data(2*cnt,2))
  end if

  !--------------------------------------------------------!
  !   Consider boundary conditions defined in this block   !
  !--------------------------------------------------------!
  do bc = 1, cgns_base(base) % block(block) % n_bnd_conds
    if(sect_name .eq. cgns_base(base) % block(block) % bnd_cond(bc) % name) then

      if(verbose) then
        print *, '#         ---------------------------------'
        print *, '#         Bnd section name:  ', sect_name
        print *, '#         ---------------------------------'
        print *, '#         Bnd section index: ', sect
        print *, '#         Bnd section type:  ', ElementTypeName(cell_type)
        print *, '#         Bnd condition mark:',   &
                 cgns_base(base) % block(block) % bnd_cond(bc) % mark
        print *, '#         Number of faces:   ', cnt
      end if

      ! Count boundary cells
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') cnt_qua = cnt_qua + cnt
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) cnt_tri = cnt_tri + cnt

      ! Update numer of boundary cells in the block
      cnt_block_bnd_cells = cnt_block_bnd_cells + cnt

      ! Allocate memory
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') n_nodes = 4
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) n_nodes = 3
      allocate(face_n(1:n_nodes, cnt))

      call Cg_Elements_Read_F(file_id,       & 
                              base_id,       &
                              block_id,      &
                              sect_id,       &
                              face_n(1,1),   &
                              parent_data,   &
                              error)          
      ! Fetch the data
      do c=1,cnt  ! I have no clue why the size has to be 2*cnt
        cell = parent_data(c,1) + cnt_cells
        dir  = parent_data(c,2)
        grid % cells_bnd_mark(dir,cell) =  &
             cgns_base(base) % block(block) % bnd_cond(bc) % mark
      end do
     
      if(verbose) then
        do c = 1,min(8,cnt)
          print '(4i7)', (face_n(n,c), n = 1, n_nodes)
        end do
        do c = 1,min(8,cnt)
          print '(a2,3i7)', 'c=', c, parent_data(c,1), parent_data(c,2)
        end do
      end if

      deallocate(face_n)

    end if
  end do

  !-------------------------------------------------!
  !   Consider three-dimensional cells / sections   !
  !-------------------------------------------------!
  if ( ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PENTA_6') .or.  &
       ( ElementTypeName(cell_type) .eq. 'TETRA_4') ) then

    if(verbose) then
      print *, '#         ---------------------------------'
      print *, '#         Cell section name: ', sect_name
      print *, '#         ---------------------------------'
      print *, '#         Cell section idx:  ', sect
      print *, '#         Cell section type: ', ElementTypeName(cell_type)
      print *, '#         Number of cells:   ', cnt
    end if

    ! Count cells in sect
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) cnt_hex = cnt_hex + cnt
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) cnt_pyr = cnt_pyr + cnt
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') cnt_wed = cnt_wed + cnt
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') cnt_tet = cnt_tet + cnt

    ! Allocate memory
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) n_nodes = 8
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) n_nodes = 5
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') n_nodes = 6
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') n_nodes = 4
    allocate(cell_n(1:n_nodes, cnt))

    call Cg_Elements_Read_F(file_id,       & 
                            base_id,       &
                            block_id,      &
                            sect_id,       &
                            cell_n(1,1),   &
                            parent_data,   &
                            error)          
    
    ! Fetch received parameters
    do c = 1, cnt
      do n = 1, n_nodes
        grid % cells_n(n, c + cnt_cells) = cell_n(n,c) + cnt_nodes
      end do
    end do

    if(verbose) then
      do c = 1, min(8,cnt)
        print '(8i7)', (cell_n(n,c), n = 1, n_nodes)
      end do
    end if
 
    deallocate(cell_n)
  end if

  end subroutine
