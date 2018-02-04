!==============================================================================!
  subroutine Cgns_Mod_Read_Section_Connections(base, block, sect)
!------------------------------------------------------------------------------!
!   Read elements connection info for current sect
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer*8 :: base, block, sect
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
  integer*8         :: cnt, bc
  integer*8         :: parent_data = 0
  integer*8              :: n_nodes, c, n
  integer*8, allocatable :: cell_n(:,:)
  integer*8, allocatable :: face_n(:,:)
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  sect_id  = sect

  ! Introduce some abbraviations 
  sect_name  = cgns_base(base) % block(block) % section(sect) % name
  cell_type  = cgns_base(base) % block(block) % section(sect) % cell_type
  first_cell = cgns_base(base) % block(block) % section(sect) % first_cell
  last_cell  = cgns_base(base) % block(block) % section(sect) % last_cell

  ! Number of cells in this section
  cnt = last_cell - first_cell + 1 ! cells in this sections

  ! Consider boundary conditions defined in this block
  do bc = 1, cgns_base(base) % block(block) % n_bnd_conds
    if(sect_name .eq. cgns_base(base) % block(block) % bnd_cond(bc) % name) then

      print *, '#         ---------------------------------'
      print *, '#         Bnd section name:  ', sect_name
      print *, '#         ---------------------------------'
      print *, '#         Bnd section index: ', sect
      print *, '#         Bnd section type:  ', ElementTypeName(cell_type)
      print *, '#         Number of faces:   ', cnt

      ! Count boundary cells
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') cnt_qua = cnt_qua + cnt
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) cnt_tri = cnt_tri + cnt

      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') n_nodes = 4
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) n_nodes = 3

      ! Allocate memory
      allocate(face_n(1:n_nodes, cnt))
      call Cg_Elements_Read_F(file_id,       & 
                              base_id,       &
                              block_id,      &
                              sect_id,       &
                              face_n(1,1),   &
                              parent_data,   &
                              error)          
      do c=1,10
        print '(4i7)', (face_n(n,c), n = 1, n_nodes)
      end do
      deallocate(face_n)

    end if
  end do

  ! Consider only three-dimensional cells / sections
  if ( ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PENTA_6') .or.  &
       ( ElementTypeName(cell_type) .eq. 'TETRA_4') ) then

    print *, '#         ---------------------------------'
    print *, '#         Cell section name: ', sect_name
    print *, '#         ---------------------------------'
    print *, '#         Cell section idx:  ', sect
    print *, '#         Cell section type: ', ElementTypeName(cell_type)
    print *, '#         Number of cells:   ', cnt

    ! Count cells in sect
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) cnt_hex = cnt_hex + cnt
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) cnt_pyr = cnt_pyr + cnt
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') cnt_wed = cnt_wed + cnt
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') cnt_tet = cnt_tet + cnt

    ! Count cells in sect
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) n_nodes = 8
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) n_nodes = 5
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') n_nodes = 6
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') n_nodes = 4
 
    ! Allocate memory
    allocate(cell_n(1:n_nodes, cnt))
    call Cg_Elements_Read_F(file_id,       & 
                            base_id,       &
                            block_id,      &
                            sect_id,       &
                            cell_n(1,1),   &
                            parent_data,   &
                            error)          
    do c=1,10
      print '(8i7)', (cell_n(n,c), n = 1, n_nodes)
    end do
    deallocate(cell_n)
  end if

  end subroutine
