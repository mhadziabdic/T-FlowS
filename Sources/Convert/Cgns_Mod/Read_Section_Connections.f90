!==============================================================================!
  subroutine Cgns_Mod_Read_Section_Connections
!------------------------------------------------------------------------------!
!   Read elements connection for current sect_id
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Get info for an element section
  ! Recieves sect_name, first_cell: last_cell, n_bnd and cell_type
  call Cg_Section_Read_F(file_id,      & ! cgns file index number
                         base_id,      & ! base index number
                         zone_id,      & ! zone index number
                         sect_id,      & ! element section index
                         sect_name,    & ! name of the Elements_t node
                         cell_type,    & ! type of element
                         first_cell,   & ! index of first element
                         last_cell,    & ! index of last element
                         n_bnd,        & ! index of last boundary element
                         iparent_flag, & ! if the parent data are defined
                         ier)            ! error status

  if (ier.ne.0) then
    print *, "# Failed to read info for section ", sect_id
    call Cg_Error_Exit_F()
  endif

  i = first_cell
  j = last_cell

  ! Read cells
  if (ElementTypeName(cell_type) .eq. 'HEXA_8')  allocate(buffer_r2(1:8, i:j))
  if (ElementTypeName(cell_type) .eq. 'PYRA_5')  allocate(buffer_r2(1:5, i:j))
  if (ElementTypeName(cell_type) .eq. 'PENTA_6') allocate(buffer_r2(1:6, i:j))
  if (ElementTypeName(cell_type) .eq. 'TETRA_4') allocate(buffer_r2(1:4, i:j))
  if (ElementTypeName(cell_type) .eq. 'QUAD_4')  allocate(buffer_r2(1:4, i:j))
  if (ElementTypeName(cell_type) .eq. 'TRI_3')   allocate(buffer_r2(1:3, i:j))

  ! Read HEXA_8/PYRA_5/PENTA_6/TETRA_4/QUAD_4/TRI_3 elements
  call Cg_Elements_Read_F( file_id,      & ! cgns file index number
                           base_id,      & ! base index number
                           zone_id,      & ! zone index number
                           sect_id,      & ! element section index
                           buffer_r2,    & ! element connectivity data
                           iparent_data, & ! for boundary or interface
                           ier)            ! error status

  if (ier.ne.0) then
    print *, "# Failed to read elemets in section ", sect_id
    call Cg_Error_Exit_F()
  endif

  if      (ElementTypeName(cell_type) .eq. 'HEXA_8') then
    j = last_hexa + j - i + 1
    i = last_hexa + 1
    hexa_connections(1:8, i:j) = buffer_r2(1:8, first_cell:last_cell) + last_hexa
    last_hexa = j
  else if (ElementTypeName(cell_type) .eq. 'PYRA_5') then
    j = last_pyra + j - i + 1
    i = last_pyra + 1
    pyra_connections(1:5, i:j) = buffer_r2(1:5, first_cell:last_cell) + last_pyra
    last_pyra = j
  else if (ElementTypeName(cell_type) .eq. 'PENTA_6') then
    j = last_pris + j - i + 1
    i = last_pris + 1
    pris_connections(1:6, i:j) = buffer_r2(1:6, first_cell:last_cell) + last_pris
    last_pris = j
  else if (ElementTypeName(cell_type) .eq. 'TETRA_4') then
    j = last_tetr + j - i + 1
    i = last_tetr + 1
    tetr_connections(1:4, i:j) = buffer_r2(1:4, first_cell:last_cell) + last_tetr
    last_tetr = j
  else if (ElementTypeName(cell_type) .eq. 'QUAD_4') then
    j = last_quad + j - i + 1
    i = last_quad + 1
    quad_connections(1:4, i:j) = buffer_r2(1:4, first_cell:last_cell) + last_quad
    last_quad = j
  else if (ElementTypeName(cell_type) .eq. 'TRI_3') then
    j = last_tria + j - i + 1
    i = last_tria + 1
    tria_connections(1:3, i:j) = buffer_r2(1:3, first_cell:last_cell) + last_tria
    last_tria = j
  end if

  !do c = first_cell, last_cell
  !  print "(A,8I9)", "buffer=", buffer_r2(:,c)
  !end do

  print *, i, ":", j, " copyed from", first_cell, ":", last_cell," type = ",ElementTypeName(cell_type)

  deallocate(buffer_r2)

  end subroutine
