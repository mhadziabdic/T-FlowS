!==============================================================================!
  subroutine Cgns_Mod_Read_Section_Info
!------------------------------------------------------------------------------!
!   Read elements connection info for current sect_id
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
    print *, "#     FAILED to read section ", sect_id, " info"
    call Cg_Error_Exit_F()
  endif

  print *, "#         section idx:  ", sect_id
  print *, "#         section name: ", sect_name
  print *, "#         section type: ", ElementTypeName(cell_type)
  print *, "#         first cell:",    first_cell
  print *, "#         last cell:",     last_cell

  i = last_cell - first_cell + 1 ! cells in this sections

  ! count cells in sect_id
  if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) n_hexa = n_hexa + i
  if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) n_pyra = n_pyra + i
  if ( ElementTypeName(cell_type) .eq. 'PENTA_6') n_pris = n_pris + i
  if ( ElementTypeName(cell_type) .eq. 'TETRA_4') n_tetr = n_tetr + i
  if ( ElementTypeName(cell_type) .eq. 'QUAD_4' ) then
    n_quad = n_quad + i
    print *, "#         This section was identified as b.c."
    bc_id(sect_id) = 1
  end if
  if ( ElementTypeName(cell_type) .eq. 'TRI_3'  ) then
    n_tria = n_tria + i
    print *, "#         This section was identified as b.c."
    bc_id(sect_id) = 1
  end if

  end subroutine
