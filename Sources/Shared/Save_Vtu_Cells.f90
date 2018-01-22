!==============================================================================!
  subroutine Save_Vtu_Cells(grid, sub, n_nodes_sub, n_cells_sub)
!------------------------------------------------------------------------------!
! Writes: name.vtu, name.faces.vtu, name.shadow.vtu                            !
!----------------------------------[Modules]-----------------------------------!
  use all_mod, only: material
  use gen_mod, only: NewN, NewC
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, n_nodes_sub, n_cells_sub
!-----------------------------------[Locals]-----------------------------------!
  integer            :: c, n, offset
  character(len=80)  :: name_out
  integer, parameter :: VTK_TETRA      = 10
  integer, parameter :: VTK_HEXAHEDRON = 12
  integer, parameter :: VTK_WEDGE      = 13
  integer, parameter :: VTK_PYRAMID    = 14
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .vtu file   !
  !                      !
  !----------------------!
  call Name_File(sub, name_out, '.vtu')
  open(9, file=name_out)
  print *, '# Creating the file:', trim(name_out)

  !-----------!
  !   Start   !
  !-----------!
  write(9,'(a)') '<?xml version="1.0"?>'
  write(9,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                 'byte_order="LittleEndian">'
  write(9,'(a)') '<UnstructuredGrid>'
  write(9,'(a,i9)',advance='no') '<Piece NumberOfPoints="', n_nodes_sub
  write(9,'(a,i9,a)')            '" NumberOfCells ="', n_cells_sub, '">'

  !-----------!
  !   Nodes   !
  !-----------!
  write(9,'(a)') '<Points>'
  write(9,'(a)') '<DataArray type="Float32" NumberOfComponents="3" ' //  &
                 'format="ascii">'
  do n = 1, grid % n_nodes
    if(NewN(n) /= 0) write(9, '(1PE15.7,1PE15.7,1PE15.7)')                &
                                grid % xn(n), grid % yn(n), grid % zn(n)
  end do
  write(9,'(a)') '</DataArray>'
  write(9,'(a)') '</Points>'

  !-----------!
  !   Cells   !
  !-----------!
  write(9,'(a)') '<Cells>'

  ! First write all cells' nodes
  write(9,'(a)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
  do c = 1, grid % n_cells
    if(NewC(c) /= 0) then
      if(grid % cells_n_nodes(c) == 8) then
        write(9,'(8I9)')                                          &
          NewN(grid % cells_n(1,c))-1, NewN(grid % cells_n(2,c))-1,   &
          NewN(grid % cells_n(4,c))-1, NewN(grid % cells_n(3,c))-1,   &
          NewN(grid % cells_n(5,c))-1, NewN(grid % cells_n(6,c))-1,   &
          NewN(grid % cells_n(8,c))-1, NewN(grid % cells_n(7,c))-1
      else if(grid % cells_n_nodes(c) == 6) then
        write(9,'(6I9)')                                          &
          NewN(grid % cells_n(1,c))-1, NewN(grid % cells_n(2,c))-1,   &
          NewN(grid % cells_n(3,c))-1, NewN(grid % cells_n(4,c))-1,   &
          NewN(grid % cells_n(5,c))-1, NewN(grid % cells_n(6,c))-1
      else if(grid % cells_n_nodes(c) == 4) then
        write(9,'(4I9)')                                          &
          NewN(grid % cells_n(1,c))-1, NewN(grid % cells_n(2,c))-1,   &
          NewN(grid % cells_n(3,c))-1, NewN(grid % cells_n(4,c))-1
      else if(grid % cells_n_nodes(c) == 5) then
        write(9,'(5I9)')                                          &
          NewN(grid % cells_n(5,c))-1, NewN(grid % cells_n(1,c))-1,   &
          NewN(grid % cells_n(2,c))-1, NewN(grid % cells_n(4,c))-1,   &
          NewN(grid % cells_n(3,c))-1
      else
        print *, '# Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), ' nodes.'
        print *, '# Exiting'
        stop 
      end if
    end if
  end do  
  write(9,'(a)') '</DataArray>'

  ! Now write all cells' offsets
  write(9,'(a)') '<DataArray type="Int32" Name="offsets" format="ascii">'
  offset = 0
  do c = 1, grid % n_cells
    if(NewC(c) /= 0) then
      offset = offset + grid % cells_n_nodes(c)
      write(9,'(i9)') offset
    end if
  end do
  write(9,'(a)') '</DataArray>'
 
  ! Now write all cells' types
  write(9,'(a)') '<DataArray type="UInt8" Name="types" format="ascii">'
  do c = 1, grid % n_cells
    if(NewC(c) /= 0) then
      if(grid % cells_n_nodes(c) == 4) write(9,'(i9)') VTK_TETRA
      if(grid % cells_n_nodes(c) == 8) write(9,'(i9)') VTK_HEXAHEDRON
      if(grid % cells_n_nodes(c) == 6) write(9,'(i9)') VTK_WEDGE
      if(grid % cells_n_nodes(c) == 5) write(9,'(i9)') VTK_PYRAMID
    end if
  end do
  write(9,'(a)') '</DataArray>'
  write(9,'(a)') '</Cells>'
 
  !---------------!
  !   Cell data   !
  !---------------!
  write(9,'(a)') '<CellData Scalars="scalars" vectors="velocity">'

  ! Materials
  write(9,'(a)') '<DataArray type="UInt8" Name="materials" format="ascii">'
  do c = 1, grid % n_cells
    if(NewC(c) /= 0) then
      write(9,'(i9)') material(c)
    end if
  end do

  write(9,'(a)') '</DataArray>'

  !------------!
  !   Footer   !
  !------------!
  write(9,'(a)') '</CellData>'
  write(9,'(a)') '</Piece>'
  write(9,'(a)') '</UnstructuredGrid>'
  write(9,'(a)') '</VTKFile>'

  close(9)

  end subroutine
