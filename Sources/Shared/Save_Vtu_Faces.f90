!==============================================================================!
  subroutine Save_Vtu_Faces(grid)
!------------------------------------------------------------------------------!
! Writes .faces.vtu file.                                                      !
!----------------------------------[Modules]-----------------------------------!
  use all_mod, only: bcmark
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, NNsub
!-----------------------------------[Locals]-----------------------------------!
  integer             :: c1, c2, n, s, offset
  character(len=80)  :: name_out
  integer, parameter :: VTK_TRIANGLE = 5
  integer, parameter :: VTK_QUAD     = 9
!==============================================================================!

  !-----------------------------------------!
  !                                         !
  !   Create boundary condition .vtu file   !
  !                                         !
  !-----------------------------------------!
  call Name_File(sub, name_out, '.faces.vtu')
  open(9, file=name_out)
  print *, '# Creating the file:', trim(name_out)

  !-----------!
  !   Start   !
  !-----------!
  write(9,'(a)') '<?xml version="1.0"?>'
  write(9,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                 'byte_order="LittleEndian">'
  write(9,'(a)') '<UnstructuredGrid>'
  write(9,'(a,i9)',advance='no') '<Piece NumberOfPoints="', grid % n_nodes
  write(9,'(a,i9,a)')            '" NumberOfCells ="', grid % n_faces, '">'

  !-----------!
  !   Nodes   !
  !-----------!
  write(9,'(a)') '<Points>'
  write(9,'(a)') '<DataArray type="Float32" NumberOfComponents="3" ' //  &
                 'format="ascii">'
  do n = 1, grid % n_nodes
    write(9, '(1PE15.7,1PE15.7,1PE15.7)')                &
               grid % xn(n), grid % yn(n), grid % zn(n)
  end do
  write(9,'(a)') '</DataArray>'
  write(9,'(a)') '</Points>'

  !-----------!
  !   Faces   !
  !-----------!
  write(9,'(a)') '<Cells>'

  ! First write all faces' nodes
  write(9,'(a)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
  do s = 1, grid % n_faces
    if(grid % faces_n_nodes(s) == 4) then
      write(9,'(4I9)')                                 &
        grid % faces_n(1,s)-1, grid % faces_n(2,s)-1,  &
        grid % faces_n(3,s)-1, grid % faces_n(4,s)-1
    else if(grid % faces_n_nodes(s) == 3) then
      write(9,'(3I9)')                                 &
        grid % faces_n(1,s)-1, grid % faces_n(2,s)-1,  &
        grid % faces_n(3,s)-1
    else
      print *, '# Unsupported cell type ',       &
                 grid % faces_n_nodes(s), ' nodes.'
      print *, '# Exiting'
      stop 
    end if
  end do  
  write(9,'(a)') '</DataArray>'

  ! Then write all faces' offsets
  write(9,'(a)') '<DataArray type="Int32" Name="offsets" format="ascii">'
  offset = 0
  do s = 1, grid % n_faces
    offset = offset + grid % faces_n_nodes(s)
    write(9,'(i9)') offset
  end do
  write(9,'(a)') '</DataArray>'
 
  ! Now write all cells' types
  write(9,'(a)') '<DataArray type="UInt8" Name="types" format="ascii">'
  do s = 1, grid % n_faces
    if(grid % faces_n_nodes(s) == 4) write(9,'(i9)') VTK_QUAD
    if(grid % faces_n_nodes(s) == 3) write(9,'(i9)') VTK_TRIANGLE
  end do
  write(9,'(a)') '</DataArray>'
  write(9,'(a)') '</Cells>'
 
  !---------------!
  !   Cell data   !
  !---------------!
  write(9,'(a)') '<CellData Scalars="scalars" vectors="velocity">'

  ! Boundary conditions
  write(9,'(a)') '<DataArray type="UInt8" Name="boundary conditions" ' //  & 
                 'format="ascii">'
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
   
    ! If boundary 
    if( c2 < 0 ) then 
      write(9,'(i9)') bcmark(c2)

    ! If inside 
    else 
      write(9,'(i9)') 0
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
