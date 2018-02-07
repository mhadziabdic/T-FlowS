!==============================================================================!
  subroutine Save_Vtu_Results(grid, name_save)
!------------------------------------------------------------------------------!
!   Writes results in VTU file format (for VisIt and Paraview)                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod
  use all_mod
  use pro_mod
  use rans_mod
  use par_mod, only: this_proc, n_proc
  use Tokenizer_Mod
  use Grid_Mod
  use Constants_Pro_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_save
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n, offset
  character(len=80) :: name_out, store_name
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: VTK_TETRA      = 10  ! cell shapes in VTK format
  integer, parameter :: VTK_HEXAHEDRON = 12  
  integer, parameter :: VTK_WEDGE      = 13
  integer, parameter :: VTK_PYRAMID    = 14
  character(len= 0)  :: IN_0 = ''           ! indentation levels 
  character(len= 2)  :: IN_1 = '  '
  character(len= 4)  :: IN_2 = '    '
  character(len= 6)  :: IN_3 = '      '
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
  character(len=12)  :: IN_6 = '            '
  character(len=14)  :: IN_7 = '              '
  character(len=16)  :: IN_8 = '                '
!==============================================================================!

  ! Store the name
  store_name = problem_name     

  problem_name = name_save  

  call Wait 

  !----------------------!
  !                      !
  !   Create .vtu file   !
  !                      !
  !----------------------!
  call Name_File(this_proc, name_out, '.vtu')
  open(9, file=name_out)
  print *, '# Creating the file: ', trim(name_out)

  !-----------!
  !   Start   !
  !-----------!
  write(9,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(9,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                         'byte_order="LittleEndian">'
  write(9,'(a,a)') IN_1, '<UnstructuredGrid>'
  write(9,'(a,a,i0.0,a,i0.0,a)')   &
                   IN_2, '<Piece NumberOfPoints="', grid % n_nodes,      &
                              '" NumberOfCells ="', grid % n_cells, '">'

  !----------!
  !          !
  !   Grid   !
  !          !
  !----------!

  !-----------!
  !   Nodes   !
  !-----------!
  write(9,'(a,a)') IN_3, '<Points>'
  write(9,'(a,a)') IN_4, '<DataArray type="Float32" NumberOfComponents' //  &
                         '="3" format="ascii">'
  do n = 1, grid % n_nodes
    write(9, '(a,1pe15.7,1pe15.7,1pe15.7)')                &
               IN_5, grid % xn(n), grid % yn(n), grid % zn(n)
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'
  write(9,'(a,a)') IN_3, '</Points>'

  !-----------!
  !   Cells   !
  !-----------!
  write(9,'(a,a)') IN_3, '<Cells>'

  ! First write all cells' nodes
  write(9,'(a,a)') IN_4, '<DataArray type="Int32" Name="connectivity"' //  &
                         ' format="ascii">'

  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) == 8) then
      write(9,'(a,8i9)')                                &
        IN_5,                                           &
        grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
        grid % cells_n(4,c)-1, grid % cells_n(3,c)-1,   &
        grid % cells_n(5,c)-1, grid % cells_n(6,c)-1,   &
        grid % cells_n(8,c)-1, grid % cells_n(7,c)-1
    else if(grid % cells_n_nodes(c) == 6) then
      write(9,'(a,6i9)')                                &
        IN_5,                                           &
        grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
        grid % cells_n(3,c)-1, grid % cells_n(4,c)-1,   &
        grid % cells_n(5,c)-1, grid % cells_n(6,c)-1
    else if(grid % cells_n_nodes(c) == 4) then
      write(9,'(a,4i9)')                                &
        IN_5,                                           &
        grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
        grid % cells_n(3,c)-1, grid % cells_n(4,c)-1
    else if(grid % cells_n_nodes(c) == 5) then
      write(9,'(a,5i9)')                                &
        IN_5,                                           &
        grid % cells_n(5,c)-1, grid % cells_n(1,c)-1,   &
        grid % cells_n(2,c)-1, grid % cells_n(4,c)-1,   &
        grid % cells_n(3,c)-1
    else
      print *, '# Unsupported cell type with ',  &
                  grid % cells_n_nodes(c), ' nodes.'
      print *, '# Exiting'
      stop 
    end if
  end do  
  write(9,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' offsets
  write(9,'(a,a)') IN_4, '<DataArray type="Int32" Name="offsets" format="ascii">'
  offset = 0
  do c = 1, grid % n_cells
    offset = offset + grid % cells_n_nodes(c)
    write(9,'(a,i9)') IN_5, offset
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'
 
  ! Now write all cells' types
  write(9,'(a,a)') IN_4, '<DataArray type="UInt8" Name="types" format="ascii">'
  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) == 8) then
      write(9,'(a,i9)') IN_5, VTK_HEXAHEDRON
    else if(grid % cells_n_nodes(c) == 6) then
      write(9,'(a,i9)') IN_5, VTK_WEDGE
    else if(grid % cells_n_nodes(c) == 4) then
      write(9,'(a,i9)') IN_5, VTK_TETRA
    else if(grid % cells_n_nodes(c) == 5) then
      write(9,'(a,i9)') IN_5, VTK_PYRAMID
    else
      print *, '# Unsupported cell type with ',  &
                  grid % cells_n_nodes(c), ' nodes.'
      print *, '# Exiting'
      stop 
    end if
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'
  write(9,'(a,a)') IN_3, '</Cells>'

  !---------------------------------!
  !                                 !
  !   Results and other cell data   !
  !                                 !
  !---------------------------------!
  write(9,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

  !---------------!
  !   Materials   !
  !---------------!
  write(9,'(a,a)') IN_4, '<DataArray type="UInt8" Name="materials"' //  &
                         ' format="ascii">'
  do c = 1, grid % n_cells
    write(9,'(a,i9)') IN_5, material(c)
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'

  !--------------!
  !   Velocity   !
  !--------------!
  write(9,'(a,a)') IN_4, '<DataArray type="Float32" Name="velocity" ' // &
                         'NumberOfComponents="3" format="ascii">'
  do c = 1, grid % n_cells
    write(9,'(a,1pe15.7,1pe15.7,1pe15.7)') IN_5, u % n(c), v % n(c), w % n(c)
  end do  
  write(9,'(a,a)') IN_4, '</DataArray>'


  !--------------!
  !   Pressure   !
  !--------------!
  write(9,'(a,a)') IN_4, '<DataArray type="Float32" Name="pressure"' //  &
                         ' format="ascii">'
  do c = 1, grid % n_cells
    write(9,'(a,1pe15.7)') IN_5, P % n(c)
  end do  
  write(9,'(a,a)') IN_4, '</DataArray>'

  !-----------------!
  !   Temperature   !
  !-----------------!
  if(HOT == YES) then
    write(9,'(a,a)') IN_4, '<DataArray type="Float32" Name="temperature"' //  &
                           ' format="ascii">'
    do c = 1, grid % n_cells
      write(9,'(a,1pe15.7)') IN_5, t % n(c)
    end do  
    write(9,'(a,a)') IN_4, '</DataArray>'
  end if 

  !----------------------------------------------------------!
  !   Turbulent quantities for single point closure models   !
  !----------------------------------------------------------!
  if(SIMULA == K_EPS.or.SIMULA==ZETA) then
    write(9,'(a,a)') '<DataArray type="Float32" Name="kinetic energy"' //  &
                   ' format="ascii">'
    do c = 1, grid % n_cells
      write(9,'(a,1pe15.7)') IN_5, kin % n(c)
    end do  
    write(9,'(a,a)') IN_4, '</DataArray>'

    write(9,'(a,a)') IN_5, '<DataArray type="Float32" Name="dissipation"' //  &
                   ' format="ascii">'
    do c = 1, grid % n_cells
      write(9,'(a,1pe15.7)') IN_5, eps % n(c)
    end do  
    write(9,'(a,a)') IN_4, '</DataArray>'
  end if

  write(9,'(a,a)') IN_3, '</CellData>'
  write(9,'(a,a)') IN_2, '</Piece>'
  write(9,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(9,'(a,a)') IN_0, '</VTKFile>'

  close(9)

  !-----------------------!
  !                       !
  !   Create .pvtu file   !
  !                       !
  !-----------------------!

  ! Create it only from processor 1, when running in parallel
  if(n_proc > 1 .and. this_proc == 1) then

    call Name_File(0, name_out, '.pvtu')
    print *, '# Creating the file: ', trim(name_out)
    open(9, file = name_out)
  
    ! Header
    write(9,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(9,'(a,a)') IN_0, '<VTKFile type="PUnstructuredGrid">'
    write(9,'(a,a)') IN_1, '<PUnstructuredGrid GhostLevel="0">'

    ! This section must be present
    write(9,'(a,a)') IN_2, '<PPoints>'
    write(9,'(a,a)') IN_3, '<PDataArray type="Float32" NumberOfComponents=' // &
                           '"3" format="ascii"/>'
    write(9,'(a,a)') IN_2, '</PPoints>'

    ! Data section is not mandatory, but very useful
    write(9,'(a,a)') IN_2, '<PCellData Scalars="scalars" vectors="velocity">'
    write(9,'(a,a)') IN_3, '<PDataArray type="UInt8" Name="materials"' //  & 
                           ' format="ascii"/>'
    write(9,'(a,a)') IN_3, '<PDataArray type="Float32" Name="velocity"' //  &
                           ' NumberOfComponents="3" format="ascii"/>'
    write(9,'(a,a)') IN_3, '<PDataArray type="Float32" Name="pressure"' //  &
                           ' format="ascii"/>'
    write(9,'(a,a)') IN_2, '</PCellData>'

    ! Write out the names of all the pieces
    do n = 1, n_proc
      call Name_File(n, name_out, '.vtu')
      write(9, '(a,a,a,a)') IN_2, '<Piece Source="', trim(name_out), '"/>'
    end do

    ! Footer
    write(9, '(a,a)'), IN_1, '</PUnstructuredGrid>'
    write(9, '(a,a)'), IN_0, '</VTKFile>'

    close(9)

  end if

  ! Restore the name
  problem_name = store_name

  end subroutine
