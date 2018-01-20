!==============================================================================!
  subroutine Save_Vtu_Results(grid, name_save)
!------------------------------------------------------------------------------!
!   Writes results in VTU file format (for VisIt and Paraview)                 !
!----------------------------------[Modules]-----------------------------------!
  use allp_mod
  use all_mod
  use pro_mod
  use rans_mod
  use par_mod, only: this_proc
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
  call Name_File(this_proc, name_out, '.vtu', len_trim('.vtu'))
  open(9, file=name_out)
  write(6, *) '# Now creating the file:', trim(name_out)

  !-----------!
  !   Start   !
  !-----------!
  write(9,'(a)') '<?xml version="1.0"?>'
  write(9,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                 'byte_order="LittleEndian">'
  write(9,'(a)') '<UnstructuredGrid>'
  write(9,'(a,i9)',advance='no') '<Piece NumberOfPoints="', grid % n_nodes
  write(9,'(a,i9,a)')            '" NumberOfCells ="', grid % n_cells, '">'

  !----------!
  !          !
  !   Grid   !
  !          !
  !----------!

  !-----------!
  !   Nodes   !
  !-----------!
  write(9,'(a)') '<Points>'
  write(9,'(a)') '<DataArray type="Float32" NumberOfComponents="3" ' //  &
                 'format="ascii">'
  do n = 1, grid % n_nodes
    write(9, '(1PE15.7, 1PE15.7, 1PE15.7)')  &
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
    if(grid % cells_n_nodes(c) == 8) then
      write(9,'(8I9)')                                  &
        grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
        grid % cells_n(4,c)-1, grid % cells_n(3,c)-1,   &
        grid % cells_n(5,c)-1, grid % cells_n(6,c)-1,   &
        grid % cells_n(8,c)-1, grid % cells_n(7,c)-1
    else if(grid % cells_n_nodes(c) == 6) then
      write(9,'(6I9)')                                  &
        grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
        grid % cells_n(3,c)-1, grid % cells_n(4,c)-1,   &
        grid % cells_n(5,c)-1, grid % cells_n(6,c)-1
    else if(grid % cells_n_nodes(c) == 4) then
      write(9,'(4I9)')                                  &
        grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
        grid % cells_n(3,c)-1, grid % cells_n(4,c)-1
    else if(grid % cells_n_nodes(c) == 5) then
      write(9,'(5I9)')                                  &
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
  write(9,'(a)') '</DataArray>'

  ! Now write all cells' offsets
  write(9,'(a)') '<DataArray type="Int32" Name="offsets" format="ascii">'
  offset = 0
  do c = 1, grid % n_cells
    offset = offset + grid % cells_n_nodes(c)
    write(9,'(i9)') offset
  end do
  write(9,'(a)') '</DataArray>'
 
  ! Now write all cells' types
  write(9,'(a)') '<DataArray type="UInt8" Name="types" format="ascii">'
  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) == 8) then  ! hexahedral cells
      write(9,'(i9)') 12
    else if(grid % cells_n_nodes(c) == 6) then  ! prismatic cells
      write(9,'(i9)') 13
    else if(grid % cells_n_nodes(c) == 4) then  ! tetrahedral cells
      write(9,'(i9)') 10
    else if(grid % cells_n_nodes(c) == 5) then  ! pyramid cells
      write(9,'(i9)') 14
    else
      print *, '# Unsupported cell type with ',  &
                  grid % cells_n_nodes(c), ' nodes.'
      print *, '# Exiting'
      stop 
    end if
  end do
  write(9,'(a)') '</DataArray>'
 
!  !---------------!
!  !   Materials   !
!  !---------------!
!  write(9,'(A10,2I5)') 'materials', grid % n_materials, 0
!  do n = 1, grid % n_materials
!    write(9,'(a)') grid % materials(n) % name
!  end do
!  do c = 1, grid % n_cells
!    if(NewC(c) /= 0) then
!      write(9,'(a)') material(c)
!    end if
!  end do        

  !-------------!
  !             !
  !   Results   !
  !             !
  !-------------!

  !--------------!
  !   Velocity   !
  !--------------!
  write(9,'(a)') '<DataArray type="Float32" Name="velocity" ' // &
                 'NumberOfComponents="3" format="ascii">'
  do c = 1, grid % n_cells
    write(9, '(1PE15.7, 1PE15.7, 1PE15.7)') u % n(c), v % n(c), w % n(c)
  end do  
  write(9,'(a)') '</DataArray>'

  write(9,'(a)') '</Cells>'
  write(9,'(a)') '</Piece>'
  write(9,'(a)') '</UnstructuredGrid>'
  write(9,'(a)') '</VTKFile>'

  !--------------!
  !   Pressure   !
  !--------------!
  write(9,'(a)') '<DataArray type="Float32" Name="pressure" format="ascii">'
  do c = 1, grid % n_cells
    write(9,*) P % n(c)
  end do  
  write(9,'(a)') '</DataArray>'

  !-----------------!
  !   Temperature   !
  !-----------------!
  if(HOT == YES) then
    write(9,'(a)') '<DataArray type="Float32" Name="temperature" ' //  &
                   'format="ascii">'
    do c = 1, grid % n_cells
      write(9,*) t % n(c)
    end do  
    write(9,'(a)') '</DataArray>'
  end if 

  !----------------------------------------------------------!
  !   Turbulent quantities for single point closure models   !
  !----------------------------------------------------------!
  if(SIMULA == K_EPS.or.SIMULA==ZETA) then
    write(9,'(a)') '<DataArray type="Float32" Name="kinetic energy" ' //  &
                   'format="ascii">'
    do c = 1, grid % n_cells
      write(9,*) kin % n(c)
    end do  
    write(9,'(a)') '</DataArray>'

    write(9,'(a)') '<DataArray type="Float32" Name="dissipation" ' //  &
                   'format="ascii">'
    do c = 1, grid % n_cells
      write(9,*) eps % n(c)
    end do  
    write(9,'(a)') '</DataArray>'
  end if

  close(9)

  ! Restore the name
  problem_name = store_name

  end subroutine
