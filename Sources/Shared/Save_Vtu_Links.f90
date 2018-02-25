!==============================================================================!
  subroutine Save_Vtu_Links(grid, sub, NNsub, NCsub, NSsub, NBCsub, NBFsub) 
!------------------------------------------------------------------------------!
!   Creates the file "name.ln.vtu" to check the cell connections.              !
!                                                                              !
!   Links between the computational cells have been introduced as aditional    !
!   cells of general type. Cell centers are introduced as aditional nodes.     !
!   Material of these links is different than from the cells, so that they     !
!   can be visualised  more easily in GMV.                                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use div_mod,  only: BuSeIn, BuReIn
  use gen_mod,  only: new_n, new_c, new_f
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, NNsub, NCsub, NSsub, NBCsub, NBFsub
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n, c, c1, c2, s, offset
  integer           :: nf_sub_non_per, nf_sub_per
  character(len=80) :: name_out
!------------------------------[Local parameters]------------------------------!
  integer,           parameter :: VTK_LINE       =  3  ! cells in VTK format
  integer,           parameter :: VTK_TETRA      = 10
  integer,           parameter :: VTK_HEXAHEDRON = 12  
  integer,           parameter :: VTK_WEDGE      = 13
  integer,           parameter :: VTK_PYRAMID    = 14
  character(len= 0), parameter :: IN_0 = ''            ! indentation levels 
  character(len= 2), parameter :: IN_1 = '  '
  character(len= 4), parameter :: IN_2 = '    '
  character(len= 6), parameter :: IN_3 = '      '
  character(len= 8), parameter :: IN_4 = '        '
  character(len=10), parameter :: IN_5 = '          '
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .gmv file   !
  !                      !
  !----------------------!
  name_out = problem_name         

  call Name_File(sub, name_out, '.links.vtu')
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
                   IN_2, '<Piece NumberOfPoints="',               &
                          NNsub + NCsub + NBCsub + NBFsub,        &
                         '" NumberOfCells ="',                    &
                          NCsub + NSsub + NBFsub, '">'  ! + NBFsub
  !-----------!
  !   Nodes   !
  !-----------!
  write(9,'(a,a)') IN_3, '<Points>'
  write(9,'(a,a)') IN_4, '<DataArray type="Float32" NumberOfComponents' //  &
                         '="3" format="ascii">'
  do n = 1, grid % n_nodes
    if(new_n(n) /= 0) write(9, '(a,1pe15.7,1pe15.7,1pe15.7)')                &
                                IN_5, grid % xn(n), grid % yn(n), grid % zn(n)
  end do
  do c = 1, grid % n_cells
    if(new_c(c) /= 0) write(9, '(a,1pe15.7,1pe15.7,1pe15.7)')                &
                                IN_5, grid % xc(c), grid % yc(c), grid % zc(c)
  end do
  do c = -1,-grid % n_bnd_cells,-1
    if(new_c(c) /= 0) write(9, '(a,1pe15.7,1pe15.7,1pe15.7)')                &
                                IN_5, grid % xc(c), grid % yc(c), grid % zc(c)
  end do
  do c = 1,NBFsub
    write(9, '(a,1pe15.7,1pe15.7,1pe15.7)') IN_5,                           &
               grid % xc(BuReIn(c)), grid % yc(BuReIn(c)), grid % zc(BuReIn(c))
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
    if(new_c(c) /= 0) then

      ! Hexahedral
      if(grid % cells_n_nodes(c) == 8) then
        write(9,'(a,8i9)')                                           &
          IN_5,                                                      &
          new_n(grid % cells_n(1,c))-1, new_n(grid % cells_n(2,c))-1,  &
          new_n(grid % cells_n(4,c))-1, new_n(grid % cells_n(3,c))-1,  &
          new_n(grid % cells_n(5,c))-1, new_n(grid % cells_n(6,c))-1,  &
          new_n(grid % cells_n(8,c))-1, new_n(grid % cells_n(7,c))-1

      ! Wedge       
      else if(grid % cells_n_nodes(c) == 6) then
        write(9,'(a,6i9)')                                           &
          IN_5,                                                      &
          new_n(grid % cells_n(1,c))-1, new_n(grid % cells_n(2,c))-1,  &
          new_n(grid % cells_n(3,c))-1, new_n(grid % cells_n(4,c))-1,  &
          new_n(grid % cells_n(5,c))-1, new_n(grid % cells_n(6,c))-1

      ! Tetrahedra  
      else if(grid % cells_n_nodes(c) == 4) then
        write(9,'(a,4i9)')                                           &
          IN_5,                                                      &
          new_n(grid % cells_n(1,c))-1, new_n(grid % cells_n(2,c))-1,  &
          new_n(grid % cells_n(3,c))-1, new_n(grid % cells_n(4,c))-1

      ! Pyramid     
      else if(grid % cells_n_nodes(c) == 5) then
        write(9,'(a,5i9)')                                           &
          IN_5,                                                      &
          new_n(grid % cells_n(1,c))-1, new_n(grid % cells_n(2,c))-1,  &
          new_n(grid % cells_n(4,c))-1, new_n(grid % cells_n(3,c))-1,  &
          new_n(grid % cells_n(5,c))-1
      else
        print *, '# Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), ' nodes.'
        print *, '# Exiting'
        stop 
      end if
    end if
  end do  

  ! Physical links; non-periodic
  nf_sub_non_per = 0 
  do s=1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if( (new_f(s) > 0) .and. (new_f(s) <= NSsub) ) then

      if( (grid % sx(s) * (grid % xc(c2)-grid % xc(c1) ) +  &
           grid % sy(s) * (grid % yc(c2)-grid % yc(c1) ) +  &
           grid % sz(s) * (grid % zc(c2)-grid % zc(c1) ))  > 0.0 ) then 

        nf_sub_non_per = nf_sub_non_per + 1

        c1 = new_c(grid % faces_c(1,s))
        c2 = new_c(grid % faces_c(2,s))
        if( c2  > 0 ) then
          write(9,'(a,2i9)') IN_5, NNsub+c1-1, NNsub+c2-1
        else
          write(9,'(a,2i9)') IN_5, NNsub+c1-1, NNsub+NCsub-c2-1
        end if
      end if

    end if
  end do  

  ! Physical links; periodic
  nf_sub_per    = 0 
  do s=1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if( (new_f(s) > 0) .and. (new_f(s) <= NSsub) ) then

      if( (grid % sx(s) * (grid % xc(c2)-grid % xc(c1) ) +  &
           grid % sy(s) * (grid % yc(c2)-grid % yc(c1) ) +  &
           grid % sz(s) * (grid % zc(c2)-grid % zc(c1) ))  < 0.0 ) then 

        nf_sub_per = nf_sub_per + 1

        c1 = new_c(grid % faces_c(1,s))
        c2 = new_c(grid % faces_c(2,s))
        if( c2  > 0 ) then
          write(9,'(a,2i9)') IN_5, NNsub+c1-1, NNsub+c2-1                                
        else
          write(9,'(a,2i9)') IN_5, NNsub+c1-1, NNsub+NCsub-c2-1
        end if
      end if

    end if
  end do  

  ! Interprocessor links
  do c = 1, NBFsub
    c1 = BuSeIn(c) 
    write(9,'(a,2i9)') IN_5, NNsub+c1-1, NNsub+NCsub+NBCsub+c-1
  end do  

  write(9,'(a,a)') IN_4, '</DataArray>'

  print *, '# Non-periodic links    :', nf_sub_non_per
  print *, '# Periodic links        :', nf_sub_per
  print *, '# Inter-processor links :', NBFsub

  ! Now write all cells' offsets
  write(9,'(a,a)') IN_4, '<DataArray type="Int32" Name="offsets" format="ascii">'
  offset = 0
  do c = 1, grid % n_cells
    if(new_c(c) /= 0) then
      offset = offset + grid % cells_n_nodes(c)
      write(9,'(a,i9)') IN_5, offset
    end if
  end do
  do c = 1, nf_sub_non_per
    offset = offset + 2
    write(9,'(a,i9)') IN_5, offset
  end do
  do c = 1, nf_sub_per
    offset = offset + 2
    write(9,'(a,i9)') IN_5, offset
  end do
  do c = 1, NBFsub
    offset = offset + 2
    write(9,'(a,i9)') IN_5, offset
  end do

  write(9,'(a,a)') IN_4, '</DataArray>'
 
  ! Now write all cells' types
  write(9,'(a,a)') IN_4, '<DataArray type="UInt8" Name="types" format="ascii">'
  do c = 1, grid % n_cells
    if(new_c(c) /= 0) then
      if(grid % cells_n_nodes(c) == 4) write(9,'(a,i9)') IN_5, VTK_TETRA
      if(grid % cells_n_nodes(c) == 8) write(9,'(a,i9)') IN_5, VTK_HEXAHEDRON
      if(grid % cells_n_nodes(c) == 6) write(9,'(a,i9)') IN_5, VTK_WEDGE
      if(grid % cells_n_nodes(c) == 5) write(9,'(a,i9)') IN_5, VTK_PYRAMID
    end if
  end do
  do c = 1, nf_sub_non_per
    write(9,'(a,i9)') IN_5, VTK_LINE
  end do
  do c = 1, nf_sub_per
    write(9,'(a,i9)') IN_5, VTK_LINE
  end do
  do c = 1, NBFsub
    write(9,'(a,i9)') IN_5, VTK_LINE
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'
  write(9,'(a,a)') IN_3, '</Cells>'
 
  !----------------!
  !   Link types   !
  !----------------!
  write(9,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'
  write(9,'(a,a)') IN_4, '<DataArray type="UInt8" Name="link type" format="ascii">'
  do c = 1, grid % n_cells
    if(new_c(c) /= 0) then
      write(9,'(a,i9)') IN_5, 0
    end if
  end do
  do c = 1, nf_sub_non_per
    write(9,'(a,i9)') IN_5, 1
  end do
  do c = 1, nf_sub_per
    write(9,'(a,i9)') IN_5, 2
  end do
  do c = 1, NBFsub
    write(9,'(a,i9)') IN_5, 3
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'
  write(9,'(a,a)') IN_3, '</CellData>'

  !------------!
  !   Footer   !
  !------------!
  write(9,'(a,a)') IN_2, '</Piece>'
  write(9,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(9,'(a,a)') IN_0, '</VTKFile>'

  close(9)

  end subroutine
