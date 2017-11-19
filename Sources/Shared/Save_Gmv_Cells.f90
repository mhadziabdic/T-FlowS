!==============================================================================!
  subroutine Save_Gmv_Cells(grid, sub, NNsub, NCsub)
!------------------------------------------------------------------------------!
! Writes: NAME.gmv, NAME.faces.gmv, NAME.shadow.gmv                            !
!----------------------------------[Modules]-----------------------------------!
  use all_mod, only: material
  use gen_mod, only: NewN, NewC
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, NNsub, NCsub
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n
  character(len=80) :: name_out
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .gmv file   !
  !                      !
  !----------------------!
  call Name_File(sub, name_out, '.gmv', len_trim('.gmv'))
  open(9, file=name_out)
  write(6, *) '# Now creating the file:', trim(name_out)

  !-----------!
  !   Start   !
  !-----------!
  write(9,'(A14)') 'gmvinput ascii'

  !-----------!
  !   Nodes   !
  !-----------!
  write(9,*) 'nodes', NNsub

  do n = 1, grid % n_nodes
    if(NewN(n) /= 0) write(9, '(1PE14.7)') grid % xn(n)
  end do
  do n = 1, grid % n_nodes
    if(NewN(n) /= 0) write(9, '(1PE14.7)') grid % yn(n)
  end do
  do n = 1, grid % n_nodes
    if(NewN(n) /= 0) write(9, '(1PE14.7)') grid % zn(n)
  end do

  !-----------!
  !   Cells   !
  !-----------!
  write(9,*) 'cells', NCsub
  do c = 1, grid % n_cells
    if(NewC(c) /= 0) then
      if(grid % cells_n_nodes(c) == 8) then
        write(9,*) 'hex 8'
        write(9,'(8I9)')                                          &
          NewN(grid % cells_n(1,c)), NewN(grid % cells_n(2,c)),   &
          NewN(grid % cells_n(4,c)), NewN(grid % cells_n(3,c)),   &
          NewN(grid % cells_n(5,c)), NewN(grid % cells_n(6,c)),   &
          NewN(grid % cells_n(8,c)), NewN(grid % cells_n(7,c))
      else if(grid % cells_n_nodes(c) == 6) then
        write(9,*) 'prism 6'
        write(9,'(6I9)')                                          &
          NewN(grid % cells_n(1,c)), NewN(grid % cells_n(2,c)),   &
          NewN(grid % cells_n(3,c)), NewN(grid % cells_n(4,c)),   &
          NewN(grid % cells_n(5,c)), NewN(grid % cells_n(6,c))
      else if(grid % cells_n_nodes(c) == 4) then
        write(9,*) 'tet 4'
        write(9,'(4I9)')                                          &
          NewN(grid % cells_n(1,c)), NewN(grid % cells_n(2,c)),   &
          NewN(grid % cells_n(3,c)), NewN(grid % cells_n(4,c))
      else if(grid % cells_n_nodes(c) == 5) then
        write(9,*) 'pyramid 5'
        write(9,'(5I9)')                                          &
          NewN(grid % cells_n(5,c)), NewN(grid % cells_n(1,c)),   &
          NewN(grid % cells_n(2,c)), NewN(grid % cells_n(4,c)),   &
          NewN(grid % cells_n(3,c))
      else
        write(*,*) '# Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), ' nodes.'
        write(*,*) '# Exiting'
        stop 
      end if 
    end if
  end do  

  !---------------!
  !   Materials   !
  !---------------!
  write(9,'(A10,2I5)') 'materials', grid % n_materials, 0
  do n = 1, grid % n_materials
    write(9,*) grid % materials(n) % name
  end do
  do c = 1, grid % n_cells
    if(NewC(c) /= 0) then
      write(9,*) material(c)
    end if
  end do        

  write(9,'(A6)') 'endgmv'            !  end the GMV file
  close(9)

  end subroutine
