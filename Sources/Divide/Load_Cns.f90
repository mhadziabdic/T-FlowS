!==============================================================================!
  subroutine Load_Cns(grid, this_proc)
!------------------------------------------------------------------------------!
!   Reads: .cns file.                                                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: this_proc  ! needed if called from Processor
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n, s
  character(len=80) :: name_in
!==============================================================================!

  !-------------------------------!
  !     Read the file with the    !
  !   connections between cells   !
  !-------------------------------!
  call Name_File(this_proc, name_in, '.cns', len_trim('.cns'))  

  open(9, file=name_in,form='unformatted')
  if(this_proc < 2) write(*,*) '# Now reading the file:', name_in

  ! Number of cells, boundary cells and sides
  read(9) grid % n_nodes
  read(9) grid % n_cells
  read(9) grid % n_bnd_cells
  read(9) grid % n_faces
  read(9) grid % n_sh

  ! Allocate memory =--> carefull, there is no checking!
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells, grid % n_bnd_cells) 
  call Grid_Mod_Allocate_Faces(grid, grid % n_faces, grid % n_sh) 

  ! Number of materials and boundary conditions
  read(9) grid % n_materials
  read(9) grid % n_boundary_conditions

  allocate(grid % materials          (grid % n_materials))
  allocate(grid % boundary_conditions(grid % n_boundary_conditions))

  ! Materials' and boundary conditions' keys
  do n = 1, grid % n_materials
    read(9) grid % materials(n) % name
  end do
  do n = 1, grid % n_boundary_conditions
    read(9) grid % boundary_conditions(n) % name
  end do

  ! Cell materials
  allocate (material(-grid % n_bnd_cells:grid % n_cells))
  read(9) (material(c), c =  1, grid % n_cells)
  read(9) (material(c), c = -1,-grid % n_bnd_cells,-1)

  ! Faces
  read(9) (grid % faces_c(1,s), s = 1, grid % n_faces)
  read(9) (grid % faces_c(2,s), s = 1, grid % n_faces)

  ! Boundary cells
  allocate (bcmark(-grid % n_bnd_cells-1:-1))
  read(9) (bcmark(c), c = -1,-grid % n_bnd_cells, -1)

  ! Boundary copy cells
  allocate (CopyC(-grid % n_bnd_cells:-1))
  read(9) (CopyC(c), c = -1,-grid % n_bnd_cells, -1)
 
  read(9) grid % n_copy
  write(*,*) '# Number of copy cells/faces: ', grid % n_copy
  allocate (CopyS(2,grid % n_copy))
  read(9) (CopyS(1,s), s = 1,grid % n_copy)
  read(9) (CopyS(2,s), s = 1,grid % n_copy)

  close(9)

  end subroutine
