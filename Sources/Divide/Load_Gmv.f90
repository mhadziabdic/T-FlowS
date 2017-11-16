!==============================================================================!
  subroutine Load_Gmv(grid)
!------------------------------------------------------------------------------!
!   Reads: name.gmv                                                            !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use Tokenizer_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n
  character(len=80) :: dum_s, name_in
!==============================================================================!

  !---------------------------------!
  !     Read .gmv file with node    !
  !   coordinates and volume nodes  !
  !---------------------------------!
  name_in = name
  name_in(len_trim(name)+1:len_trim(name)+4) = '.gmv'
  write(*,*) '# Now reading the file:', name_in
  open(9, file=name_in)

  ! Read the number of nodes
  call Tokenizer_Mod_Read_Line(9)
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(2),*)   grid % n_nodes 

  !--------------------!
  !   Alocate memory   ! 
  !--------------------!
  write(*,'(A25)')       '# Allocating memory for: ' 
  write(*,'(A1,I8,A6)')  '#', grid % n_nodes,     ' nodes' 
  write(*,'(A1,I8,A6)')  '#', grid % n_cells,     ' cells' 
  write(*,'(A1,I8,A15)') '#', grid % n_bnd_cells, ' boundary cells'         
  write(*,'(A1,I8,A11)') '#', grid % n_faces,     ' cell faces' 

  ! Allocate nodes because you will soon be reading them
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes) 

  ! Variables defined in all_mod.h90:
  allocate (delta(-grid % n_bnd_cells:grid % n_cells));  delta=0.0
  allocate (WallDs(grid % n_faces));                     WallDs=0.0
  allocate (f(grid % n_faces));                          f=0.0

  ! Variables declared in gen_mod.h90:
  allocate (NewC(-grid % n_bnd_cells-1:grid % n_cells)); NewC = 0 
  allocate (NewS( grid % n_faces));                      NewS = 0
  allocate (NewN( grid % n_nodes));                      NewN = 0 

  ! Variables declared in div.h90:
  allocate (ix(-grid % n_bnd_cells:grid % n_cells));  ix=0
  allocate (iy(-grid % n_bnd_cells:grid % n_cells));  iy=0
  allocate (iz(-grid % n_bnd_cells:grid % n_cells));  iz=0
  allocate (iin(-grid % n_bnd_cells:grid % n_cells)); iin=0
  allocate (criter(grid % n_cells));   criter=0

  allocate (proces(grid % n_cells)); proces=0
  allocate (BuSeIn(grid % n_faces)); BuSeIn=0
  allocate (BuReIn(grid % n_faces)); BuReIn=0
  allocate (BufPos(grid % n_faces)); BufPos=0

  write(*,*) '# Allocation successfull !'

  ! Read node coordinates
  do n = 1, grid % n_nodes
    read(9,*) grid % xn(n)
  end do
  do n = 1, grid % n_nodes
    read(9,*) grid % yn(n)
  end do
  do n = 1, grid % n_nodes
    read(9,*) grid % zn(n)
  end do

  ! Read cell nodes 
  call Tokenizer_Mod_Read_Line(9)  ! cells, number of cells
  do c = 1, grid % n_cells  
    call Tokenizer_Mod_Read_Line(9)
    read(line % tokens(1), *) dum_s

    ! Hexahedral cells
    if(dum_s  ==  'hex') then
      grid % cells_n_nodes(c) = 8
      call Tokenizer_Mod_Read_Line(9)
      read(line % whole, *)                           &
           grid % cells_n(1,c), grid % cells_n(2,c),  &
           grid % cells_n(4,c), grid % cells_n(3,c),  &
           grid % cells_n(5,c), grid % cells_n(6,c),  &
           grid % cells_n(8,c), grid % cells_n(7,c) 

    ! Prismatic cells
    else if(dum_s  ==  'prism') then
      grid % cells_n_nodes(c) = 6
      call Tokenizer_Mod_Read_Line(9)
      read(line % whole, *)                         &
         grid % cells_n(1,c), grid % cells_n(2,c),  &
         grid % cells_n(3,c), grid % cells_n(4,c),  &
         grid % cells_n(5,c), grid % cells_n(6,c)

    ! Tetrahedral cells
    else if(dum_s  ==  'tet') then
      grid % cells_n_nodes(c) = 4
      call Tokenizer_Mod_Read_Line(9)
      read(line % whole, *)                         &
         grid % cells_n(1,c), grid % cells_n(2,c),  &
         grid % cells_n(3,c), grid % cells_n(4,c)

    ! Pyramid cells
    else if(dum_s  ==  'pyramid') then
      grid % cells_n_nodes(c) = 5
      call Tokenizer_Mod_Read_Line(9)
      read(line % whole, *)                         &
         grid % cells_n(5,c), grid % cells_n(1,c),  &
         grid % cells_n(2,c), grid % cells_n(4,c),  &
         grid % cells_n(3,c)

    ! Unsupported cells
    else
      write(*,*) '# Unsupported cell type: ', dum_s
      write(*,*) '# Exiting'
      stop
    end if
  end do

  ! Read cell materials
  call Tokenizer_Mod_Read_Line(9)  ! materials, number of materials 
  read(line % tokens(2),*) grid % n_materials
  do n = 1, grid % n_materials
    call Tokenizer_Mod_Read_Line(9)  
    grid % materials(n) % name = line % tokens(1)
  end do 
  do c = 1, grid % n_cells
    read(9,*) material(c)
  end do 

  close(9)

  end subroutine
