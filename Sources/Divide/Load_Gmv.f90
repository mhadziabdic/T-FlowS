!==============================================================================!
  subroutine Load_Gmv(grid)
!------------------------------------------------------------------------------!
!   Reads: name.gmv                                                            !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod 
  use div_mod
  use par_mod
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
  call ReadC(9,inp,tn,ts,te)
  call ReadC(9,inp,tn,ts,te)
  read(inp(ts(2):te(2)),*)   grid % n_nodes   ! number of nodes

  !--------------------!
  !   Alocate memory   ! 
  !--------------------!
  write(*,'(A25)')       '# Allocating memory for: ' 
  write(*,'(A1,I8,A6)')  '#', grid % n_nodes,     ' nodes' 
  write(*,'(A1,I8,A6)')  '#', grid % n_cells,     ' cells' 
  write(*,'(A1,I8,A15)') '#', grid % n_bnd_cells, ' boundary cells'         
  write(*,'(A1,I8,A11)') '#', grid % n_faces,     ' cell faces' 

  ! Variables defined in all_mod.h90:
  allocate (delta(-grid % n_bnd_cells:grid % n_cells));  delta=0.0
  allocate (WallDs(grid % n_faces));                     WallDs=0.0
  allocate (f(grid % n_faces));                          f=0.0

  ! Variables declared in gen_mod.h90:
  allocate (NewC(-grid % n_bnd_cells-1:grid % n_cells)); NewC = 0 
  allocate (NewS( grid % n_faces));                      NewS = 0
  allocate (NewN( grid % n_nodes));                      NewN = 0 

  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes) 
  call Grid_Mod_Allocate_Cells(grid, grid % n_bnd_cells, grid % n_cells) 
  call Grid_Mod_Allocate_Faces(grid, grid % n_faces) 


  ! Variables declared in div.h90:
  allocate (ix(-grid % n_bnd_cells:grid % n_cells));  ix=0
  allocate (iy(-grid % n_bnd_cells:grid % n_cells));  iy=0
  allocate (iz(-grid % n_bnd_cells:grid % n_cells));  iz=0
  allocate (iin(-grid % n_bnd_cells:grid % n_cells)); iin=0
  allocate (criter(grid % n_cells));   criter=0

  ! Variables declared in par_mod.h90:
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
  call ReadC(9,inp,tn,ts,te)  ! cells, number of cells
  do c = 1, grid % n_cells  !->>> ovo bi se moglo napravit neformatirano
    call ReadC(9,inp,tn,ts,te)
    read(inp(ts(1):te(1)),*) dum_s
    if(dum_s  ==  'hex') then
      grid % cells_n_nodes(c) = 8
      call ReadC(9,inp,tn,ts,te)
      read(inp,*)                                     &
           grid % cells_n(1,c), grid % cells_n(2,c),  &
           grid % cells_n(4,c), grid % cells_n(3,c),  &
           grid % cells_n(5,c), grid % cells_n(6,c),  &
           grid % cells_n(8,c), grid % cells_n(7,c) 
    else if(dum_s  ==  'prism') then
      grid % cells_n_nodes(c) = 6
      call ReadC(9,inp,tn,ts,te)
      read(inp,*)                                   &
         grid % cells_n(1,c), grid % cells_n(2,c),  &
         grid % cells_n(3,c), grid % cells_n(4,c),  &
         grid % cells_n(5,c), grid % cells_n(6,c)
    else if(dum_s  ==  'tet') then
      grid % cells_n_nodes(c) = 4
      call ReadC(9,inp,tn,ts,te)
      read(inp,*)                                   &
         grid % cells_n(1,c), grid % cells_n(2,c),  &
         grid % cells_n(3,c), grid % cells_n(4,c)
    else if(dum_s  ==  'pyramid') then
      grid % cells_n_nodes(c) = 5
      call ReadC(9,inp,tn,ts,te)
      read(inp,*)                                   &
         grid % cells_n(5,c), grid % cells_n(1,c),  &
         grid % cells_n(2,c), grid % cells_n(4,c),  &
         grid % cells_n(3,c)
    else
      write(*,*) '# Unsupported cell type: ', dum_s
      write(*,*) '# Exiting'
      stop
    end if
  end do

  ! Read cell materials
  call ReadC(9,inp,tn,ts,te)  ! materials, number of materials 
  read(inp(ts(2):te(2)),*) grid % n_materials
  do n = 1, grid % n_materials
    call ReadC(9,inp,tn,ts,te)  
    grid % materials(n) % name = inp(ts(1):te(1))
  end do 
  do c = 1, grid % n_cells
    read(9,*) material(c)
  end do 

  close(9)

  end subroutine
