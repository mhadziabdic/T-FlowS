!==============================================================================!
  subroutine Load_Neu() 
!------------------------------------------------------------------------------!
!   Reads the Fluents (Gambits) neutral file format.                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use gen_mod 
  use neu_mod 
  use par_mod 
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=130)  :: name_in
  integer             :: i, j, n_blocks, n_bnd_sect, dum1, dum2
  integer,allocatable :: temp(:)
  integer             :: c, dir
!==============================================================================!

  name_in = name
  name_in(len_trim(name)+1:len_trim(name)+4) = '.neu'

  open(9,file=name_in)
  write(*,*) '# Now reading the file: ', name_in

  ! Skip first 6 lines
  do i=1,6
    call Read_Line(9,inp,tn,ts,te)
  end do 

  ! Read the line which contains usefull information  
  call Read_Line(9,inp,tn,ts,te)
  read(inp(ts(1):te(1)),*) NN  
  read(inp(ts(2):te(2)),*) NC
  read(inp(ts(3):te(3)),*) n_blocks
  read(inp(ts(4):te(4)),*) n_bnd_sect

  write(*,*) '# Total number of nodes:  ',            NN
  write(*,*) '# Total number of cells:  ',            NC
  write(*,*) '# Total number of blocks: ',            n_blocks
  write(*,*) '# Total number of boundary sections: ', n_bnd_sect

  ! Count the boundary cells
  NbC = 0
  do 
    call Read_Line(9,inp,tn,ts,te)
    if( inp(ts(1):te(1)) == 'BOUNDARY' ) then
      do j=1,n_bnd_sect
        if(j>1) call Read_Line(9,inp,tn,ts,te) ! BOUNDARY CONDITIONS
        call Read_Line(9,inp,tn,ts,te)
        read(inp(ts(3):te(3)),*) dum1  
        NBc = NBc + dum1 
        do i=1,dum1
          read(9,*) c, dum2, dir
        end do
        call Read_Line(9,inp,tn,ts,te)         ! ENDOFSECTION
      end do
      write(*,*) '# Total number of boundary cells: ', NbC
      go to 1
    end if
  end do 

1 rewind(9)

  ! Skip first 7 lines
  do i=1,7
    call Read_Line(9,inp,tn,ts,te)
  end do 

  ! Allocate memory =--> carefull, there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, NN) 
  call Grid_Mod_Allocate_Cells(grid, NbC, NC) 
  call Grid_Mod_Allocate_Faces(grid, NC*5) 

  allocate(material(-NbC:NC));     material=0 
  allocate(BCtype(NC,6));          BCtype=0
  allocate(BCmark(-NbC-1:-1));     BCmark=0
! WARNING:  allocate(CellN(-NbC-1:NC,-1:8)); CellN=0
  allocate(SideC(0:2,NC*5));       SideC=0    
  n_copy = 1000000 
  allocate(CopyC(-NbC:-1));  CopyC=0
  allocate(CopyS(2,n_copy));  CopyS=0

  allocate(NewN(NN));        NewN=0  
  allocate(NewC(-NbC-1:NC)); NewC=0  
  allocate(NewS(NC*5));      NewS=0  

  allocate(proces(0:NC)); proces=0

  allocate(temp(NC)); temp=0

  ! Skip one line 
  call Read_Line(9,inp,tn,ts,te)

  !--------------------------------!
  !   Read the nodal coordinates   !
  !--------------------------------!
  call Read_Line(9,inp,tn,ts,te)          ! NODAL COORDINATES
  do i=1,NN
    call Read_Line(9,inp,tn,ts,te)
    read(inp(ts(2):te(2)),*) grid % xn(i) 
    read(inp(ts(3):te(3)),*) grid % yn(i)
    read(inp(ts(4):te(4)),*) grid % zn(i)
  end do
  call Read_Line(9,inp,tn,ts,te)          ! ENDOFSECTION

  !-----------------------------!
  !   Read nodes of each cell   !
  !-----------------------------!
  call Read_Line(9,inp,tn,ts,te)          ! ELEMENTS/CELLS
  do i=1,NC
    read(9,'(I8,1X,I2,1X,I2,1X,7I8:/(15X,7I8:))') dum1, dum2, &
           grid % cells_n_nodes(i),                           &
          (grid % cells_n(j,i), j = 1,grid % cells_n_nodes(i))
  end do
  call Read_Line(9,inp,tn,ts,te)          ! ENDOFSECTION

  !---------------------------------!
  !   Read block (material?) data   !
  !---------------------------------!
  allocate(grid % materials(n_blocks))
  Nmat = 1
  grid % materials(Nmat) % name = "FLUID"

  do j=1,n_blocks
    call Read_Line(9,inp,tn,ts,te)        ! ELEMENT GROUP
    call Read_Line(9,inp,tn,ts,te)
    read(inp(ts(4):te(4)),'(I10)') dum1  
    call Read_Line(9,inp,tn,ts,te)        ! block*
    call Read_Line(9,inp,tn,ts,te)        ! 0
    read(9,'(10I8)') (temp(i), i=1,dum1)
    do i=1,dum1
      material(temp(i)) = j
    end do
    call Read_Line(9,inp,tn,ts,te)        ! ENDOFSECTION
  end do

  !-------------------------!
  !   Boundary conditions   !
  !-------------------------!
  Nbnd = n_bnd_sect
  allocate(grid % boundary_conditions(n_bnd_sect))

  do j=1,n_bnd_sect
    call Read_Line(9,inp,tn,ts,te)        ! BOUNDARY CONDITIONS
    call Read_Line(9,inp,tn,ts,te)
    call To_Upper_Case(  inp(ts(1):te(1))  )
    write(*,*)  inp(ts(1):te(1)), j 
    grid % boundary_conditions(j) % name = inp(ts(1):te(1))
    read(inp(ts(3):te(3)),'(I8)') dum1  
    do i=1,dum1
      read(9,*) c, dum2, dir
      BCtype(c,dir) = j 
    end do
    call Read_Line(9,inp,tn,ts,te)        ! ENDOFSECTION
  end do

  write(*,*) '#==================================================='
  write(*,*) '# Found the following boundary conditions:'
  write(*,*) '#---------------------------------------------------'
  do j=1,n_bnd_sect
    write(*,*) '# ', grid % boundary_conditions(j) % name
  end do
  write(*,*) '#---------------------------------------------------'

  close(9)

  end subroutine Load_Neu
