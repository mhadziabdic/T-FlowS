!==============================================================================!
  subroutine Load_Cgns(grid) 
!------------------------------------------------------------------------------!
!   Reads the Fluents (Gambits) neutral file format.                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use gen_mod 
  use neu_mod 
  use Grid_Mod
  use Tokenizer_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include "cgns_io_f.h"
  include "cgnslib_f.h"
!-----------------------------------[Locals]-----------------------------------!
  character(len=80)       :: name_in
  integer                 :: i, j, n_blocks, n_bnd_sect, dum1, dum2
  integer,allocatable     :: temp(:)
  integer                 :: c, n, dir
  integer                 :: index_file, index_base, index_zone, index_sect, ier
  integer                 :: n_sect, n_bnd
  integer, dimension(3,3) :: isize
  integer, dimension(3)   :: irmin, irmax
  integer                 :: istart, iend, itype, iparent_data, iparent_flag
  character(len=80)       :: name_zone, name_sect
!==============================================================================!

  name_in = problem_name
  name_in(len_trim(problem_name)+1:len_trim(problem_name)+5) = '.cgns'

  !--------------------!
  !   Open CGNS file   !
  !--------------------!
  print *, '# Reading the file: ', trim(name_in)
  call Cg_Open_F(name_in,         &  ! file name
                 CGIO_MODE_READ,  &  ! mode
                 index_file,      &  ! file handle
                 ier)

  ! Assume there is only one zone and one base, 
  ! but in the future it should be checked
  index_base = 1
  index_zone = 1

  ! Read size of the zone (number of nodes, cells, boundary cells???)
  call Cg_Zone_Read_F(index_file,      &  ! file handle
                      index_base,      &  ! assumed to be 1
                      index_zone,      &  ! assumed to be 1
                      name_zone,       &  ! returns name of the zone
                      isize,           &  ! number of nodes I hope
                      ier)
  print *, '# Zone name: ', name_zone

  grid % n_nodes     = isize(1,1)
  grid % n_cells     = isize(2,1)
  grid % n_bnd_cells = isize(3,1) ! ???

  ! Allocate memory =--> carefull, there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes) 
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells, grid % n_bnd_cells) 

  !---------------------------------!
  !   Read x, y and z coordinates   !
  !---------------------------------!

  ! Variables irmin and irmax have to be set to proper 
  ! ranges to fetch information about nodes and cells
  irmin(1) = 1
  irmin(2) = 1
  irmin(3) = 1
  irmax(1) = isize(1,1)
  irmax(2) = isize(2,1)
  irmax(3) = isize(3,1)

  ! Do actual reading
  call Cg_Coord_Read_F(index_file,     &  ! file handle
                       index_base,     &  ! assumed to be 1
                       index_zone,     &  ! assumed to be 1
                       'CoordinateX',  &  ! take "x" coordinate
                       RealDouble,     &  ! data type
                       irmin, irmax,   &  ! range?
                       grid % xn,      &  ! where to store 
                       ier)
  call Cg_Coord_Read_F(index_file,     &  ! file handle
                       index_base,     &  ! assumed to be 1
                       index_zone,     &  ! assumed to be 1
                       'CoordinateY',  &  ! take "x" coordinate
                       RealDouble,     &  ! data type
                       irmin, irmax,   &  ! range?
                       grid % yn,      &  ! where to store 
                       ier)
  call Cg_Coord_Read_F(index_file,     &  ! file handle
                       index_base,     &  ! assumed to be 1
                       index_zone,     &  ! assumed to be 1
                       'CoordinateZ',  &  ! take "x" coordinate
                       RealDouble,     &  ! data type
                       irmin, irmax,   &  ! range?
                       grid % zn,      &  ! where to store 
                       ier)

! do n = 1, grid % n_nodes
!   print *, grid % xn(n), grid % yn(n), grid % zn(n)
! end do

  !----------------!
  !   Read cells   !
  !----------------!

  ! Read number of sections
  call Cg_Nsections_F(index_file,  &
                      index_base,  &
                      index_zone,  &
                      n_sect,  &
                      ier)
  print *, '# Number of sections: ', n_sect

  ! Browse through all sections to read elements
  do index_sect = 1, n_sect

    call Cg_Section_Read_F(index_file,    &
                           index_base,    &
                           index_zone,    &
                           index_sect,    &  ! loop counter
                           name_sect,     &  ! name of the section
                           itype,         &  ! element type index
                           istart, iend,  &  ! section range
                           n_bnd,         &  
                           iparent_flag,  &
                           ier)
    print *, '# Reading section data...'
    print *, '#    section name: ', name_sect
    print *, '#    section type: ', ElementTypeName(itype)
    print *, '#    range:        ', istart, iend             

    ! Read cells' nodes for the whole section
    call Cg_Elements_Read_F(index_file,      &
                            index_base,      &
                            index_zone,      &
                            index_sect,      &
                            grid % cells_n,  &
                            iparent_data,    &
                            ier)
  end do 

  do c = 1, grid % n_cells
    print *, '# Cell: ', c
    print *, '# Nodes: ', (grid % cells_n(n,c), n=1,8)
  end do

stop

  end subroutine
