!==============================================================================!
  subroutine Load_Cgns(grid)
!------------------------------------------------------------------------------!
!   Reads the Fluents (Gambits) neutral file format.                           !
!------------------------------------------------------------------------------!
!   https://cgns.github.io/CGNS_docs_current/midlevel/structural.html
!   |-> mesh_info          !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use neu_mod
  use Grid_Mod
  use Tokenizer_Mod
  use cgns_mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  character(len=80)       :: name_in

!==============================================================================!

  name_in = problem_name

  name_in(len_trim(problem_name)+1:len_trim(problem_name)+5) = ".cgns"

  file_name = name_in

  !--------------------------------------!
  !   First run: just read info from DB  !
  !--------------------------------------!

  ! Open a CGNS file (->file_id)
  call Cgns_Mod_Open_File

  ! Read number of CGNS bases in file_id (->n_bases)
  call Cgns_Mod_Read_Number_Of_Bases_In_File

  !------------------------------!
  !   Browse through all bases   !
  !------------------------------!
  do base_id = 1, n_bases

    ! Print CGNS base information in base_id
    call Cgns_Mod_Print_Base_Info

    ! Read number of zones in base_id(->n_zones)
    call Cgns_Mod_Read_Number_Of_Zones_In_Base

    !------------------------------!
    !   Browse through all zones   !
    !------------------------------!
    do zone_id = 1, n_zones

     ! Read zone_id information (-> n_nodes, n_cells)
      call Cgns_Mod_Read_Zone_Info

      ! Read type of zone_id (->structured or unstructured)
      call Cgns_Mod_Read_Zone_Type

      ! Read number of element sections (->n_sects)
      call Cgns_Mod_Read_Number_Of_Element_Sections

      !-----------------------------------------!
      !   Browse through all element sections   !
      !-----------------------------------------!
      do sect_id = 1, n_sects

        ! Read info for an element section (-> cell_type, first_cell, last_cell)
        call Cgns_Mod_Read_Section_Info

      end do ! elements sections

!     we need to compare n_bc from ZoneBC and though elements
!     if they are different -> stop program

    end do ! zones

  end do ! bases

  !-----------------------!
  !   First run: results  !
  !-----------------------!

  print *, "# First run finished!"
  print *, "# Found nodes: ",           n_nodes
  print *, "# Found cells: ",           n_cells
  print *, "# Found hex cells: ",       n_hexa
  print *, "# Found pyramids cells: ",  n_pyra
  print *, "# Found prism cells: ",     n_pris
  print *, "# Found tetra cells: ",     n_tetr
  print *, "# Found triangles faces: ", n_tria
  print *, "# Found quads faces: ",     n_quad
  if (n_quad + n_tria .eq. 0) then
    print *, "# No b.c. faces were found !"
    stop
  end if
  print *, "# Found b.c. faces: ",      n_quad + n_tria

  !-------------------------------!
  !   Allocate memory for arrays  !
  !-------------------------------!

  allocate( x_coord(1:n_nodes) )
  allocate( y_coord(1:n_nodes) )
  allocate( z_coord(1:n_nodes) )
  allocate( bc_mark(1:n_nodes) )

  allocate( hexa_connections(1:8, 1:n_hexa) )
  allocate( pyra_connections(1:5, 1:n_pyra) )
  allocate( pris_connections(1:6, 1:n_pris) )
  allocate( tetr_connections(1:4, 1:n_tetr) )
  allocate( quad_connections(1:4, 1:n_quad) )
  allocate( tria_connections(1:3, 1:n_tria) )

  !------------------------------------!
  !   Second run: read arrays from DB  !
  !------------------------------------!

  print *, "# Filling arrays.."

  !------------------------------!
  !   Browse through all bases   !
  !------------------------------!
  do base_id = 1, n_bases

    !------------------------------!
    !   Browse through all zones   !
    !------------------------------!
    do zone_id = 1, n_zones

      !---------------------------!
      !   Read coordinates block  !
      !---------------------------!

      ! Reads number of coordinates arrays from zone_id (->n_coords)
      call Cgns_Mod_Read_Number_Of_Coordinates

      ! Read x, y and z coordinates
      do coord_id = 1, n_coords

        ! Reads coord_name of coord_id(->coord_name)
        call Cgns_Mod_Print_Coordinate_Info

        ! Read grid coordinates (-> x_, y_, z_coord)
        call Cgns_Mod_Read_Coordinate_Array

      end do ! coordinates

      !---------------------!
      !   Read cells block  !
      !---------------------!

      ! Browse through all sections to read elements
      do sect_id = 1, n_sects

        ! Read element data (count HEXA_8/PYRA_5/PENTA_6/TETRA_4/QUAD_4/TRI_3)
        call Cgns_Mod_Read_Section_Connections

        ! Mark nodes with sect_id, if b.c.
        call Cgns_Mod_Mark_Bound_Cond

      end do ! elements sections

    end do ! zones

  end do ! bases

  ! use this to check bc_mark.
  !do c = 1, n_nodes
  !  if (bc_mark(c) .ne. 0) print *, "# b.c. mark =", bc_mark(c), " x=", x_coord(c), " y=", y_coord(c), " z=", z_coord(c)
  !end do

  grid % n_nodes     = n_nodes
  grid % n_cells     = n_cells
  grid % n_bnd_cells = n_tria + n_quad

  ! Allocate memory =--> carefull, there is no checking!
  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)
  call Grid_Mod_Allocate_Cells(grid, grid % n_cells, grid % n_bnd_cells)


  !--------------------------------------!
  !  Conversion to T-FlowS B.C. format   !
  !--------------------------------------!

  ! Reorder elements connectivity according to neu_mod

  ! HEXA_8 (-> f8n)
  do c = lbound(hexa_connections,dim=2), &
    ubound(hexa_connections,dim=2) - lbound(hexa_connections,dim=2) + 1

    i = hexa_connections(4, c)
    hexa_connections(4, c) = hexa_connections(3, c)
    hexa_connections(3, c) = i
    i = hexa_connections(8, c)
    hexa_connections(8, c) = hexa_connections(7, c)
    hexa_connections(7, c) = i
  end do

  ! PYRA_5 (-> f5n)
  do c =  lbound(pyra_connections,dim=2), &
    ubound(pyra_connections,dim=2) - lbound(pyra_connections,dim=2) + 1

    i = pyra_connections(4, c)
    pyra_connections(4, c) = pyra_connections(3, c)
    pyra_connections(3, c) = i
  end do

  ! PENTA_6 (-> f6n)
  ! no changes

  ! TETRA_4 (-> f4n)
  ! no changes

stop

end subroutine