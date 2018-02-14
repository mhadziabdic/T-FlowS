!==============================================================================!
  subroutine Cgns_Mod_Write_Block_Info_Seq(base, block)
!------------------------------------------------------------------------------!
!   Gets n_bases from base node                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id             ! base index number    
  integer           :: block_id            ! block index number
  character(len=80) :: block_name          ! name of the block
  integer           :: block_mesh_info(3)  ! n_nodes, n_cells, and ...
                                           ! ... n_b_nodes(if sorted)
  integer           :: error
!==============================================================================!

  ! Set input parameters
  base_id         = base
  block_id        = block
  block_name      = trim(cgns_base(base) % block(block) % name)
  block_mesh_info = cgns_base(base) % block(block) % mesh_info

    ! Create and/or write to a zone node 
    call Cg_Zone_Write_F(file_id,         &
                         base_id,         &
                         block_name,      &
                         block_mesh_info, &
                         Unstructured,    &
                         block_id,        &
                         error)

  if (error .ne. 0) then
    print *, "# Failed to write block info"
    call Cg_Error_Exit_F()
  endif

  if(verbose) then
    print *, '#     =========================================='
    print *, "#     Block name:                ",  &
             cgns_base(base) % block(block) % name
    print *, '#     =========================================='
    print *, "#     Block index:               ", block    
    print *, "#     Nodes:                     ",  &
             cgns_base(base) % block(block) % mesh_info(1)
    print *, "#     Cells:                     ",  &
             cgns_base(base) % block(block) % mesh_info(2)
    print *, "#     Boundary nodes(if sorted): ",  &
             cgns_base(base) % block(block) % mesh_info(3)
  end if

  if (cgns_base(base) % block(block) % mesh_info(3) .ne. 0) then
    print *, "# Boundary condition nodes != 0 -> Unsupported"
    stop
  endif

  end subroutine
