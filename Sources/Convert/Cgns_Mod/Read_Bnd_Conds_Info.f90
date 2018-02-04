!==============================================================================!
  subroutine Cgns_Mod_Read_Bnd_Conds_Info(base, block, bc)
!------------------------------------------------------------------------------!
!   Reads boundary condition info.                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer*8 :: base, block, bc
!-----------------------------------[Locals]-----------------------------------!
  integer*8            :: base_id         ! base index number
  integer*8            :: block_id        ! block index number
  integer*8            :: bc_id           ! block index number
  character(len=80)    :: bc_name         ! name of the boundary condition
  integer*8            :: bc_type         ! boundary condition type
  integer*8            :: bc_ptset_type   ! boundary node/cell placement
  integer*8            :: bc_n_nodes      ! boundary nodes or cells
  integer*8            :: NormalIndex(3)  
  integer*8            :: NormalListFlag
  integer*8            :: bc_data_type     
  integer*8            :: bc_n_datasets
  integer*8            :: error
  integer*8, parameter :: one = 1         ! go figure :-(
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  bc_id    = bc       

  ! Position yourself at ZoneBC
  call Cg_Goto_F(file_id,     &  ! cgns file index number
                 base,        &  ! base index number
                 error,       &  ! error status
                 'Zone_t',    &  ! node of block type
                 block,       &  ! block index number
                 'ZoneBC_t',  &  ! search for node "ZoneBC_t"
                 one,         &  ! ???
                 'end')          ! indicates end of call
  if (error.ne.0) then
    print *,"# Failed to navigate ro ZoneBC node"
    call Cg_Error_Exit_F()
  endif

  ! Get boundary condition info
  call Cg_Boco_Info_F(file_id,         &
                      base_id,         &
                      block_id,        &
                      bc_id,           &
                      bc_name,         &
                      bc_type,         &
                      bc_ptset_type,   &
                      bc_n_nodes,      &
                      NormalIndex,     &
                      NormalListFlag,  &
                      bc_data_type,    &
                      bc_n_datasets,   &
                      error)             
  if (error .ne. 0) then
    print *,"# Failed to read boundary conditions info"
    call Cg_Error_Exit_F()
  endif
 
  ! Fetch received parameters
  cgns_base(base) % block(block) % bnd_cond(bc) % name    = trim(bc_name)
  cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes = bc_n_nodes
  cgns_base(base) % block(block) % bnd_cond(bc) % mark    = bc     

  if(verbose) then
    print *, '#       ----------------------------------------'
    print *, "#       Boundary condition name:   ",   &
             trim(cgns_base(base) % block(block) % bnd_cond(bc) % name)
    print *, '#       ----------------------------------------'
    print *, "#       Boundary condition index:  ", bc
    print *, "#       Boundary condition nodes:  ",   &
             cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes
    print *, "#       Boundary condition Extent: ",   &
             PointSetTypeName(bc_ptset_type)
  end if

  end subroutine
