!==============================================================================!
  subroutine Cgns_Mod_Read_Bnd_Conds_Data(base, block, bc)
!------------------------------------------------------------------------------!
!   Reads boundary condition info.                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer*8 :: base, block, bc
!-----------------------------------[Locals]-----------------------------------!
  integer*8 :: ier
!==============================================================================!
!   Description of arguments for CGNS function call
!
!   Cg_Boco_Read_F(file_id,         &  ! cgns file index number
!                  base,            &  ! base index number
!                  block,           &  ! block index number
!                  bnd_cond,        &  ! boundary condition index number
!                  buffer_r1,       &  ! array of point or element indices 
!                                      ! ... defining bnd. cond. region
!                  NormalListFlag,  &  ! list of vectors normal to the 
!                                      ! ... boundary condition patch 
!                                      ! ... pointing into the interior 
!                                      ! ... of the block
!                  ier)                ! error status
!------------------------------------------------------------------------------!

  ! Read boundary condition data and normals
  call Cg_Boco_Read_F(file_id,                 &
                      base,                    &
                      block,                   &
                      cgns_bnd_cond(bc) % id,  &
                      buffer_r1,               &
                      NormalListFlag,          &
                      ier)
  print *, "sdfsdf"

  print *, "# ---------------------------"

  end subroutine
