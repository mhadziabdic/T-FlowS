!==============================================================================!
  subroutine Cgns_Mod_Read_Number_Of_Coordinates_In_Block(base, block)
!------------------------------------------------------------------------------!
!   Reads number of coordinates arrays from block_id
!------------------------------------------------------------------------------!
  implicit none
  integer*8 :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer*8 :: ier
!==============================================================================!

  ! Get number of coordinate arrays (1 for unstructure, 3 for structured)
  call Cg_Ncoords_F(file_id,   &  !  cgns file index number
                    base,      &  !  base index number
                    block,     &  !  block index number
                    n_coords,  &  !  number of coordinate arrays for block
                    ier)          !  error status

  if (ier .ne. 0) then
    print *, "# Failed to get number of coordinate arrays"
    call Cg_Error_Exit_F()
  endif

  print *, "# Number of coordinate arrays for block:", n_coords

  end subroutine
