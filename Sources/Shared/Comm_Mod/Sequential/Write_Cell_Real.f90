!==============================================================================!
  subroutine Comm_Mod_Write_Cell_Real(fh, buffer, disp)
!------------------------------------------------------------------------------!
!   Sequential version of writing a "distributed" cell-based array.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh         ! file handle
  real    :: buffer(:)
  integer :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Write "distributed" cell data 
  do c = 1, nc_t
    write(9) buffer(c)
  end do

  disp = disp + nc_t * SIZE_REAL

  end subroutine
