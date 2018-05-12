!==============================================================================!
  subroutine Comm_Mod_Write_Bnd_Real(fh, buffer, disp)
!------------------------------------------------------------------------------!
!   Sequential version of writing a "distributed" boundary-cell-based array.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh         ! file handle
  real    :: buffer(:)
  integer :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Write "distributed" boundary cell data 
  do c = 1, nb_t
    write(9) buffer(c)
  end do

  disp = disp + nb_t * SIZE_REAL

  end subroutine