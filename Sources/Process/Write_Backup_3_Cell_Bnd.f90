!==============================================================================!
  subroutine Write_Backup_3_Cell_Bnd(fh, disp, var_name, com1, com2, com3)
!------------------------------------------------------------------------------!
!   Writes a vector variable with boundary cells to backup file.               !
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: fh, disp
  character(len=*) :: var_name
  real             :: com1(-nb_s:nc_s), com2(-nb_s:nc_s), com3(-nb_s:nc_s)
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: vn
  integer           :: vs  ! variable size
!==============================================================================!

  ! Vector with boundaries
  vn = var_name;                      call Comm_Mod_Write_Text(fh, vn, disp)
  vs = 3 * (nc_t + nb_t) * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Cell_Real(fh, com1(1:nc_s),   disp)
  call Comm_Mod_Write_Cell_Real(fh, com2(1:nc_s),   disp)
  call Comm_Mod_Write_Cell_Real(fh, com3(1:nc_s),   disp)
  call Comm_Mod_Write_Bnd_Real (fh, com1(-nb_s:-1), disp)
  call Comm_Mod_Write_Bnd_Real (fh, com2(-nb_s:-1), disp)
  call Comm_Mod_Write_Bnd_Real (fh, com3(-nb_s:-1), disp)

  end subroutine
