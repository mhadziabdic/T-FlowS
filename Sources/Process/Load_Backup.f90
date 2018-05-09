!==============================================================================!
  subroutine Load_Backup(grid, time_step)
!------------------------------------------------------------------------------!
!   Read backup files. name.backup                                             !
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
! use Const_Mod
  use Flow_Mod
! use les_mod
  use Comm_Mod
! use rans_mod
! use Tokenizer_Mod
  use Grid_Mod
! use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: time_step
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_out, vn
  integer           :: fh, error, vs, disp, c
!==============================================================================!

  !------------------------!
  !   Create backup file   !
  !------------------------!
  call Name_File(0, name_out, '.backup')

  ! Open backup file
  call Comm_Mod_Open_File_Read(fh, name_out)

  ! Create new types
  call Comm_Mod_Create_New_Types()

  disp = 0

  !-----------------------------------------------!
  !   Skip three coordinates for the time being   !
  !-----------------------------------------------! 
  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  if(vn == 'x_coordinate') disp = disp + vs
  
  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  if(vn == 'y_coordinate') disp = disp + vs

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  if(vn == 'z_coordinate') disp = disp + vs

  !---------------!
  !               !
  !   Load data   !
  !               !
  !---------------! 

  ! Velocity 
  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, u % n( 1: nc_s), disp)
  call Comm_Mod_Read_Bnd_Real (fh, u % n(-nb_s:-1), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, v % n( 1: nc_s), disp)
  call Comm_Mod_Read_Bnd_Real (fh, v % n(-nb_s:-1), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, w % n( 1: nc_s), disp)
  call Comm_Mod_Read_Bnd_Real (fh, w % n(-nb_s:-1), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, u % o( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, v % o( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, w % o( 1: nc_s), disp)


  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, u % a( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, v % a( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, w % a( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, u % a_o( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, v % a_o( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, w % a_o( 1: nc_s), disp)


  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, u % d_o( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, v % d_o( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, w % d_o( 1: nc_s), disp)

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  end subroutine
