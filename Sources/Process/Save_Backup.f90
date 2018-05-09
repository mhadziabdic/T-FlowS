!==============================================================================!
  subroutine Save_Backup(grid, time_step)
!------------------------------------------------------------------------------!
!   Writes backup files. name.backup                                          !
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
  integer           :: fh, error, vs, disp
!==============================================================================!

  ! Name backup file 
  call Name_File(0, name_out, '.backup')

  ! Open backup file
  call Comm_Mod_Open_File_Write(fh, name_out)

  ! Create new types
  call Comm_Mod_Create_New_Types()

  disp = 0

  !-----------------------------------------------------------------------!
  !   Save cell-centre coordinates.  Could be useful for interpolations   !
  !-----------------------------------------------------------------------! 
  vn = 'x_coordinate';     call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t + nb_t) * 8;  call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, grid % xc( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, grid % xc(-nb_s:-1), disp)

  vn = 'y_coordinate';     call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t + nb_t) * 8;  call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, grid % yc( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, grid % yc(-nb_s:-1), disp)

  vn = 'z_coordinate';     call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t + nb_t) * 8;  call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, grid % zc( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, grid % zc(-nb_s:-1), disp)

  !---------------!
  !               !
  !   Save data   !
  !               !
  !---------------! 

  ! Velocity 
  vn = 'u_velocity';      call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t + nb_t) * 8; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, u % n( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, u % n(-nb_s:-1), disp)

  vn = 'v_velocity';      call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t + nb_t) * 8; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, v % n( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, v % n(-nb_s:-1), disp)

  vn = 'w_velocity';      call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t + nb_t) * 8; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, w % n( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, w % n(-nb_s:-1), disp)


  vn = 'u_velocity_old';  call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;        call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, u % o( 1: nc_s), disp)

  vn = 'v_velocity_old';  call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;        call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, v % o( 1: nc_s), disp)

  vn = 'w_velocity_old';  call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;        call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, w % o( 1: nc_s), disp)


  vn = 'u_velocity_advection'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;             call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, u % a( 1: nc_s), disp)

  vn = 'v_velocity_advection'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;             call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, v % a( 1: nc_s), disp)

  vn = 'w_velocity_advection'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;             call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, w % a( 1: nc_s), disp)

  vn = 'u_velocity_advection_old'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;                 call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, u % a_o( 1: nc_s), disp)

  vn = 'v_velocity_advection_old'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;                 call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, v % a_o( 1: nc_s), disp)

  vn = 'w_velocity_advection_old'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;                 call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, w % a_o( 1: nc_s), disp)


  vn = 'u_velocity_diffusion_old'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;                 call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, u % d_o( 1: nc_s), disp)

  vn = 'v_velocity_diffusion_old'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;                 call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, v % d_o( 1: nc_s), disp)

  vn = 'w_velocity_diffusion_old'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t) * 8;                 call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, w % d_o( 1: nc_s), disp)


  ! Close backup file
  call Comm_Mod_Close_File(fh)

  end subroutine
