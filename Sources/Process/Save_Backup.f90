!==============================================================================!
  subroutine Save_Backup(grid, time_step, name_save)
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
  type(Grid_Type)  :: grid
  integer          :: time_step
  character(len=*) :: name_save
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_out, store_name, vn
  integer           :: fh, error, vs, disp, s, c1, c2, c
!==============================================================================!

  store_name = problem_name

  problem_name = name_save

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
  vn = 'x_coordinate';            call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t + nb_t) * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, grid % xc( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, grid % xc(-nb_s:-1), disp)

  vn = 'y_coordinate';            call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t + nb_t) * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, grid % yc( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, grid % yc(-nb_s:-1), disp)

  vn = 'z_coordinate';            call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t + nb_t) * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, grid % zc( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, grid % zc(-nb_s:-1), disp)

  !---------------!
  !               !
  !   Save data   !
  !               !
  !---------------! 

  ! Time
  vn = 'time_step'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = SIZE_INT;    call Comm_Mod_Write_Int (fh, vs, disp)
  vs = time_step;   call Comm_Mod_Write_Int (fh, vs, disp)

  ! Velocity 
  vn = 'velocity';                    call Comm_Mod_Write_Text(fh, vn, disp)
  vs = 3 * (nc_t + nb_t) * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, u % n( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, v % n( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, w % n( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, u % n(-nb_s:-1), disp)
  call Comm_Mod_Write_Bnd_Real (fh, v % n(-nb_s:-1), disp)
  call Comm_Mod_Write_Bnd_Real (fh, w % n(-nb_s:-1), disp)

  vn = 'velocity_old';       call Comm_Mod_Write_Text(fh, vn, disp)
  vs = 3 * nc_t * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, u % o( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, v % o( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, w % o( 1: nc_s), disp)

  vn = 'velocity_advection'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = 3 * nc_t * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, u % a( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, v % a( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, w % a( 1: nc_s), disp)

  vn = 'velocity_advection_old'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = 3 * nc_t * SIZE_REAL;     call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, u % a_o( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, v % a_o( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, w % a_o( 1: nc_s), disp)

  vn = 'velocity_diffusion_old'; call Comm_Mod_Write_Text(fh, vn, disp)
  vs = 3 * nc_t * SIZE_REAL;     call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, u % d_o( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, v % d_o( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, w % d_o( 1: nc_s), disp)

  ! Pressure
  vn = 'pressure';                    call Comm_Mod_Write_Text(fh, vn, disp)
  vs = 2 * (nc_t + nb_t) * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Cell_Real(fh, p  % n( 1: nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, pp % n( 1: nc_s), disp)
  call Comm_Mod_Write_Bnd_Real (fh, p  % n(-nb_s:-1), disp)
  call Comm_Mod_Write_Bnd_Real (fh, pp % n(-nb_s:-1), disp)

  !----------------------------------------------------!
  !   Mass flow rates (ask Egor if name is correct?)   !
  !----------------------------------------------------!

  ! Change the sign of flux where necessary
  do s = nf_s + 1, nf_s + nbf_s
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    flux(s) = flux(s) * buf_face_sgn(s-nf_s)
  end do

!TEST  ! Test face buffers
!TEST  do s = 1, grid % n_faces
!TEST    flux(s) = 1000000.0 * grid % xf(s) +  &
!TEST                 1000.0 * grid % yf(s) +  &
!TEST                    1.0 * grid % zf(s)
!TEST    print '(a6,i4.4,a4,f18.3)', ' flux(', s, ') = ', flux(s)
!TEST  end do

  ! Perform the actual saving
  vn = 'mass_flow_rate';   call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nf_t) * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)
  call Comm_Mod_Write_Face_Real(fh, flux( 1: nf_s + nbf_s), disp)

  ! Fix signs for the face fluxes back
  do s = nf_s + 1, nf_s + nbf_s
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    flux(s) = flux(s) * buf_face_sgn(s-nf_s)
  end do

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  problem_name = store_name

  end subroutine
