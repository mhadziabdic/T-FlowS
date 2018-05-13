!==============================================================================!
  subroutine Load_Backup(grid, time_step, restart)
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
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: time_step
  logical         :: restart 
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_in, vn, answer
  integer           :: fh, error, vs, disp, c, c1, c2, s
!==============================================================================!

  ! Full name is specified in control file
  call Control_Mod_Load_Backup_Name(name_in)

  answer = name_in
  call To_Upper_Case(answer)
  if(answer == 'SKIP') then
    restart = .false.
    return 
  end if

  ! Open backup file
  call Comm_Mod_Open_File_Read(fh, name_in)

  ! Create new types
  call Comm_Mod_Create_New_Types()

  disp = 0

  !-----------------------------------------------!
  !   Skip three coordinates for the time being   !
  !-----------------------------------------------! 
  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn)
  if(vn == 'x_coordinate') disp = disp + vs
  
  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn)
  if(vn == 'y_coordinate') disp = disp + vs

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn)
  if(vn == 'z_coordinate') disp = disp + vs

  !---------------!
  !               !
  !   Load data   !
  !               !
  !---------------! 

  ! Time step
  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn), 'with size', vs
  call Comm_Mod_Read_Int (fh, time_step, disp)

  ! Velocity 
  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, u % n( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, v % n( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, w % n( 1: nc_s), disp)
  call Comm_Mod_Read_Bnd_Real (fh, u % n(-nb_s:-1), disp)
  call Comm_Mod_Read_Bnd_Real (fh, v % n(-nb_s:-1), disp)
  call Comm_Mod_Read_Bnd_Real (fh, w % n(-nb_s:-1), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, u % o( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, v % o( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, w % o( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, u % a( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, v % a( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, w % a( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, u % a_o( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, v % a_o( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, w % a_o( 1: nc_s), disp)

  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, u % d_o( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, v % d_o( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, w % d_o( 1: nc_s), disp)

  ! Pressure
  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Cell_Real(fh, p  % n( 1: nc_s), disp)
  call Comm_Mod_Read_Cell_Real(fh, pp % n( 1: nc_s), disp)
  call Comm_Mod_Read_Bnd_Real (fh, p  % n(-nb_s:-1), disp)
  call Comm_Mod_Read_Bnd_Real (fh, pp % n(-nb_s:-1), disp)

  !----------------------------------------------------!
  !   Mass flow rates (ask Egor if name is correct?)   !
  !----------------------------------------------------!

  ! Perform actual reading
  call Comm_Mod_Read_Text(fh, vn, disp)
  call Comm_Mod_Read_Int (fh, vs, disp)
  if(this_proc < 2) print *, '# Found variable: ', trim(vn)
  call Comm_Mod_Read_Face_Real(fh, flux( 1: nf_s + nbf_s), disp)

  ! Fix signs for the face fluxes back
  do s = nf_s + 1, nf_s + nbf_s
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    flux(s) = flux(s) * buf_face_sgn(s-nf_s)
  end do

!TEST  do s = 1, nf_s
!TEST    print '(a6,i4.4,a4,f18.3)', ' flux(', face_map(s)+1, ') = ', flux(s)
!TEST  end do
!TEST  do s = 1, nbf_s
!TEST    print '(a6,i4.4,a4,f18.3)', ' flux(', buf_face_map(buf_face_ord(s))+1, ') = ', flux(nf_s+s)
!TEST  end do

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  end subroutine
