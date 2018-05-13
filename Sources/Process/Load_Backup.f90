!==============================================================================!
  subroutine Load_Backup(grid, time_step, restart)
!------------------------------------------------------------------------------!
!   Read backup files. name.backup                                             !
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use Const_Mod
  use Flow_Mod
! use les_mod
  use Comm_Mod
! use rans_mod
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
  integer           :: fh, error, vs, d, c, s
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

  ! Initialize displacement
  d = 0

  !-----------------------------------------------!
  !   Skip three coordinates for the time being   !
  !-----------------------------------------------! 
  call Read_Backup_3_Cell_Bnd(fh, d, 'coordinates',  &
                              grid % xc(-nb_s:nc_s),  &
                              grid % yc(-nb_s:nc_s),  &
                              grid % zc(-nb_s:nc_s))

  !---------------!
  !               !
  !   Load data   !
  !               !
  !---------------! 

  ! Time step
  call Read_Backup_1_Int(fh, d, 'time_step', time_step)

  !--------------!
  !   Velocity   !
  !--------------!
  call Read_Backup_3_Cell_Bnd(fh, d, 'vel',      u % n(-nb_s:nc_s),  &
                                                 v % n(-nb_s:nc_s),  &
                                                 w % n(-nb_s:nc_s))
  call Read_Backup_3_Cell(fh, d, 'vel_old',      u % o  (1:nc_s),    &
                                                 v % o  (1:nc_s),    &
                                                 w % o  (1:nc_s))
  call Read_Backup_3_Cell(fh, d, 'vel_adv',      u % a  (1:nc_s),    &
                                                 v % a  (1:nc_s),    &
                                                 w % a  (1:nc_s))
  call Read_Backup_3_Cell(fh, d, 'vel_adv_old',  u % a_o(1:nc_s),    &
                                                 v % a_o(1:nc_s),    &
                                                 w % a_o(1:nc_s))
  call Read_Backup_3_Cell(fh, d, 'vel_diff_old', u % d_o(1:nc_s),    &
                                                 v % d_o(1:nc_s),    &
                                                 w % d_o(1:nc_s))

  !--------------------------------------!
  !   Pressure and pressure correction   !
  !--------------------------------------!
  call Read_Backup_1_Cell_Bnd(fh, d, 'press',       p  % n(-nb_s:nc_s))
  call Read_Backup_1_Cell_Bnd(fh, d, 'press_corr',  pp % n(-nb_s:nc_s))

  !----------------------------------------------------!
  !   Mass flow rates (ask Egor if name is correct?)   !
  !----------------------------------------------------!
  call Read_Backup_1_Face(fh, d, 'mass_flow_rate', flux(1:nf_s+nbf_s))

  !--------------!
  !   Etnhalpy   !
  !--------------!
  if(heat_transfer == YES) then
    call Read_Backup_1_Cell_Bnd(fh, d, 'temp',      t  % n(-nb_s:nc_s))
    call Read_Backup_1_Cell_Bnd(fh, d, 'heat_flux', t  % q(-nb_s:nc_s))
    call Read_Backup_1_Cell(fh, d, 'temp_old',      t  % o  (1:nc_s))
    call Read_Backup_1_Cell(fh, d, 'temp_adv',      t  % a  (1:nc_s))
    call Read_Backup_1_Cell(fh, d, 'temp_adv_old',  t  % a_o(1:nc_s))
    call Read_Backup_1_Cell(fh, d, 'temp_diff_old', t  % d_o(1:nc_s))
    call Read_Backup_1_Cell(fh, d, 'temp_cros',     t  % c  (1:nc_s))
    call Read_Backup_1_Cell(fh, d, 'temp_cros_old', t  % c_o(1:nc_s))
  end if

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  end subroutine

!TEST  do s = 1, nf_s
!TEST    print '(a6,i4.4,a4,f18.3)', ' flux(', face_map(s)+1, ') = ', flux(s)
!TEST  end do
!TEST  do s = 1, nbf_s
!TEST    print '(a6,i4.4,a4,f18.3)', ' flux(', buf_face_map(buf_face_ord(s))+1, ') = ', flux(nf_s+s)
!TEST  end do

