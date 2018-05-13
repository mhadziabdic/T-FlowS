!==============================================================================!
  subroutine Save_Backup(grid, time_step, name_save)
!------------------------------------------------------------------------------!
!   Writes backup files. name.backup                                          !
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use Const_Mod
  use Flow_Mod
! use les_mod
  use Comm_Mod
! use rans_mod
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
  integer           :: fh, error, vs, d, s, c
!==============================================================================!

  store_name = problem_name

  problem_name = name_save

  ! Name backup file 
  call Name_File(0, name_out, '.backup')

  ! Open backup file
  call Comm_Mod_Open_File_Write(fh, name_out)

  ! Create new types
  call Comm_Mod_Create_New_Types()

  ! Initialize displacement
  d = 0

  !-----------------------------------------------------------------------!
  !   Save cell-centre coordinates.  Could be useful for interpolations   !
  !-----------------------------------------------------------------------! 
  call Write_Backup_3_Cell_Bnd(fh, d, 'coordinates',   & 
                               grid % xc(-nb_s:nc_s),  &
                               grid % yc(-nb_s:nc_s),  &
                               grid % zc(-nb_s:nc_s))
  !---------------!
  !               !
  !   Save data   !
  !               !
  !---------------! 

  ! Time
  call Write_Backup_1_Int(fh, d, 'time_step', time_step)

  !--------------!
  !   Velocity   !
  !--------------!
  call Write_Backup_3_Cell_Bnd(fh, d, 'vel',      u % n(-nb_s:nc_s),  &
                                                  v % n(-nb_s:nc_s),  &
                                                  w % n(-nb_s:nc_s))
  call Write_Backup_3_Cell(fh, d, 'vel_old',      u % o  (1:nc_s),    &
                                                  v % o  (1:nc_s),    &
                                                  w % o  (1:nc_s))
  call Write_Backup_3_Cell(fh, d, 'vel_adv',      u % a  (1:nc_s),    &
                                                  v % a  (1:nc_s),    &
                                                  w % a  (1:nc_s))
  call Write_Backup_3_Cell(fh, d, 'vel_adv_old',  u % a_o(1:nc_s),    &
                                                  v % a_o(1:nc_s),    &
                                                  w % a_o(1:nc_s))
  call Write_Backup_3_Cell(fh, d, 'vel_diff_old', u % d_o(1:nc_s),    &
                                                  v % d_o(1:nc_s),    &
                                                  w % d_o(1:nc_s))
  !--------------------------------------!
  !   Pressure and pressure correction   !
  !--------------------------------------!
  call Write_Backup_1_Cell_Bnd(fh, d, 'press',       p  % n(-nb_s:nc_s))
  call Write_Backup_1_Cell_Bnd(fh, d, 'press_corr',  pp % n(-nb_s:nc_s))

  !----------------------------------------------------!
  !   Mass flow rates (ask Egor if name is correct?)   !
  !----------------------------------------------------!
  call Write_Backup_1_Face(fh, d, 'mass_flow_rate', flux(1:nf_s+nbf_s))

  !--------------!
  !   Etnhalpy   !
  !--------------!
  if(heat_transfer == YES) then
    call Write_Backup_1_Cell_Bnd(fh, d, 'temp',      t  % n(-nb_s:nc_s))
    call Write_Backup_1_Cell_Bnd(fh, d, 'heat_flux', t  % q(-nb_s:nc_s))
    call Write_Backup_1_Cell(fh, d, 'temp_old',      t  % o  (1:nc_s))
    call Write_Backup_1_Cell(fh, d, 'temp_adv',      t  % a  (1:nc_s))
    call Write_Backup_1_Cell(fh, d, 'temp_adv_old',  t  % a_o(1:nc_s))
    call Write_Backup_1_Cell(fh, d, 'temp_diff_old', t  % d_o(1:nc_s))
    call Write_Backup_1_Cell(fh, d, 'temp_cros',     t  % c  (1:nc_s))
    call Write_Backup_1_Cell(fh, d, 'temp_cros_old', t  % c_o(1:nc_s))
  end if

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  problem_name = store_name

  end subroutine

!TEST  ! Test face buffers
!TEST  do s = 1, grid % n_faces
!TEST    flux(s) = 1000000.0 * grid % xf(s) +  &
!TEST                 1000.0 * grid % yf(s) +  &
!TEST                    1.0 * grid % zf(s)
!TEST    print '(a6,i4.4,a4,f18.3)', ' flux(', s, ') = ', flux(s)
!TEST  end do
