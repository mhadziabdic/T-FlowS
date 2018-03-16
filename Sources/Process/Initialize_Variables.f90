!==============================================================================!
  subroutine Initialize_Variables(grid)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use les_mod
  use Comm_Mod
  use rans_mod
  use Grid_Mod
  use Bulk_Mod
  use User_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------!
  integer :: Key_Ind
!-----------------------------------[Locals]-----------------------------------!
  integer           :: i, c, c1, c2, m, s, n, found, nks, nvs
  integer :: n_wall, n_inflow, n_outflow, n_symmetry, n_heated_wall, n_convect
  character(len=80) :: keys(128)
  real              :: vals(0:128) ! Note that they start from zero!
  real              :: s_tot

  ! Default values for initial conditions 
  real, parameter   :: u_def   = 0.0,  v_def   = 0.0,  w_def    = 0.0
  real, parameter   :: p_def   = 0.0,  t_def   = 0.0,  q_def    = 0.0
  real, parameter   :: kin_def = 0.0,  eps_def = 0.0,  f22_def  = 0.0
  real, parameter   :: vis_def = 0.0,  v2_def  = 0.0,  zeta_def = 0.0
  real, parameter   :: uu_def  = 0.0,  vv_def  = 0.0,  ww_def   = 0.0
  real, parameter   :: uv_def  = 0.0,  uw_def  = 0.0,  vw_def   = 0.0
!==============================================================================!

  area  = 0.0
  print *, 'grid % n_materials: ', grid % n_materials

  ! Found the line where boundary condition defintion is defined
  call Control_Mod_Position_At_One_Key('INITIAL_CONDITION',       &
                                       found,                     &
                                       .true.)

  ! Found the section with intial condions
  if(found == YES) then

print *, 'FOUND INITIAL CONDITIONS'
    call Control_Mod_Read_Strings_On('VARIABLES', keys,    nks, .true.)
    call Control_Mod_Read_Real_Array_On('VALUES', vals(1), nvs, .true.)

    ! Check validity of the input
    if(nks .eq. 0 .or. nvs .eq. 0) then
      print '(2a)', '# Critical, for initial condition: ',        &
                    ' no values or variables have been provided' 
      stop
    end if
    if(nks .ne. nvs) then
      print '(2a)', '# Critical for initial conditions, number of values ',  &
                    ' is not the same as number of provided variable names' 
      stop
    end if
 
    ! Input is valid, turn keys to upper case
    do i = 1, nks
      call To_Upper_Case(keys(i))
    end do

    do n = 1, grid % n_materials
      do c = 1, grid % n_cells

        u % mean(c) = 0.0
        v % mean(c) = 0.0
        w % mean(c) = 0.0

        vals(0) = u_def;  u % n(c) = vals(Key_Ind('U', keys, nks))
        vals(0) = v_def;  v % n(c) = vals(Key_Ind('V', keys, nks))
        vals(0) = w_def;  w % n(c) = vals(Key_Ind('W', keys, nks))

        u % o(c)  = u % n(c)
        u % oo(c) = u % n(c)
        v % o(c)  = v % n(c)
        v % oo(c) = v % n(c)
        w % o(c)  = w % n(c)
        w % oo(c) = w % n(c)

        if(heat_transfer == YES) then
          vals(0) = t_def;  t % n(c) = vals(Key_Ind('T', keys, nks))
          t % o(c)  = t % n(c)
          t % oo(c) = t % n(c)
          Tinf      = t % n(c)
        end if 

        if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
           turbulence_model == HANJALIC_JAKIRLIC) then
          vals(0) = uu_def;  uu  % n(c) = vals(Key_Ind('UU',  keys, nks))
          vals(0) = vv_def;  vv  % n(c) = vals(Key_Ind('VV',  keys, nks))
          vals(0) = ww_def;  ww  % n(c) = vals(Key_Ind('WW',  keys, nks))
          vals(0) = uv_def;  uv  % n(c) = vals(Key_Ind('UV',  keys, nks))
          vals(0) = uw_def;  uw  % n(c) = vals(Key_Ind('UW',  keys, nks))
          vals(0) = vw_def;  vw  % n(c) = vals(Key_Ind('VW',  keys, nks))
          vals(0) = eps_def; eps % n(c) = vals(Key_Ind('EPS', keys, nks))
          uu % o(c)  = uu % n(c)
          uu % oo(c) = uu % n(c)
          vv % o(c)  = vv % n(c)
          vv % oo(c) = vv % n(c)
          ww % o(c)  = ww % n(c)
          ww % oo(c) = ww % n(c)
          uv % o(c)  = uv % n(c)
          uv % oo(c) = uv % n(c)
          uw % o(c)  = uw % n(c)
          uw % oo(c) = uw % n(c)
          vw % o(c)  = vw % n(c)
          vw % oo(c) = vw % n(c)
          if(turbulence_model == REYNOLDS_STRESS_MODEL) then
            vals(0) = f22_def; f22 % n(c) = vals(Key_Ind('F22', keys, nks))
            f22 % o(c)  = f22 % n(c)
            f22 % oo(c) = f22 % n(c)
          end if
        end if
  
        if(turbulence_model == K_EPS .or.  &
           turbulence_model == HYBRID_PITM) then
          vals(0) = kin_def; kin % n(c) = vals(Key_Ind('KIN', keys, nks))
          vals(0) = eps_def; eps % n(c) = vals(Key_Ind('EPS', keys, nks))
          kin % o(c)  = kin % n(c)
          kin % oo(c) = kin % n(c)
          eps % o(c)  = eps % n(c)
          eps % oo(c) = eps % n(c)
          u_tau(c)  = 0.047
          y_plus(c) = 30.0
        end if
  
        if(turbulence_model == K_EPS_ZETA_F  .or.  & 
           turbulence_model == HYBRID_K_EPS_ZETA_F) then
          vals(0) = kin_def;  kin  % n(c) = vals(Key_Ind('KIN',  keys, nks))
          vals(0) = eps_def;  eps  % n(c) = vals(Key_Ind('EPS',  keys, nks))
          vals(0) = zeta_def; zeta % n(c) = vals(Key_Ind('ZETA', keys, nks))
          vals(0) = f22_def;  f22  % n(c) = vals(Key_Ind('F22',  keys, nks))
          kin  % o(c)  = kin  % n(c)
          kin  % oo(c) = kin  % n(c)
          eps  % o(c)  = eps  % n(c)
          eps  % oo(c) = eps  % n(c)
          zeta % o(c)  = zeta % n(c)
          zeta % oo(c) = zeta % n(c)
          f22  % o(c)  = f22  % n(c)
          f22  % oo(c) = f22  % n(c)
          u_tau(c)  = 0.047
          y_plus(c) = 30.0
        end if
  
        if(turbulence_model == SPALART_ALLMARAS .or.  &
           turbulence_model == DES_SPALART) then      
          vals(0) = vis_def; vis % n(c) = vals(Key_Ind('VIS', keys, nks))
          vis % o(c)  = vis % n(c)
          vis % oo(c) = vis % n(c)
        end if

      end do   ! through cells
    end do   !end do n=1,grid % n_materials

  end if

  call User_Mod_Initialize(grid)

!@if(TGV == YES) then
!@  do c = 1, grid % n_cells
!@    u % n(c)  = -sin(grid % xc(c))*cos(grid % yc(c))
!@    u % o(c)  = -sin(grid % xc(c))*cos(grid % yc(c))
!@    u % oo(c) = -sin(grid % xc(c))*cos(grid % yc(c))
!@    v % n(c)  =  cos(grid % xc(c))*sin(grid % yc(c))
!@    v % o(c)  =  cos(grid % xc(c))*sin(grid % yc(c))
!@    v % oo(c) =  cos(grid % xc(c))*sin(grid % yc(c))
!@    w % n(c)  = 0.0
!@    w % o(c)  = 0.0
!@    w % oo(c) = 0.0
!@    P % n(c)  = 0.25*(cos(2*grid % xc(c)) + cos(2*grid % yc(c)))
!@  end do
!@end if

  !---------------------------------!
  !      Calculate the inflow       !
  !   and initializes the flux(s)   ! 
  !   at both inflow and outflow    !
  !---------------------------------!
  n_wall        = 0
  n_inflow      = 0
  n_outflow     = 0
  n_symmetry    = 0
  n_heated_wall = 0
  n_convect     = 0
  do m = 1, grid % n_materials
    bulk(m) % mass_in = 0.0
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2  < 0) then 
        flux(s) = density*( u % n(c2) * grid % sx(s) + &
                            v % n(c2) * grid % sy(s) + &
                            w % n(c2) * grid % sz(s) )
                                       
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == INFLOW) then
          if(grid % material(c1) == m) then
            bulk(m) % mass_in = bulk(m) % mass_in - flux(s) 
          end if
          s_tot = sqrt(  grid % sx(s)**2  &
                       + grid % sy(s)**2  &
                       + grid % sz(s)**2)
          area = area  + s_tot
        endif
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALL)      &
          n_wall        = n_wall        + 1 
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == INFLOW)    &
          n_inflow      = n_inflow      + 1  
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW)   &
          n_outflow     = n_outflow     + 1 
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == SYMMETRY)  &
          n_symmetry    = n_symmetry    + 1 
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == WALLFL)    &
          n_heated_wall = n_heated_wall + 1 
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT)   &
          n_convect     = n_convect     + 1 
      else
        flux(s) = 0.0 
      end if
    end do
    call Comm_Mod_Global_Sum_Int(n_wall)
    call Comm_Mod_Global_Sum_Int(n_inflow)
    call Comm_Mod_Global_Sum_Int(n_outflow)
    call Comm_Mod_Global_Sum_Int(n_symmetry)
    call Comm_Mod_Global_Sum_Int(n_heated_wall)
    call Comm_Mod_Global_Sum_Int(n_convect)
    call Comm_Mod_Global_Sum_Real(bulk(m) % mass_in)
    call Comm_Mod_Global_Sum_Real(area)
  end do                  

  !----------------------!
  !   Initializes time   ! 
  !----------------------!
  if(this_proc  < 2) then
    if(n_inflow .gt. 0) then
      print *, '# Mass inflow =', bulk(1) % mass_in
      print *, '# Average inflow velocity =', bulk(1) % mass_in / area
    end if
    print *, '# number of faces on the wall        : ', n_wall
    print *, '# number of inflow faces             : ', n_inflow
    print *, '# number of outflow faces            : ', n_outflow
    print *, '# number of symetry faces            : ', n_symmetry
    print *, '# number of faces on the heated wall : ', n_heated_wall
    print *, '# number of convective outflow faces : ', n_convect
    print *, '# Variables initialized !'
  end if

  end subroutine
