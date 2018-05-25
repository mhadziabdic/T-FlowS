!==============================================================================!
  subroutine Load_Boundary_Conditions(grid, in_out)
!------------------------------------------------------------------------------!
!   Reads: .bnd file                                                           !
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use rans_mod
  use Comm_Mod, only: this_proc
  use Tokenizer_Mod
  use Grid_Mod
  use User_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: in_out
!----------------------------------[Calling]-----------------------------------!
  real    :: Distance
  integer :: Key_Ind
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, m, l, k, i, n, n_points, nks, nvs, found, us
! integer           :: grid_bnd_type
  character(len=80) :: name_prof(128)
  real              :: wi, dist_min, x, y, z, xp, dist
  real, allocatable :: prof(:,:)
  logical           :: here
  character(len=80) :: bc_type_name, dumms
  integer           :: bc_type_tag
  character(len=80) :: keys(128)
  real              :: vals(0:128)           ! note that they start from zero!
  integer           :: types_per_color(128)  ! how many types in each color
  character(len=80) :: types_names(128)      ! name of each type
  logical           :: types_file(128)       ! type specified in a file?
  integer           :: c_types               ! counter types
  character(len=4)  :: q_name = 'Q_00'

  ! Default values for boundary conditions
  real, parameter :: u_def   = 0.0,  v_def    = 0.0,  w_def   = 0.0
  real, parameter :: p_def   = 0.0,  t_def    = 0.0,  q_def   = 0.0
  real, parameter :: kin_def = 0.0,  eps_def  = 0.0,  f22_def = 0.0
  real, parameter :: vis_def = 0.0,  zeta_def = 0.0
  real, parameter :: uu_def  = 0.0,  vv_def   = 0.0,  ww_def  = 0.0
  real, parameter :: uv_def  = 0.0,  uw_def   = 0.0,  vw_def  = 0.0
  real, parameter :: c_def   = 0.0
!==============================================================================!

  !-----------------------------------!
  !   Add one type for buffer cells   !
  !-----------------------------------!
  n = grid % n_bnd_cond
  grid % bnd_cond % type(n+1) = BUFFER
  do c = -1,-grid % n_bnd_cells,-1
    if(grid % bnd_cond % color(c) == BUFFER) then
      grid % bnd_cond % color(c) = n+1
    end if
  end do

  !----------------------------------------------------------------!
  !   Count number of types per boundary condition, total number   !
  !        of types specified, and also extract their names        !
  !----------------------------------------------------------------!
  types_per_color(:) = 0
  types_file(:)      = .false.
  c_types            = 0 

  do n = 1, grid % n_bnd_cond
    call Control_Mod_Position_At_Two_Keys('BOUNDARY_CONDITION',       &
                                          grid % bnd_cond % name(n),  &
                                          found,                      &
                                          .false.)
    if(found == YES) then 
1     continue

      ! Try to read next 'TYPE' in the control file
      call Control_Mod_Read_Char_Item_On('TYPE', 'VOID', bc_type_name, .false.) 

      ! Get out of the loop if you fail
      if(bc_type_name .eq. 'VOID') goto 2

      ! Skip following two lines
      call Control_Mod_Read_Char_Item_On('VARIABLES', 'VOID', dumms, .false.) 
      call Control_Mod_Read_Char_Item_On('VALUES',    'VOID', dumms, .false.) 

      types_per_color(n) = types_per_color(n) + 1
      c_types = c_types + 1
      types_names(c_types) = bc_type_name

      ! If dumms is 'VOID', it didn't find 'VALUES' 
      ! meaning that the keyword 'FILE' was specified
      if(dumms .eq. 'VOID') then
        types_file(c_types) = .true.
      end if

      goto 1
    end if

2 continue 

  end do

  c_types = 0 
  do n = 1, grid % n_bnd_cond
    do l = 1, types_per_color(n) 
      c_types = c_types + 1
      print *, 'TYPE NAME: ', types_names(c_types)
      if( types_file(c_types) ) then
        print *, 'SPECIFIED IN A FILE'
      else
        print *, 'SPECIFIED BY VALUES'
      end if
    end do
  end do

  !------------------------------------------------!
  !                                                !
  !                                                !
  !   Read boundary conditions from control file   !
  !                                                !
  !                                                !
  !------------------------------------------------!
  c_types = 0 

  do n = 1, grid % n_bnd_cond

    ! Position yourself well
    call Control_Mod_Position_At_Two_Keys('BOUNDARY_CONDITION',       &
                                          grid % bnd_cond % name(n),  &
                                          found,                      &
                                          .false.)
    do l = 1, types_per_color(n) 

      ! Update the counter
      c_types = c_types + 1

      !---------------------------------------------!
      !                                             !
      !   Read first line which is common for all   !
      !                                             !
      !---------------------------------------------!
      call Control_Mod_Read_Char_Item_On ('TYPE', 'WALL',  bc_type_name, .false.) 
      call To_Upper_Case(bc_type_name)

      ! Copy boundary conditions which were given for the grid
      if( bc_type_name == 'INFLOW') then 
        bc_type_tag = INFLOW
        grid % bnd_cond % type(n) = INFLOW
        PER_BC = NO
      else if( bc_type_name == 'WALL') then 
        bc_type_tag = WALL
        grid % bnd_cond % type(n) = WALL
      else if( bc_type_name == 'OUTFLOW') then 
        bc_type_tag = OUTFLOW
        grid % bnd_cond % type(n) = OUTFLOW
      else if( bc_type_name == 'SYMMETRY') then 
        bc_type_tag = SYMMETRY
        grid % bnd_cond % type(n) = SYMMETRY
      else if( bc_type_name == 'HEAT_FLUX') then 
        bc_type_tag = WALLFL  
        grid % bnd_cond % type(n) = WALLFL
      else if( bc_type_name == 'CONVECTIVE') then 
        bc_type_tag = CONVECT
        grid % bnd_cond % type(n) = CONVECT
      else if( bc_type_name == 'PRESSURE') then 
        bc_type_tag = PRESSURE  
        grid % bnd_cond % type(n) = PRESSURE
      else
        if(this_proc < 2)  &
          print *, '# Load_Boundary_Conditions: '//          &
                   '# Unknown boundary condition type: ',  &
                   bc_type_name
        stop  
      end if

      !----------------------------------------------!
      !                                              !
      !   Read second line which is common for all   !
      !                                              !
      !----------------------------------------------!
      call Control_Mod_Read_Strings_On('VARIABLES', keys, nks, .false.) 
      do i = 1, nks
        call To_Upper_Case(keys(i))
      end do

      !---------------------------------------------!
      !                                             !
      !   Boundary values are specified in a list   !
      !                                             !
      !---------------------------------------------!
      if( .not. types_file(c_types) ) then
        call Control_Mod_Read_Real_Array_On('VALUES', vals(1), nvs, .false.) 

        !--------------------------------------------------!
        !   Distribute boundary values to boundary cells   !
        !--------------------------------------------------!
         
        do c = -1,-grid % n_bnd_cells,-1
          if(grid % bnd_cond % color(c) .eq. n .and. in_out) then

            ! For velocity and pressure
            vals(0) = u_def;  u % n(c) = vals(Key_Ind('U', keys, nks))
            vals(0) = v_def;  v % n(c) = vals(Key_Ind('V', keys, nks))
            vals(0) = w_def;  w % n(c) = vals(Key_Ind('W', keys, nks))
            vals(0) = p_def;  p % n(c) = vals(Key_Ind('P', keys, nks))

            ! Temperature
            if(heat_transfer == YES) then
              t % bnd_cell_type(c) = bc_type_tag
              if(bc_type_tag .eq. WALLFL) then
                vals(0) = q_def;  t % q(c) = vals(Key_Ind('Q', keys, nks))
              else
                vals(0) = t_def;  t % n(c) = vals(Key_Ind('T', keys, nks))
              end if
            end if

            ! For user scalars
            do us = 1, n_user_scalars
              user_scalar(us) % bnd_cell_type(c) = bc_type_tag
              if(bc_type_tag .eq. WALLFL) then
                write(q_name(3:4),'(i2.2)') us  
                vals(0) = q_def   
                user_scalar(us) % q(c) = vals(Key_Ind(q_name, keys, nks))
              else
                write(c_name(3:4),'(i2.2)') us  
                vals(0) = c_def   
                user_scalar(us) % n(c) = vals(Key_Ind(c_name, keys, nks))
              end if
            end do

            ! For turbulence models
            if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
               turbulence_model == HANJALIC_JAKIRLIC) then
              vals(0) = uu_def;  uu  % n(c) = vals(Key_Ind('UU',  keys, nks))
              vals(0) = vv_def;  vv  % n(c) = vals(Key_Ind('VV',  keys, nks))
              vals(0) = ww_def;  ww  % n(c) = vals(Key_Ind('WW',  keys, nks))
              vals(0) = uv_def;  uv  % n(c) = vals(Key_Ind('UV',  keys, nks))
              vals(0) = uw_def;  uw  % n(c) = vals(Key_Ind('UW',  keys, nks))
              vals(0) = vw_def;  vw  % n(c) = vals(Key_Ind('VW',  keys, nks))
              vals(0) = eps_def; eps % n(c) = vals(Key_Ind('EPS', keys, nks))
              if(turbulence_model == REYNOLDS_STRESS_MODEL) then
                vals(0) = f22_def; f22 % n(c) = vals(Key_Ind('F22', keys, nks))
              end if
            end if

            if(turbulence_model == K_EPS) then
              vals(0)   = kin_def; kin % n(c) = vals(Key_Ind('KIN', keys, nks))
              vals(0)   = eps_def; eps % n(c) = vals(Key_Ind('EPS', keys, nks))
              u_tau(c)  = 0.047
              y_plus(c) = 30.0
            end if

            if(turbulence_model == K_EPS_ZETA_F     .or.  &
               turbulence_model == HYBRID_K_EPS_ZETA_F) then
              vals(0) = kin_def;  kin  % n(c) = vals(Key_Ind('KIN',  keys, nks))
              vals(0) = eps_def;  eps  % n(c) = vals(Key_Ind('EPS',  keys, nks))
              vals(0) = zeta_def; zeta % n(c) = vals(Key_Ind('ZETA', keys, nks))
              vals(0) = f22_def;  f22  % n(c) = vals(Key_Ind('F22',  keys, nks))
            end if

            if(turbulence_model == SPALART_ALLMARAS .or.  &
               turbulence_model == DES_SPALART) then
              vals(0) = vis_def; vis % n(c) = vals(Key_Ind('VIS', keys, nks))
            end if
          end if ! in_out
        end do

      !---------------------------------------------!
      !                                             !
      !   Boundary values are specified in a file   !
      !                                             !
      !---------------------------------------------!
      else
        call Control_Mod_Read_Strings_On('FILE', name_prof, nvs, .false.) 

        open(9, file=name_prof(1))
        if(this_proc < 2) print *, '# Reading the file: ', trim(name_prof(1))
        call Tokenizer_Mod_Read_Line(9)
        read(line % tokens(1),*) n_points  ! number of points

        allocate(prof(n_points, 0:nks))

        !----------------------------------!
        !   Read the entire profile file   !
        !----------------------------------!
        do m = 1, n_points
          call Tokenizer_Mod_Read_Line(9)
          do i = 1, nks
            read(line % tokens(i), *) prof(m, i) 
          end do
        end do
        close(9)

        !------------------------!
        !   A plane is defined   !
        !------------------------!
        if(keys(1) == 'X' .and. keys(2) == 'Y' .or.  &
           keys(1) == 'X' .and. keys(2) == 'Z' .or.  &
           keys(1) == 'Y' .and. keys(2) == 'Z') then       

          ! Set the closest point
          do c = -1,-grid % n_bnd_cells,-1

            ! if in_out true set boundary values, otherwise just type
            if(grid % bnd_cond % color(c) == n .and. in_out) then

              dist_min = HUGE
              do m = 1, n_points

                i = Key_Ind('X', keys, nks); prof(m,0) = 0.0;  x = prof(m,i)
                i = Key_Ind('Y', keys, nks); prof(m,0) = 0.0;  y = prof(m,i)
                i = Key_Ind('Z', keys, nks); prof(m,0) = 0.0;  z = prof(m,i)

                if(keys(1) == 'Y' .and. keys(2) == 'Z') then
                  dist= Distance(y,            z,            0.0,  &
                                 grid % yc(c), grid % zc(c), 0.0)

                else if(keys(1) == 'X' .and. keys(2) == 'Z') then
                  dist= Distance(x,            z,            0.0,  &
                                 grid % xc(c), grid % zc(c), 0.0)

                else if(keys(1) == 'X' .and. keys(2) == 'Y') then
                  dist= Distance(x,            y,            0.0,  &
                                 grid % xc(c), grid % yc(c), 0.0)

                end if

                ! Store closest point in k
                if(dist < dist_min) then
                  dist_min = dist
                  k = m
                end if

              end do

              ! For velocity and pressure
              i=Key_Ind('U',keys,nks); prof(k,0) = u_def; u % n(c) = prof(k,i)
              i=Key_Ind('V',keys,nks); prof(k,0) = v_def; v % n(c) = prof(k,i)
              i=Key_Ind('W',keys,nks); prof(k,0) = w_def; w % n(c) = prof(k,i)
              i=Key_Ind('P',keys,nks); prof(k,0) = p_def; p % n(c) = prof(k,i)

              ! For temperature
              if(heat_transfer == YES) then
                t % bnd_cell_type(c) = bc_type_tag
                if(bc_type_tag .eq. WALLFL) then
                  i=Key_Ind('Q',keys,nks); prof(k,0)=q_def; t % q(c)=prof(k,i)
                else
                  i=Key_Ind('T',keys,nks); prof(k,0)=t_def; t % n(c)=prof(k,i)
                end if
              end if

              ! For user scalars
              do us = 1, n_user_scalars
                user_scalar(us) % bnd_cell_type(c) = bc_type_tag
                if(bc_type_tag .eq. WALLFL) then
                  write(q_name(3:4),'(i2.2)') us  ! set user scalar name
                  i = Key_Ind(q_name, keys, nks)
                  prof(k,0) = q_def
                  user_scalar(us) % q(c) = prof(k,i)
                else
                  write(c_name(3:4),'(i2.2)') us  ! set user scalar name
                  i = Key_Ind(c_name, keys, nks)
                  prof(k,0) = c_def
                  user_scalar(us) % n(c) = prof(k,i)
                end if
              end do

              ! For turbulence models
              if(turbulence_model == K_EPS) then
                i=Key_Ind('KIN',keys,nks); prof(k,0)=kin_def; kin%n(c)=prof(k,i)
                i=Key_Ind('EPS',keys,nks); prof(k,0)=eps_def; eps%n(c)=prof(k,i)
              end if

              if(turbulence_model == K_EPS_ZETA_F) then
                i = Key_Ind('KIN',  keys, nks)
                prof(k,0) = kin_def
                kin % n(c) = prof(k, i)

                i = Key_Ind('EPS',  keys, nks)
                prof(k,0) = eps_def
                eps % n(c) = prof(k, i)

                i = Key_Ind('ZETA', keys, nks)
                prof(k,0) = zeta_def
                zeta % n(c) = prof(k, i)

                i = Key_Ind('F22',  keys, nks)
                prof(k,0) = f22_def
                f22 % n(c) = prof(k, i)
              end if

              if(turbulence_model == SPALART_ALLMARAS .or.  &
                 turbulence_model == DES_SPALART) then
                i=Key_Ind('VIS',keys,nks); prof(k,0)=vis_def; vis%n(c)=prof(k,i)
              end if

              if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
                 turbulence_model == HANJALIC_JAKIRLIC) then
                i=Key_Ind('UU', keys,nks); prof(k,0)=uu_def;  uu %n(c)=prof(k,i)
                i=Key_Ind('VV', keys,nks); prof(k,0)=vv_def;  vv %n(c)=prof(k,i)
                i=Key_Ind('WW', keys,nks); prof(k,0)=ww_def;  ww %n(c)=prof(k,i)
                i=Key_Ind('UV', keys,nks); prof(k,0)=uv_def;  uv %n(c)=prof(k,i)
                i=Key_Ind('UW', keys,nks); prof(k,0)=uw_def;  uw %n(c)=prof(k,i)
                i=Key_Ind('VW', keys,nks); prof(k,0)=vw_def;  vw %n(c)=prof(k,i)
                i=Key_Ind('EPS',keys,nks); prof(k,0)=eps_def; eps%n(c)=prof(k,i)

                if(turbulence_model == REYNOLDS_STRESS_MODEL) then
                  i = Key_Ind('F22', keys, nks); 
                  prof(k,0) = f22_def; f22 % n(c) = prof(k,i)
                end if
              end if        
            end if      !end if(grid % bnd_cond % color(c) == n .and. in_out)
          end do        !end do c = -1,-grid % n_bnd_cells,-1

        !----------------------------!
        !   A plane is not defined   !
        !----------------------------!
        else  ! dir == "XPL" ...
           
          do c = -1,-grid % n_bnd_cells,-1

            ! If in_out is set to true, set boundary values,
            ! otherwise, just the TypeBC remains set.
            if(grid % bnd_cond % color(c) == n .and. in_out) then
          
              do m = 1, n_points-1
                here = .false. 

                i = Key_Ind(keys(1), keys, nks); 
                prof(m,   0) = 0.0;  
                prof(m+1, 0) = 0.0;  
                x  = prof(m,i)
                xp = prof(m+1,i)

                ! Compute the weight factors
                if( keys(1) == 'X' .and.  &
                    grid % xc(c) >= x .and. grid % xc(c) <= xp ) then
                  wi = ( xp - grid % xc(c) ) / (xp - x)
                  here = .true.
                else if( keys(1) == 'Y' .and.  & 
                         grid % yc(c) >= x .and. grid % yc(c) <= xp ) then
                  wi = ( xp - grid % yc(c) ) / (xp - x)
                  here = .true.
                else if( keys(1) == 'Z' .and.  & 
                         grid % zc(c) >= x .and. grid % zc(c) <= xp ) then
                  wi = ( xp - grid % zc(c) ) / (xp - x)
                  here = .true.

                ! Beware; for cylindrical coordinates you have "inversion"
                else if( (keys(1) == 'RX' .and.  &
                     sqrt(grid % yc(c)**2 + grid % zc(c)**2) >= xp .and.      &
                     sqrt(grid % yc(c)**2 + grid % zc(c)**2) <= x) ) then
                  wi = ( xp - sqrt(grid % yc(c)**2 + grid % zc(c)**2) ) / (xp-x)
                  here = .true.
                else if( (keys(1) == 'RY' .and.  &
                     sqrt(grid % xc(c)**2 + grid % zc(c)**2) >= xp .and.      &
                     sqrt(grid % xc(c)**2 + grid % zc(c)**2) <= x) ) then
                  wi = ( xp - sqrt(grid % xc(c)**2 + grid % zc(c)**2) ) / (xp-x)
                  here = .true.
                else if( (keys(1) == 'RZ' .and.  &
                     sqrt(grid % xc(c)**2 + grid % yc(c)**2) >= xp .and.      &
                     sqrt(grid % xc(c)**2 + grid % yc(c)**2) <= x) ) then
                  wi = ( xp - sqrt(grid % xc(c)**2 + grid % yc(c)**2) ) / (xp-x)
                  here = .true.
                end if

                ! Interpolate the profiles     
                if(here) then

                  ! For velocity and pressure
                  i = Key_Ind('U',keys,nks); 
                  u % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  i = Key_Ind('V',keys,nks); 
                  v % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  i = Key_Ind('W',keys,nks); 
                  w % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  i = Key_Ind('P',keys,nks); 
                  p % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                  ! For temperature
                  if(heat_transfer == YES) then
                    t % bnd_cell_type(c) = bc_type_tag
                    if(bc_type_tag .eq. WALLFL) then
                      prof(m,   0) = q_def 
                      prof(m+1, 0) = q_def
                      i = Key_Ind('Q',keys,nks); 
                      t % q(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                    else
                      prof(m,   0) = t_def 
                      prof(m+1, 0) = t_def
                      i = Key_Ind('T',keys,nks); 
                      t % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                    end if
                  end if

                  ! For user scalars
                  do us = 1, n_user_scalars
                    user_scalar(us) % bnd_cell_type(c) = bc_type_tag
                    if(bc_type_tag .eq. WALLFL) then
                      prof(m,   0) = q_def 
                      prof(m+1, 0) = q_def
                      write(q_name(3:4),'(i2.2)') us  ! set variable name
                      i = Key_Ind(q_name, keys, nks)
                      user_scalar(us) % q(c) = wi*prof(m,i)+(1.-wi)*prof(m+1,i)
                    else
                      prof(m,   0) = c_def 
                      prof(m+1, 0) = c_def
                      write(c_name(3:4),'(i2.2)') us  ! set variable name
                      i = Key_Ind(c_name, keys, nks)
                      user_scalar(us) % n(c) = wi*prof(m,i)+(1.-wi)*prof(m+1,i)
                    end if
                  end do

                  ! For turbulence models
                  if(turbulence_model == K_EPS) then

                    prof(m,   0) = kin_def 
                    prof(m+1, 0) = kin_def
                    i = Key_Ind('KIN',keys,nks); 
                    kin % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    prof(m,   0) = eps_def 
                    prof(m+1, 0) = eps_def
                    i = Key_Ind('EPS',keys,nks); 
                    eps % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  end if

                  if(turbulence_model == K_EPS_ZETA_F  .or.  &
                     turbulence_model == HYBRID_K_EPS_ZETA_F) then

                    prof(m,   0) = kin_def 
                    prof(m+1, 0) = kin_def
                    i = Key_Ind('KIN',keys,nks); 
                    kin % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    prof(m,   0) = eps_def 
                    prof(m+1, 0) = eps_def
                    i = Key_Ind('EPS',keys,nks); 
                    eps % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    prof(m,   0) = zeta_def 
                    prof(m+1, 0) = zeta_def
                    i = Key_Ind('ZETA',keys,nks); 
                    zeta % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    prof(m,   0) = f22_def 
                    prof(m+1, 0) = f22_def
                    i = Key_Ind('F22',keys,nks); 
                    f22 % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  end if

                  if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
                     turbulence_model == HANJALIC_JAKIRLIC) then

                    prof(m,   0) = uu_def 
                    prof(m+1, 0) = uu_def
                    i = Key_Ind('UU',keys,nks); 
                    uu % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    prof(m,   0) = vv_def 
                    prof(m+1, 0) = vv_def
                    i = Key_Ind('VV',keys,nks); 
                    vv % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    prof(m,   0) = ww_def 
                    prof(m+1, 0) = ww_def
                    i = Key_Ind('WW',keys,nks); 
                    ww % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    prof(m,   0) = uv_def 
                    prof(m+1, 0) = uv_def
                    i = Key_Ind('UV',keys,nks); 
                    uv % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    prof(m,   0) = uw_def 
                    prof(m+1, 0) = uw_def
                    i = Key_Ind('UW',keys,nks); 
                    uw % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    prof(m,   0) = vw_def 
                    prof(m+1, 0) = vw_def
                    i = Key_Ind('VW',keys,nks); 
                    vw % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    prof(m,   0) = eps_def 
                    prof(m+1, 0) = eps_def
                    i = Key_Ind('EPS',keys,nks); 
                    eps % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                    if(turbulence_model == REYNOLDS_STRESS_MODEL) then
                      i = Key_Ind('F22', keys, nks); 
                      prof(k,0) = f22_def; f22 % n(c) = prof(k,i)
                    end if
                  end if        

                  if(turbulence_model == SPALART_ALLMARAS .or.  &
                     turbulence_model == DES_SPALART) then
                    prof(m,   0) = vis_def 
                    prof(m+1, 0) = vis_def
                    i = Key_Ind('VIS',keys,nks); 
                    vis % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  end if

                end if  ! (here)
              end do  ! m = 1, n_points-1
            end if  ! if(in_out)
          end do  ! c = -1,-grid % n_bnd_cells,-1
        end if  ! plane is defined?
        close(9)
      end if  ! boundary defined in a file
    end do

  end do

  end subroutine
