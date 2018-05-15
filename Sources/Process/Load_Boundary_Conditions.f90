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
  integer           :: c, m, k, i, n, n_points, nks, nvs, found
  character(len=80) :: name_prof(128)
  real              :: wi, dist_min, x, y, z, xp, dist
  real, allocatable :: prof(:,:)
  logical           :: here
  character(len=80) :: bc_type
  character(len=80) :: keys(128)
  real              :: vals(0:128)     ! Note that they start from zero!

  ! Default values for boundary conditions
  real, parameter   :: u_def   = 0.0,  v_def   = 0.0,  w_def    = 0.0
  real, parameter   :: p_def   = 0.0,  t_def   = 0.0,  q_def    = 0.0
  real, parameter   :: kin_def = 0.0,  eps_def = 0.0,  f22_def  = 0.0
  real, parameter   :: vis_def = 0.0,  zeta_def = 0.0
  real, parameter   :: uu_def  = 0.0,  vv_def  = 0.0,  ww_def   = 0.0
  real, parameter   :: uv_def  = 0.0,  uw_def  = 0.0,  vw_def   = 0.0
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

  !------------------------------------------------!
  !                                                !
  !   Read boundary conditions from control file   !
  !                                                !
  !------------------------------------------------!
  do n = 1, grid % n_bnd_cond

    ! Find the line where boundary condition defintion is defined
    call Control_Mod_Position_At_Two_Keys('BOUNDARY_CONDITION',       &
                                          grid % bnd_cond % name(n),  &
                                          found,                      &
                                          .true.)

    !---------------------------------------!
    !   Boundary condition has been found   !
    !---------------------------------------!
    if(found == YES) then 

      call Control_Mod_Read_Char_Item_On('TYPE', 'WALL', bc_type, .true.)
      call To_Upper_Case(bc_type)

      call Control_Mod_Read_Strings_On('VARIABLES', keys, nks, .true.)

      ! Try to see if there is file specified
      call Control_Mod_Read_Strings_On('FILE', name_prof, nvs, .true.)

      ! File was specified
      if(nvs == 1) then
        if(this_proc < 2) &
          print *, '# Values specified in the file: ', trim(name_prof(nvs))

      ! If file was not specified, so must have been the values
      else if(nvs == 0) then 
        if(this_proc < 2) &
          print *, '# Values are not specified in a file'

        ! Go back to key and read all sections before values
        call Control_Mod_Position_At_Two_Keys('BOUNDARY_CONDITION',       &
                                              grid % bnd_cond % name(n),  &
                                              found,                      &
                                             .true.)
        call Control_Mod_Read_Char_Item_On('TYPE', 'WALL', bc_type, .true.)
        call To_Upper_Case(bc_type)
        call Control_Mod_Read_Strings_On('VARIABLES', keys,    nks, .false.)
        call Control_Mod_Read_Real_Array_On('VALUES', vals(1), nvs, .false.)

        ! Check validity of the input
        ! if(nks .eq. 0 .or. nvs .eq. 0) then
        !   print '(3a)', '# Critical, for boundary condition ',        &
        !                 trim(grid % bnd_cond % name(n)),              &
        !                 ' no values or variables have been provided' 
        !   stop
        ! end if
        if(nks .ne. nvs) then
          print '(3a)', '# Critical, number of values for boundary condition ',  &
                        trim(grid % bnd_cond % name(n)),                         &
                        ' is not the same as number of provided variable names' 
          stop
        end if
 
      end if

      ! Input is valid, turn keys to upper case
      do i = 1, nks
        call To_Upper_Case(keys(i))
      end do

    !-------------------------------------------!
    !   Boundary condition has not been found   !
    !-------------------------------------------!
    else   
      print '(3a)', '# Critical, boundary condition ',  &
                      trim(grid % bnd_cond % name(n)),                         &
                      ' has not been found.  Exiting!' 
      stop
    end if
    
    if( bc_type == 'INFLOW') then 
      grid % bnd_cond % type(n) = INFLOW
      PER_BC = NO
    else if( bc_type == 'WALL') then 
      grid % bnd_cond % type(n) = WALL
    else if( bc_type == 'OUTFLOW') then 
      grid % bnd_cond % type(n) = OUTFLOW
    else if( bc_type == 'SYMMETRY') then 
      grid % bnd_cond % type(n)=SYMMETRY
    else if( bc_type == 'HEAT_FLUX') then 
      grid % bnd_cond % type(n) = WALLFL
    else if( bc_type == 'CONVECTIVE') then 
      grid % bnd_cond % type(n) = CONVECT
    else if( bc_type == 'PRESSURE') then 
      grid % bnd_cond % type(n) = PRESSURE
    else
      if(this_proc < 2)  &
        print *, '# Load_Boundary_Conditions: '//          &
                   '# Unknown boundary condition type: ',  &
                   bc_type
      stop  
    end if

    !------------------------------------------------------!
    !   Boundary condition is given by a single constant   !
    !------------------------------------------------------!
    if(nks == nvs) then 

      do c = -1,-grid % n_bnd_cells,-1
        if(grid % bnd_cond % color(c) == n) then

          ! If in_out is set to true, set boundary values,
          ! otherwise, just the ypTeBC remains set.
          if(in_out) then

            vals(0) = u_def;  u % n(c) = vals(Key_Ind('U', keys, nks))
            vals(0) = v_def;  v % n(c) = vals(Key_Ind('V', keys, nks))
            vals(0) = w_def;  w % n(c) = vals(Key_Ind('W', keys, nks))
            vals(0) = p_def;  p % n(c) = vals(Key_Ind('P', keys, nks))

            if(heat_transfer == YES) then
              if(grid % bnd_cond % type(n) .eq. WALLFL) then
                vals(0) = q_def;  t % q(c) = vals(Key_Ind('Q', keys, nks))
              else
                vals(0) = t_def;  t % n(c) = vals(Key_Ind('T', keys, nks))
              endif
            end if  ! for heat_transfer == YES

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

        end if 
      end do

    !------------------------------------------------!
    !   Boundary condition is prescribed in a file   !
    !------------------------------------------------!
    else

      open(9, file=name_prof(1))
      if(this_proc < 2) print *, '# Reading the file: ', trim(name_prof(1))
      call Tokenizer_Mod_Read_Line(9)
      read(line % tokens(1),*) n_points  ! number of points

      allocate(prof(n_points, 0:nks))

      ! Read the entire profile file
      do m = 1, n_points
        call Tokenizer_Mod_Read_Line(9)
        do i = 1, nks
          read(line % tokens(i), *) prof(m, i) 
        end do
      end do
      close(9)

      ! A plane is defined
      if(keys(1) == 'X' .and. keys(2) == 'Y' .or.  &
         keys(1) == 'X' .and. keys(2) == 'Z' .or.  &
         keys(1) == 'Y' .and. keys(2) == 'Z') then       

        ! Set the closest point
        do c = -1,-grid % n_bnd_cells,-1
          if(grid % bnd_cond % color(c) == n) then

            if(in_out) then  ! if true set boundary values, otherwise just type
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

              i=Key_Ind('U',keys,nks); prof(k,0) = u_def; u % n(c) = prof(k,i)
              i=Key_Ind('V',keys,nks); prof(k,0) = v_def; v % n(c) = prof(k,i)
              i=Key_Ind('W',keys,nks); prof(k,0) = w_def; w % n(c) = prof(k,i)

              if(heat_transfer == YES) then
                i=Key_Ind('T',keys,nks); prof(k,0)=t_def; t%n(c)=prof(k,i)
              end if

              if(turbulence_model == K_EPS) then
                i=Key_Ind('KIN',keys,nks); prof(k,0)=kin_def; kin%n(c)=prof(k,i)
                i=Key_Ind('EPS',keys,nks); prof(k,0)=eps_def; eps%n(c)=prof(k,i)
              end if

              if(turbulence_model == K_EPS_ZETA_F) then
                i=Key_Ind('KIN',  keys, nks); prof(k,0) = kin_def;  kin  % n(c) = prof(k,i)
                i=Key_Ind('EPS',  keys, nks); prof(k,0) = eps_def;  eps  % n(c) = prof(k,i)
                i=Key_Ind('ZETA', keys, nks); prof(k,0) = zeta_def; zeta % n(c) = prof(k,i)
                i=Key_Ind('F22',  keys, nks); prof(k,0) = f22_def;  f22  % n(c) = prof(k,i)
              end if

              if(turbulence_model == DES_SPALART) then
                i=Key_Ind('VIS',keys,nks); prof(k,0)=vis_def; vis%n(c)=prof(k,i)
              end if

              if(turbulence_model == REYNOLDS_STRESS_MODEL) then
                i=Key_Ind('UU', keys,nks); prof(k,0)=uu_def;  uu %n(c)=prof(k,i)
                i=Key_Ind('VV', keys,nks); prof(k,0)=vv_def;  vv %n(c)=prof(k,i)
                i=Key_Ind('WW', keys,nks); prof(k,0)=ww_def;  ww %n(c)=prof(k,i)
                i=Key_Ind('UV', keys,nks); prof(k,0)=uv_def;  uv %n(c)=prof(k,i)
                i=Key_Ind('UW', keys,nks); prof(k,0)=uw_def;  uw %n(c)=prof(k,i)
                i=Key_Ind('VW', keys,nks); prof(k,0)=vw_def;  vw %n(c)=prof(k,i)
                i=Key_Ind('F22',keys,nks); prof(k,0)=f22_def; f22%n(c)=prof(k,i)
                i=Key_Ind('EPS',keys,nks); prof(k,0)=eps_def; eps%n(c)=prof(k,i)
              end if        
            end if    !end if(in_out)
          end if      !end if(grid % bnd_cond % color(c) == n)
        end do        !end do c = -1,-grid % n_bnd_cells,-1

      ! A plane is not defined
      else  ! dir == "XPL" ...
           
        do c = -1,-grid % n_bnd_cells,-1
          if(grid % bnd_cond % color(c) == n) then
          
            ! If in_out is set to true, set boundary values,
            ! otherwise, just the TypeBC remains set.
            if(in_out) then
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
                  i = Key_Ind('U',keys,nks); 
                  u % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  i = Key_Ind('V',keys,nks); 
                  v % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  i = Key_Ind('W',keys,nks); 
                  w % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)

                  if(heat_transfer == YES) then
                    prof(m,   0) = t_def 
                    prof(m+1, 0) = t_def
                    i = Key_Ind('T',keys,nks); 
                    t % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  end if

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

                  if(turbulence_model == SPALART_ALLMARAS) then
                    prof(m,   0) = vis_def 
                    prof(m+1, 0) = vis_def
                    i = Key_Ind('VIS',keys,nks); 
                    vis % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  end if

                  if(turbulence_model == DES_SPALART) then
                    prof(m,   0) = vis_def 
                    prof(m+1, 0) = vis_def
                    i = Key_Ind('VIS',keys,nks); 
                    vis % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                  end if

                end if
              end do
            end if  ! if(in_out)
          end if 
        end do
      end if
      close(9)
    end if
  end do 

  end subroutine
