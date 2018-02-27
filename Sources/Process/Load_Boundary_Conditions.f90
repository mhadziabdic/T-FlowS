!==============================================================================!
  subroutine Load_Boundary_Conditions(grid, in_out)
!------------------------------------------------------------------------------!
!   Reads: .bnd file                                                           !
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
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
  integer           :: c, n, n_points, n_initial_cond, nks, n_vals
  integer           :: m, c1, bc, mt, i
  character(len=80) :: name_bou, name_prof(128), bc_name, mt_name
  integer           :: type_bnd_cond(128)
  real              :: wi
  real              :: Mres
  real, allocatable :: prof(:,:)
  logical           :: here
  integer           :: found 
  character         :: name_ini(128)*80
  character(len=80) :: bc_type
  character(len=80) :: keys(128)
  real              :: vals(0:128)     ! Note that they start from zero!
  real              :: x, y, z, xp, dist

  ! Default values for boundary conditions (how about initial?)
  real, parameter   :: u_def   = 0.0,  v_def   = 0.0,  w_def   = 0.0
  real, parameter   :: p_def   = 0.0,  t_def   = 0.0,  q_def   = 0.0
  real, parameter   :: kin_def = 0.0,  eps_def = 0.0,  f22_def = 0.0
  real, parameter   :: vis_def = 0.0,  v2_def  = 0.0
  real, parameter   :: uu_def  = 0.0,  vv_def  = 0.0,  ww_def  = 0.0
  real, parameter   :: uv_def  = 0.0,  uw_def  = 0.0,  vw_def  = 0.0
!==============================================================================!

  !--------------------------------------------!
  !   Read the file with boundary conditions   !
  !--------------------------------------------!
  name_bou = problem_name
  name_bou(len_trim(problem_name)+1:len_trim(problem_name)+4) = '.bnd'
  open(9, file=name_bou)
  if(this_proc < 2) print *, '# Reading the file: ', trim(name_bou)

  !-------------------------!
  !   Phisical properties   !
  !-------------------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) grid % n_materials
  do mt = 1,grid % n_materials

    call Tokenizer_Mod_Read_Line(9)
    call To_Upper_Case(  line % tokens(1)  )
    call To_Upper_Case(  line % tokens(2)  )
    read(line % tokens(1),*) mt_name

    ! Find material index
    do i=1, grid % n_materials 
      if(mt_name == grid % materials(i) % name) n=i      
    end do

    if( line % tokens(2) == 'FLUID') then 
      StateMat(n)=FLUID
    else if( line % tokens(2) == 'SOLID') then 
      StateMat(n)=SOLID
    else 
      if(this_proc < 2)  &
        print *, '# Load_Boundary_Conditions: Unknown material state'
      stop  
    end if
    read(line % tokens(3),*) viscosity
    read(line % tokens(4),*) density
    if(heat_transfer == YES) read(line % tokens(5),*) conductivity
    if(heat_transfer == YES) read(line % tokens(6),*) capacity 
  end do
  
  !-----------------------------------------------------!
  !   Boundary conditions 1 - read them from the file   !
  !-----------------------------------------------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) grid % n_bnd_cond
  print *, '# Found ', grid % n_bnd_cond, ' boundary conditions'

  do bc = 1, grid % n_bnd_cond  ! number of boundary conditions

    call Tokenizer_Mod_Read_Line(9)
    call To_Upper_Case(  line % tokens(1)  )
    call To_Upper_Case(  line % tokens(2)  )
    call To_Upper_Case(  line % tokens(3)  )
    read(line % tokens(1),*) bc_name

    ! Find b.c. index
    n = -1
    do i=1, grid % n_bnd_cond 
      if(bc_name .eq. grid % bnd_cond % name(i)) n = i      
    end do
    if( n == -1 ) then
      print *, '# Critical, failed to find boundary condition ', bc_name
      print *, '# Exiting!'
      stop
    end if 

    if( line % tokens(2) == 'INFLOW') then 
      type_bnd_cond(n)=INFLOW
      PER_BC = NO
    else if( line % tokens(2) == 'WALL') then 
      type_bnd_cond(n)=WALL
    else if( line % tokens(2) == 'OUTFLOW') then 
      type_bnd_cond(n)=OUTFLOW
    else if( line % tokens(2) == 'SYMMETRY') then 
      type_bnd_cond(n)=SYMMETRY
    else if( line % tokens(2) == 'WALLFLUX') then 
      type_bnd_cond(n)=WALLFL
    else if( line % tokens(2) == 'CONVECTIVE') then 
      type_bnd_cond(n)=CONVECT
    else if( line % tokens(2) == 'PRESSURE') then 
      type_bnd_cond(n)=PRESSURE
    else
      if(this_proc < 2)  &
        print *, '# Load_Boundary_Conditions: '//        &
                   '# Unknown boundary condition type: ',  &
                   line % tokens(2)
      stop  
    end if
    if( line % tokens(3)  ==  'FILE') then
      read(line % tokens(4),'(A80)') name_prof(n)
    else
      read(line % tokens(3),*) u % bound(n)
      read(line % tokens(4),*) v % bound(n)
      read(line % tokens(5),*) w % bound(n)
      if(type_bnd_cond(n)==PRESSURE) then
        read(line % tokens(6),*) P % bound(n)
        if(heat_transfer == YES) then 
          read(line % tokens(7),*) T % bound(n)
          if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
             turbulence_model == HANJALIC_JAKIRLIC) then
            read(line % tokens(8),*)   uu % bound(n)
            read(line % tokens(9),*)   vv % bound(n)
            read(line % tokens(10),*) ww % bound(n)
            read(line % tokens(11),*) uv % bound(n)
            read(line % tokens(12),*) uw % bound(n)
            read(line % tokens(13),*) vw % bound(n)
            read(line % tokens(14),*) eps% bound(n)
            if(turbulence_model == REYNOLDS_STRESS_MODEL)  &
              read(line % tokens(15),*) f22 % bound(n)
          end if
          if(turbulence_model == K_EPS) then
            read(line % tokens(8),*) kin % bound(n)
            read(line % tokens(9),*) eps % bound(n)
          end if
          if(turbulence_model == K_EPS_V2 .or.  &
             turbulence_model == K_EPS_ZETA_F     .or.  &
             turbulence_model == HYBRID_K_EPS_ZETA_F) then
            read(line % tokens(8),*) kin % bound(n)
            read(line % tokens(9),*) eps % bound(n)
            read(line % tokens(10),*) v2  % bound(n)
            read(line % tokens(11),*) f22 % bound(n)
          end if
          if(turbulence_model == SPALART_ALLMARAS) then
            read(line % tokens(8),*) VIS % bound(n)
          end if
          if(turbulence_model == DES_SPALART) then
            read(line % tokens(8),*) VIS % bound(n)
          end if
        else  ! heat_transfer .ne. YES
          if(turbulence_model == REYNOLDS_STRESS_MODEL .or.   &
             turbulence_model == HANJALIC_JAKIRLIC) then
            read(line % tokens(7),*)   uu % bound(n)
            read(line % tokens(8),*)   vv % bound(n)
            read(line % tokens(9),*) ww % bound(n)
            read(line % tokens(10),*) uv % bound(n)
            read(line % tokens(11),*) uw % bound(n)
            read(line % tokens(12),*) vw % bound(n)
            read(line % tokens(13),*) eps% bound(n)
            if(turbulence_model == REYNOLDS_STRESS_MODEL)  &
              read(line % tokens(14),*) f22 % bound(n)
          end if
          if(turbulence_model == K_EPS) then
            read(line % tokens(7),*) kin % bound(n)
            read(line % tokens(8),*) eps % bound(n)
          end if
          if(turbulence_model == K_EPS_V2 .or.  &
             turbulence_model == K_EPS_ZETA_F     .or.  &
             turbulence_model == HYBRID_K_EPS_ZETA_F) then
            read(line % tokens(7),*) kin % bound(n)
            read(line % tokens(8),*) eps % bound(n)
            read(line % tokens(9),*) v2   % bound(n)
            read(line % tokens(10),*) f22 % bound(n)
          end if
          if(turbulence_model == SPALART_ALLMARAS) then
            read(line % tokens(7),*) VIS % bound(n)
          end if
          if(turbulence_model == DES_SPALART) then
            read(line % tokens(7),*) VIS % bound(n)
          end if
        end if  ! heat_transfer == YES
        name_prof=''
      else   ! type_bnd_cond .ne. PRESSURE
        if(heat_transfer == YES) then 
          read(line % tokens(6),*) T % bound(n)
          if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
             turbulence_model == HANJALIC_JAKIRLIC) then
            read(line % tokens(7),*) uu % bound(n)
            read(line % tokens(8),*) vv % bound(n)
            read(line % tokens(9),*) ww % bound(n)
            read(line % tokens(10),*) uv % bound(n)
            read(line % tokens(11),*) uw % bound(n)
            read(line % tokens(12),*) vw % bound(n)
            read(line % tokens(13),*) eps% bound(n)
            if(turbulence_model == REYNOLDS_STRESS_MODEL)  &
              read(line % tokens(14),*) f22 % bound(n)
          end if
          if(turbulence_model == K_EPS) then
            read(line % tokens(7),*) kin % bound(n)
            read(line % tokens(8),*) eps % bound(n)
          end if
          if(turbulence_model == K_EPS_V2 .or.  &
             turbulence_model == K_EPS_ZETA_F     .or.  &
             turbulence_model == HYBRID_K_EPS_ZETA_F) then
            read(line % tokens(7),*) kin % bound(n)
            read(line % tokens(8),*) eps % bound(n)
            read(line % tokens(9),*) v2  % bound(n)
            read(line % tokens(10),*) f22 % bound(n)
          end if
          if(turbulence_model == SPALART_ALLMARAS) then
            read(line % tokens(7),*) VIS % bound(n)
          end if
          if(turbulence_model == DES_SPALART) then
            read(line % tokens(7),*) VIS % bound(n)
          end if
        else  ! HOT .ne. YES
          if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
             turbulence_model == HANJALIC_JAKIRLIC) then
            read(line % tokens(6),*) uu % bound(n)
            read(line % tokens(7),*) vv % bound(n)
            read(line % tokens(8),*) ww % bound(n)
            read(line % tokens(9),*) uv % bound(n)
            read(line % tokens(10),*) uw % bound(n)
            read(line % tokens(11),*) vw % bound(n)
            read(line % tokens(12),*) eps% bound(n)
            if(turbulence_model == REYNOLDS_STRESS_MODEL)  &
              read(line % tokens(13),*) f22 % bound(n)
          end if
          if(turbulence_model == K_EPS) then
            read(line % tokens(6),*) kin % bound(n)
            read(line % tokens(7),*) eps % bound(n)
          end if
          if(turbulence_model == K_EPS_V2 .or.  &
             turbulence_model == K_EPS_ZETA_F     .or.  &
             turbulence_model == HYBRID_K_EPS_ZETA_F) then
            read(line % tokens(6),*) kin % bound(n)
            read(line % tokens(7),*) eps % bound(n)
            read(line % tokens(8),*) v2   % bound(n)
            read(line % tokens(9),*) f22 % bound(n)
          end if
          if(turbulence_model == SPALART_ALLMARAS) then
            read(line % tokens(6),*) VIS % bound(n)
          end if
          if(turbulence_model == DES_SPALART) then
            read(line % tokens(6),*) VIS % bound(n)
          end if
        end if  ! heat_transfer == YES
        name_prof=''
      end if  ! type_bnd_cond == PRESSURE
    end if    
  end do      

  !------------------------!
  !   Initial conditions   !
  !------------------------!
  call Tokenizer_Mod_Read_Line(9)
  read(line % tokens(1), *) n_initial_cond
  print *, '# Number of initial conditions: ', n_initial_cond
  if(n_initial_cond > grid % n_materials) then
    if(this_proc < 2)   &
      print *, 'Warning: there are more initial conditions then materials'
  end if

  do n=1,n_initial_cond
    call Tokenizer_Mod_Read_Line(9)
    call To_Upper_Case(line % tokens(2))

    ! Initial conditions given in GMV file
    if(line % tokens(2) == 'FILE') then
      read(line % tokens(3),'(A80)') name_ini(n)
      print *, '# Load_Boundary_Conditions: material ', n,  &
                 '; init. cond. given by file: ', name_ini(n)
    else
      name_ini(n) = ''

      ! Initial conditions given by constant
      read(line % tokens(2),*) u % init(n)
      read(line % tokens(3),*) v % init(n)
      read(line % tokens(4),*) w % init(n)
 
      if(heat_transfer == YES) then
        read(line % tokens(5),*) T % init(n)
        if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
           turbulence_model == HANJALIC_JAKIRLIC) then
          read(line % tokens(6),*) uu % init(n)
          read(line % tokens(7),*) vv % init(n)
          read(line % tokens(8),*) ww % init(n)
          read(line % tokens(9),*) uv % init(n)
          read(line % tokens(10),*) uw % init(n)
          read(line % tokens(11),*) vw % init(n)
          read(line % tokens(12),*) eps% init(n)
          if(turbulence_model == REYNOLDS_STRESS_MODEL)  &
            read(line % tokens(13),*) f22 % init(n)
        end if
        if(turbulence_model == K_EPS) then
          read(line % tokens(6),*) kin % init(n)
          read(line % tokens(7),*) eps % init(n)
        end if
        if(turbulence_model == K_EPS_V2 .or.  &
           turbulence_model == K_EPS_ZETA_F     .or.  &
           turbulence_model == HYBRID_K_EPS_ZETA_F) then
          read(line % tokens(6),*) kin % init(n)
          read(line % tokens(7),*) eps % init(n)
          read(line % tokens(8),*) v2   % init(n)
          read(line % tokens(9),*) f22 % init(n)
        end if
        if(turbulence_model == SPALART_ALLMARAS) then
          read(line % tokens(6),*) VIS % init(n)
        end if
        if(turbulence_model == DES_SPALART) then
          read(line % tokens(6),*) VIS % init(n)
        end if
      else ! HOT /= YES
        if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
           turbulence_model == HANJALIC_JAKIRLIC) then
          read(line % tokens(5),*) uu % init(n)
          read(line % tokens(6),*) vv % init(n)
          read(line % tokens(7),*) ww % init(n)
          read(line % tokens(8),*) uv % init(n)
          read(line % tokens(9),*) uw % init(n)
          read(line % tokens(10),*) vw % init(n)
          read(line % tokens(11),*) eps% init(n)
          if(turbulence_model == REYNOLDS_STRESS_MODEL)  &
            read(line % tokens(12),*) f22 % init(n)
        end if
        if(turbulence_model == K_EPS) then
          read(line % tokens(5),*) kin % init(n)
          read(line % tokens(6),*) eps % init(n)
        end if
        if(turbulence_model == K_EPS_V2 .or.  &
           turbulence_model == K_EPS_ZETA_F     .or.  &
           turbulence_model == HYBRID_K_EPS_ZETA_F) then
          read(line % tokens(5),*) kin % init(n)
          read(line % tokens(6),*) eps % init(n)
          read(line % tokens(7),*) v2   % init(n)
          read(line % tokens(8),*) f22 % init(n)
        end if
        if(turbulence_model == SPALART_ALLMARAS) then
          read(line % tokens(5),*) VIS % init(n)
        end if
        if(turbulence_model == DES_SPALART) then
          read(line % tokens(5),*) VIS % init(n)
        end if
      end if
    end if
  end do  

  close(9)

  !--------------------------------------------------------------!
  !   Store types of boundary conditions in the grid structure   !
  !--------------------------------------------------------------!
  do n = 1, grid % n_bnd_cond
    grid % bnd_cond % type(n) = type_bnd_cond(n)
print *, type_bnd_cond(n)
  end do

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

do n = 1, grid % n_bnd_cond
print *, 'A', n, type_bnd_cond(n)
end do


  !----------------------------------------------------------------------!
  !                                                                      !
  !   Boundary conditions 2 - distribute them over computational cells   !
  !                                                                      !
  !----------------------------------------------------------------------!
  do n = 1, grid % n_bnd_cond

print *, '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print *, '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print *, 'B', n, type_bnd_cond(n)

    ! Found the line where boundary condition defintion is defined
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

print *, '# Found BOUNDARY_CONDITION', trim(grid % bnd_cond % name(n))
print *, '# and it is of TYPE       ', trim(bc_type)

      call Control_Mod_Read_Strings_On('VARIABLES', keys, nks, .true.)

      ! Try to see if there is file specified
      call Control_Mod_Read_Strings_On('FILE', name_prof, n_vals, .true.)

      ! File was specified
      if(n_vals == 1) then
        print *, '# Values specified in the file: ', trim(name_prof(n_vals))

      ! If file was not specified, so must have been the values
      else if(n_vals == 0) then 
        print *, '# Values are not specified in a file'

        ! Go back to key and read all sections before values
        call Control_Mod_Position_At_Two_Keys('BOUNDARY_CONDITION',       &
                                              grid % bnd_cond % name(n),  &
                                              found,                      &
                                             .true.)
        call Control_Mod_Read_Char_Item_On('TYPE', 'WALL', bc_type, .true.)
        call To_Upper_Case(bc_type)
        call Control_Mod_Read_Strings_On('VARIABLES', keys,      nks, .true.)
        call Control_Mod_Read_Real_Array_On('VALUES', vals(1), n_vals, .true.)

        ! Check validity of the input
        if(nks .eq. 0 .or. n_vals .eq. 0) then
          print '(3a)', '# Critical, for boundary condition ',        &
                        trim(grid % bnd_cond % name(n)),              &
                        ' no values or variables have been provided' 
          stop
        end if

        if(nks .ne. n_vals) then
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
      type_bnd_cond(n) = INFLOW
      PER_BC = NO
    else if( bc_type == 'WALL') then 
      type_bnd_cond(n) = WALL
    else if( bc_type == 'OUTFLOW') then 
      type_bnd_cond(n) = OUTFLOW
    else if( bc_type == 'SYMMETRY') then 
      type_bnd_cond(n)=SYMMETRY
    else if( bc_type == 'HEAT_FLUX') then 
      type_bnd_cond(n) = WALLFL
    else if( bc_type == 'CONVECTIVE') then 
      type_bnd_cond(n) = CONVECT
    else if( bc_type == 'PRESSURE') then 
      type_bnd_cond(n) = PRESSURE
    else
      if(this_proc < 2)  &
        print *, '# Load_Boundary_Conditions: '//        &
                   '# Unknown boundary condition type: ',  &
                   bc_type
      stop  
    end if
    grid % bnd_cond % type(n) = type_bnd_cond(n)

    ! Boundary condition is given by a single constant
    if(nks == n_vals) then 

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
              if(type_bnd_cond(n) .eq. WALLFL) then
                vals(0) = q_def;  t % q(c) = vals(Key_Ind('T', keys, nks))
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
              vals(0) = kin_def; kin % n(c) = vals(Key_Ind('KIN', keys, nks))
              vals(0) = eps_def; eps % n(c) = vals(Key_Ind('EPS', keys, nks))
              Uf(c)   = 0.047
              Ynd(c)  = 30.0
            end if

            if(turbulence_model == K_EPS_V2 .or.  &
               turbulence_model == K_EPS_ZETA_F     .or.  &
               turbulence_model == HYBRID_K_EPS_ZETA_F) then
              vals(0) = kin_def; kin % n(c) = vals(Key_Ind('KIN', keys, nks))
              vals(0) = eps_def; eps % n(c) = vals(Key_Ind('EPS', keys, nks))
              vals(0) = f22_def; f22 % n(c) = vals(Key_Ind('F22', keys, nks))
              vals(0) = v2_def;  v2  % n(c) = vals(Key_Ind('V2',  keys, nks))
            end if

            if(turbulence_model == SPALART_ALLMARAS) then
              vals(0) = vis_def; vis % n(c) = vals(Key_Ind('VIS', keys, nks))
            end if
            if(turbulence_model == DES_SPALART) then
              vals(0) = vis_def; vis % n(c) = vals(Key_Ind('VIS', keys, nks))
            end if
          end if ! in_out

        end if 
      end do

    ! Boundary condition is prescribed in a file 
    else

      open(9, file=name_prof(1))
      if(this_proc < 2) print *, '# Reading the file: ', trim(name_prof(1))
      call Tokenizer_Mod_Read_Line(9)
      read(line % tokens(1),*) n_points                  ! number of points

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

            if(in_out) then    !if .true. set boundary values, otherwise, just set TypeBC
              Mres = HUGE
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

                if(dist < Mres) then
                  Mres = dist
                  c1 = m
                end if

              end do

              i=Key_Ind('U',keys,nks); prof(m,0) = u_def; u % n(c) = prof(m,i)
              i=Key_Ind('V',keys,nks); prof(m,0) = v_def; v % n(c) = prof(m,i)
              i=Key_Ind('W',keys,nks); prof(m,0) = w_def; w % n(c) = prof(m,i)

              if(heat_transfer == YES) then
                i=Key_Ind('T',keys,nks); prof(m,0)=t_def; t%n(c)=prof(m,i)
              end if

              if(turbulence_model == K_EPS) then
                i=Key_Ind('KIN',keys,nks); prof(m,0)=kin_def; kin%n(c)=prof(m,i)
                i=Key_Ind('EPS',keys,nks); prof(m,0)=eps_def; eps%n(c)=prof(m,i)
              end if

              if(turbulence_model == K_EPS_V2 .or.  &
                 turbulence_model == K_EPS_ZETA_F) then
                i=Key_Ind('KIN',keys,nks); prof(m,0)=kin_def; kin%n(c)=prof(m,i)
                i=Key_Ind('EPS',keys,nks); prof(m,0)=eps_def; eps%n(c)=prof(m,i)
                i=Key_Ind('V2', keys,nks); prof(m,0)=v2_def;  v2 %n(c)=prof(m,i)
                i=Key_Ind('F22',keys,nks); prof(m,0)=f22_def; f22%n(c)=prof(m,i)
              end if

              if(turbulence_model == DES_SPALART) then
                i=Key_Ind('VIS',keys,nks); prof(m,0)=vis_def; vis%n(c)=prof(m,i)
              end if

              if(turbulence_model == REYNOLDS_STRESS_MODEL) then
                i=Key_Ind('UU', keys,nks); prof(m,0)=uu_def;  uu %n(c)=prof(m,i)
                i=Key_Ind('VV', keys,nks); prof(m,0)=vv_def;  vv %n(c)=prof(m,i)
                i=Key_Ind('WW', keys,nks); prof(m,0)=ww_def;  ww %n(c)=prof(m,i)
                i=Key_Ind('UV', keys,nks); prof(m,0)=uv_def;  uv %n(c)=prof(m,i)
                i=Key_Ind('UW', keys,nks); prof(m,0)=uw_def;  uw %n(c)=prof(m,i)
                i=Key_Ind('VW', keys,nks); prof(m,0)=vw_def;  vw %n(c)=prof(m,i)
                i=Key_Ind('F22',keys,nks); prof(m,0)=f22_def; f22%n(c)=prof(m,i)
                i=Key_Ind('EPS',keys,nks); prof(m,0)=eps_def; eps%n(c)=prof(m,i)
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

                  if(turbulence_model == K_EPS_V2      .or.  &
                     turbulence_model == K_EPS_ZETA_F  .or.  &
                     turbulence_model == HYBRID_K_EPS_ZETA_F) then
                    prof(m,   0) = kin_def 
                    prof(m+1, 0) = kin_def
                    i = Key_Ind('KIN',keys,nks); 
                    kin % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                    prof(m,   0) = eps_def 
                    prof(m+1, 0) = eps_def
                    i = Key_Ind('EPS',keys,nks); 
                    eps % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
                    prof(m,   0) = v2_def 
                    prof(m+1, 0) = v2_def
                    i = Key_Ind('V2',keys,nks); 
                    v2 % n(c) = wi * prof(m, i) + (1.-wi) * prof(m+1, i)
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
