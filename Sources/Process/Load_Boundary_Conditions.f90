!==============================================================================!
  subroutine Load_Boundary_Conditions(grid, in_out)
!------------------------------------------------------------------------------!
!   Reads: .bnd file                                                           !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Flow_Mod
  use rans_mod
  use par_mod
  use Tokenizer_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: in_out
!----------------------------------[Calling]-----------------------------------!
  real :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n, n_points, n_initial_cond, s
  integer           :: m, c1, bc, mt, i
  character(len=80) :: name_bou, name_prof(128), dir, bc_name, mt_name
  integer           :: type_bnd_cond(128)
  real              :: xyz(8192)
  real              :: wi
  real              :: x1(65536), x2(65536), Mres
  logical           :: here
  character         :: name_ini(128)*80
!==============================================================================!

  !--------------------------------------------!
  !   Read the file with boundary conditions   !
  !--------------------------------------------!
  name_bou = problem_name
  name_bou(len_trim(problem_name)+1:len_trim(problem_name)+4) = '.bnd'
  open(9, file=name_bou)
  if(this_proc < 2) print *, '# Reading the file: ', name_bou

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
        name_prof(n)=''
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
        name_prof(n)=''
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

  !----------------------------------------------------------------------!
  !   Boundary conditions 2 - distribute them over computational cells   !
  !----------------------------------------------------------------------!
  do n = 1, grid % n_bnd_cond

    grid % bnd_cond % type(n) = type_bnd_cond(n)

    ! Boundary condition is given by a single constant
    if(name_prof(n) == '') then 
      do c = -1,-grid % n_bnd_cells,-1
        if(grid % bnd_cond % color(c) == n) then
!new:   if(Grid_Mod_Bnd_Cond_Type(grid,c) == type_bnd_cond(n)) then

          ! If in_out is set to true, set boundary values,
          ! otherwise, just the ypTeBC remains set.
          if(in_out) then
            u % n(c) = u % bound(n) 
            v % n(c) = v % bound(n)
            w % n(c) = w % bound(n)
            P % n(c) = P % bound(n) 
            if(heat_transfer == YES) then
              if(type_bnd_cond(n) .eq. WALLFL) then
                T % q(c) =  T % bound(n)
              else
                T % n(c) =  T % bound(n)
              endif
            end if  ! for heat_transfer == YES
            if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
               turbulence_model == HANJALIC_JAKIRLIC) then
              uu % n(c) = uu % bound(n)
              vv % n(c) = vv % bound(n)
              ww % n(c) = ww % bound(n)
              uv % n(c) = uv % bound(n)
              uw % n(c) = uw % bound(n)
              vw % n(c) = vw % bound(n)
              eps % n(c) = eps % bound(n)
              if(turbulence_model == REYNOLDS_STRESS_MODEL) then
                f22 % n(c)   = f22 % bound(n)
              end if
            end if
            if(turbulence_model == K_EPS) then
              kin % n(c) = kin % bound(n)
              eps % n(c) = eps % bound(n)
              Uf(c)        = 0.047
              Ynd(c)       = 30.0
            end if
            if(turbulence_model == K_EPS_V2 .or.  &
               turbulence_model == K_EPS_ZETA_F     .or.  &
               turbulence_model == HYBRID_K_EPS_ZETA_F) then
              kin % n(c)   = kin % bound(n)
              eps % n(c)   = eps % bound(n)
              f22 % n(c)   = f22 % bound(n)
              v2  % n(c)   = v2  % bound(n)
            end if
            if(turbulence_model == SPALART_ALLMARAS) then
              VIS % n(c)   = VIS % bound(n)
            end if
            if(turbulence_model == DES_SPALART) then
              VIS % n(c)   = VIS % bound(n)
            end if
          end if
        end if 
      end do

    ! Boundary condition is prescribed in a file 
    else
      open(9, file=name_prof(n))
      if(this_proc < 2) print *, '# Reading the file: ', name_prof(n)
      call Tokenizer_Mod_Read_Line(9)
      read(line % tokens(1),*) n_points                  ! number of points
      call Tokenizer_Mod_Read_Line(9)
      read(line % tokens(1),*) dir  ! direction
      call To_Upper_Case(dir)
      if(dir=="XPL" .or. dir=="YPL" .or. dir=="ZPL") then
        do m=1,n_points
          call Tokenizer_Mod_Read_Line(9)
          read(line % tokens(1),*) x1(m)
          read(line % tokens(2),*) x2(m)
          read(line % tokens(3),*) u % pro(m)
          read(line % tokens(4),*) v % pro(m)
          read(line % tokens(5),*) w % pro(m)
          if(turbulence_model == REYNOLDS_STRESS_MODEL) then
            read(line % tokens(6),*)  uu % pro(m)
            read(line % tokens(7),*)  vv % pro(m)
            read(line % tokens(8),*)  ww % pro(m)
            read(line % tokens(9),*)  uv % pro(m)
            read(line % tokens(10),*) uw % pro(m)
            read(line % tokens(11),*) vw % pro(m)
            read(line % tokens(12),*) f22 % pro(m)
            read(line % tokens(13),*) eps % pro(m)
          end if
        end do  

        ! Set the closest point
        do c = -1,-grid % n_bnd_cells,-1
          if(grid % bnd_cond % color(c) == n) then
!new:     if(Grid_Mod_Bnd_Cond_Type(grid,c) == type_bnd_cond(n)) then
            if(in_out) then    !if .true. set boundary values, otherwise, just set TypeBC
              Mres = HUGE
              do s=1,n_points
                if(dir=="XPL") then
                  if(Distance(x1(s),x2(s),0.0,grid % yc(c),grid % zc(c),0.0) < Mres) then
                    Mres = Distance(x1(s),x2(s),0.0,grid % yc(c),grid % zc(c),0.0)
                    c1 = s
                  end if
                else if(dir=="YPL") then
                  if(Distance(x1(s),x2(s),0.0,grid % xc(c),grid % zc(c),0.0) < Mres) then
                    Mres = Distance(x1(s),x2(s),0.0,grid % xc(c),grid % zc(c),0.0)
                    c1 = s
                  end if
                else if(dir=="ZPL") then
                  if(Distance(x1(s),x2(s),0.0,grid % xc(c),grid % yc(c),0.0) < Mres) then
                    Mres = Distance(x1(s),x2(s),0.0,grid % xc(c),grid % yc(c),0.0)
                    c1 = s
                  end if
                end if
              end do
              u%n(c) = u % pro(c1)
              v%n(c) = v % pro(c1)
              w%n(c) = w % pro(c1)
              if(heat_transfer == YES) T%n(c) = T%pro(c1)
              if(turbulence_model == K_EPS) then
                kin%n(c) = kin%pro(c1)
                eps%n(c) = eps%pro(c1)
              end if
              if(turbulence_model == K_EPS_V2 .or.  &
                 turbulence_model == K_EPS_ZETA_F) then
                kin%n(c) = kin%pro(c1)
                eps%n(c) = eps%pro(c1)
                v2 %n(c) = v2 %pro(c1)
                f22%n(c) = f22%pro(c1)
              end if
              if(turbulence_model == DES_SPALART) then
                VIS%n(c) = VIS%pro(c1)
              end if
              if(turbulence_model == REYNOLDS_STRESS_MODEL) then
                uu%n(c) = uu % pro(c1)
                vv%n(c) = vv % pro(c1)
                ww%n(c) = ww % pro(c1)
                uv%n(c) = uv % pro(c1)
                uw%n(c) = uw % pro(c1)
                vw%n(c) = vw % pro(c1)
                f22%n(c) = f22 % pro(c1)
                eps%n(c) = eps % pro(c1)
              end if        
            end if    !end if(in_out)
          end if      !end if(grid % bnd_cond % color(c) == n)
        end do        !end do c = -1,-grid % n_bnd_cells,-1
      else  ! dir == "XPL" ...
        do m=1,n_points
          call Tokenizer_Mod_Read_Line(9)
          read(line % tokens(1),*) xyz(m)
          read(line % tokens(2),*) u % pro(m)
          read(line % tokens(3),*) v % pro(m)
          read(line % tokens(4),*) w % pro(m)
          if(heat_transfer == YES) then
            read(line % tokens(5),*) T % pro(m)
            if(turbulence_model == K_EPS) then
              read(line % tokens(6),*) kin % pro(m)
              read(line % tokens(7),*) eps % pro(m)
            end if
            if(turbulence_model == K_EPS_V2 .or.  &
               turbulence_model == K_EPS_ZETA_F     .or.  &
               turbulence_model == HYBRID_K_EPS_ZETA_F) then
              read(line % tokens(6),*) kin % pro(m)
              read(line % tokens(7),*) eps % pro(m)
              read(line % tokens(8),*) v2  % pro(m)
              read(line % tokens(9),*) f22 % pro(m)
            end if
            if(turbulence_model == SPALART_ALLMARAS) then
              read(line % tokens(6),*) VIS % pro(m)
            end if
            if(turbulence_model == DES_SPALART) then
              read(line % tokens(6),*) VIS % pro(m)
            end if
          else
            if(turbulence_model == K_EPS) then
              read(line % tokens(5),*) kin % pro(m)
              read(line % tokens(6),*) eps % pro(m)
            end if
            if(turbulence_model == K_EPS_V2 .or.  &
               turbulence_model == K_EPS_ZETA_F     .or.  &
               turbulence_model == HYBRID_K_EPS_ZETA_F) then
              read(line % tokens(5),*) kin % pro(m)
              read(line % tokens(6),*) eps % pro(m)
              read(line % tokens(7),*) v2  % pro(m)
              read(line % tokens(8),*) f22 % pro(m)
            end if
            if(turbulence_model == SPALART_ALLMARAS) then
              read(line % tokens(5),*) VIS % pro(m)
            end if
            if(turbulence_model == DES_SPALART) then
              read(line % tokens(5),*) VIS % pro(m)
            end if
            if(turbulence_model == REYNOLDS_STRESS_MODEL) then
              read(line % tokens(5),*) uu % pro(m)
              read(line % tokens(6),*) vv % pro(m)
              read(line % tokens(7),*) ww % pro(m)
              read(line % tokens(8),*) uv % pro(m)
              read(line % tokens(9),*) uw % pro(m)
              read(line % tokens(10),*) vw % pro(m)
              read(line % tokens(11),*) f22% pro(m)
              read(line % tokens(12),*) eps% pro(m)
            end if
          end if  
        end do
           
        do c = -1,-grid % n_bnd_cells,-1
          if(grid % bnd_cond % color(c) == n) then
!new:     if(Grid_Mod_Bnd_Cond_Type(grid,c) == type_bnd_cond(n)) then
          
            ! If in_out is set to true, set boundary values,
            ! otherwise, just the TypeBC remains set.
            if(in_out) then
              do m=1,n_points-1
                here = .false. 

                ! Compute the weight factors
                if( (dir == 'X' .or. dir == 'x') .and.                  &
                   grid % xc(c) >= xyz(m) .and. grid % xc(c) <= xyz(m+1) ) then
                  wi = ( xyz(m+1)-grid % xc(c) ) / ( xyz(m+1) - xyz(m) )
                  here = .true.
                else if( (dir == 'Y' .or. dir == 'y') .and.             &
                     grid % yc(c) >= xyz(m) .and. grid % yc(c) <= xyz(m+1) ) then
                  wi = ( xyz(m+1)-grid % yc(c) ) / ( xyz(m+1) - xyz(m) )
                    here = .true.
                else if( (dir == 'Z' .or. dir == 'z') .and.             &
                     grid % zc(c) >= xyz(m) .and. grid % zc(c) <= xyz(m+1) ) then
                  wi = ( xyz(m+1)-grid % zc(c) ) / ( xyz(m+1) - xyz(m) )
                  here = .true.
                else if( (dir == 'RX' .or. dir == 'rx') .and.           &
                     sqrt(grid % yc(c)*grid % yc(c)+grid % zc(c)*grid % zc(c)) >= xyz(m) .and.      &
                     sqrt(grid % yc(c)*grid % yc(c)+grid % zc(c)*grid % zc(c)) <= xyz(m+1) ) then
                  wi = ( xyz(m+1) - sqrt(grid % yc(c)*grid % yc(c)+grid % zc(c)*grid % zc(c)) )     &
                     / ( xyz(m+1) - xyz(m) )
                  here = .true.
                else if( (dir == 'RY' .or. dir == 'ry') .and.           &
                     sqrt(grid % xc(c)*grid % xc(c)+grid % zc(c)*grid % zc(c)) >= xyz(m) .and.      &
                     sqrt(grid % xc(c)*grid % xc(c)+grid % zc(c)*grid % zc(c)) <= xyz(m+1) ) then
                  wi = ( xyz(m+1) - sqrt(grid % xc(c)*grid % xc(c)+grid % zc(c)*grid % zc(c)) )     &
                     / ( xyz(m+1) - xyz(m) )
                  here = .true.
                else if( (dir == 'RZ' .or. dir == 'rz') .and.           &
                     sqrt(grid % xc(c)*grid % xc(c)+grid % yc(c)*grid % yc(c)) <= xyz(m) .and.      &
                     sqrt(grid % xc(c)*grid % xc(c)+grid % yc(c)*grid % yc(c)) >= xyz(m+1) ) then
                    wi = ( xyz(m+1) - sqrt(grid % xc(c)*grid % xc(c)+grid % yc(c)*grid % yc(c)) )   &
                     / ( xyz(m+1) - xyz(m) )
                  here = .true.
                end if

                ! Interpolate the profiles     
                if(here) then
                  u % n(c) = wi*u % pro(m) + (1.-wi)*u % pro(m+1)
                  v % n(c) = wi*v % pro(m) + (1.-wi)*v % pro(m+1)
                  w % n(c) = wi*w % pro(m) + (1.-wi)*w % pro(m+1)
                  if(heat_transfer == YES)  &
                    T % n(c) = wi*T % pro(m) + (1.-wi)*T % pro(m+1)
                  if(turbulence_model == K_EPS) then
                    kin % n(c) = wi*kin % pro(m) + (1.-wi)*kin % pro(m+1)
                    eps % n(c) = wi*eps % pro(m) + (1.-wi)*eps % pro(m+1)
                  end if
                  if(turbulence_model == K_EPS_V2      .or.  &
                     turbulence_model == K_EPS_ZETA_F  .or.  &
                     turbulence_model == HYBRID_K_EPS_ZETA_F) then
                    kin % n(c) = wi*kin % pro(m) + (1.-wi)*kin % pro(m+1)
                    eps % n(c) = wi*eps % pro(m) + (1.-wi)*eps % pro(m+1)
                    f22 % n(c) = wi*f22 % pro(m) + (1.-wi)*f22 % pro(m+1)
                    v2  % n(c) = wi*v2  % pro(m)  + (1.-wi)*v2 % pro(m+1)
                  end if
                  if(turbulence_model == SPALART_ALLMARAS) then
                    VIS % n(c) = wi*VIS % pro(m) + (1.-wi)*VIS % pro(m+1)
                  end if
                  if(turbulence_model == DES_SPALART) then
                    VIS % n(c) = wi*VIS % pro(m) + (1.-wi)*VIS % pro(m+1)
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
