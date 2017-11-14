!==============================================================================!
  subroutine Load_Boundary_Conditions(grid, in_out)
!------------------------------------------------------------------------------!
!   Reads: name.b                                                              !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use par_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: in_out
!----------------------------------[Calling]-----------------------------------!
  real :: Distance
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, n, n_points, n_initial_cond, s
  integer           :: m, c1, c2, bc, mt, i
  character(len=80) :: name_bou, name_prof(128), dir, bc_name, mt_name
  integer           :: typBou(128)
  real              :: xyz(10024)
  real              :: wi
  real              :: x1(55555), x2(55555), Mres
  logical           :: here
!==============================================================================!

  !--------------------------------------------!
  !   Read the file with boundary conditions   !
  !--------------------------------------------!
  name_bou = name
  name_bou(len_trim(name)+1:len_trim(name)+2) = '.b'
  open(9, file=name_bou)
  if(this_proc < 2) write(*,*) '# Now reading the file:', name_bou

  !-------------------------!
  !   Phisical properties   !
  !-------------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) grid % n_materials
  do mt = 1,grid % n_materials

    call ReadC(9,inp,tn,ts,te)
    call To_Upper_Case(  inp(ts(1):te(1))  )
    call To_Upper_Case(  inp(ts(2):te(2))  )
    read(inp(ts(1):te(1)),*) mt_name

    ! Find material index
    do i=1, grid % n_materials 
      if(mt_name == grid % materials(i) % name) n=i      
    end do

    if( inp(ts(2):te(2))  ==  'FLUID') then 
      StateMat(n)=FLUID
    else if( inp(ts(2):te(2))  ==  'SOLID') then 
      StateMat(n)=SOLID
    else 
      if(this_proc < 2) write(*,*) 'Load_Boundary_Conditions: Unknown material state'
      stop  
    end if
    read(inp(ts(3):te(3)),*) VISc
    read(inp(ts(4):te(4)),*) DENc(n)
    if(HOT==YES) read(inp(ts(5):te(5)),*) CONc(n)
    if(HOT==YES) read(inp(ts(6):te(6)),*) CAPc(n)
  end do
  
  !-----------------------------------------------------!
  !   Boundary conditions 1 - read them from the file   !
  !-----------------------------------------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) grid % n_boundary_conditions

  do bc = 1,grid % n_boundary_conditions  ! number of boundary conditions

    call ReadC(9,inp,tn,ts,te)
    call To_Upper_Case(  inp(ts(1):te(1))  )
    call To_Upper_Case(  inp(ts(2):te(2))  )
    call To_Upper_Case(  inp(ts(3):te(3))  )
    read(inp(ts(1):te(1)),*) bc_name

    ! Find b.c. index
    do i=1, grid % n_boundary_conditions 
      if(bc_name == grid % boundary_conditions(i) % name) n=i      
    end do

    if( inp(ts(2):te(2)) == 'INFLOW') then 
      typBou(n)=INFLOW
      PER_BC = NO
    else if( inp(ts(2):te(2)) == 'WALL') then 
      typBou(n)=WALL
    else if( inp(ts(2):te(2)) == 'OUTFLOW') then 
      typBou(n)=OUTFLOW
    else if( inp(ts(2):te(2)) == 'SYMMETRY') then 
      typBou(n)=SYMMETRY
    else if( inp(ts(2):te(2)) == 'WALLFLUX') then 
      typBou(n)=WALLFL
    else if( inp(ts(2):te(2)) == 'CONVECTIVE') then 
      typBou(n)=CONVECT
    else if( inp(ts(2):te(2)) == 'PRESSURE') then 
      typBou(n)=PRESSURE
    else
      if(this_proc < 2) write(*,*) 'Load_Boundary_Conditions: Unknown boundary condition type: ', inp(ts(2):te(2))
      stop  
    end if
    if( inp(ts(3):te(3))  ==  'FILE') then
      read(inp(ts(4):te(4)),'(A80)') name_prof(n)
      write(*,*) 'n =            ', n
      write(*,*) 'name_prof(n) = ', name_prof(n)
    else
      read(inp(ts(3):te(3)),*) U % bound(n)
      read(inp(ts(4):te(4)),*) V % bound(n)
      read(inp(ts(5):te(5)),*) W % bound(n)
      if(typBou(n)==PRESSURE) then
        read(inp(ts(6):te(6)),*) P % bound(n)
        if(HOT==YES) then 
          read(inp(ts(7):te(7)),*) T % bound(n)
          if(SIMULA==EBM.or.SIMULA==HJ) then
            read(inp(ts(8):te(8)),*)   uu % bound(n)
            read(inp(ts(9):te(9)),*)   vv % bound(n)
            read(inp(ts(10):te(10)),*) ww % bound(n)
            read(inp(ts(11):te(11)),*) uv % bound(n)
            read(inp(ts(12):te(12)),*) uw % bound(n)
            read(inp(ts(13):te(13)),*) vw % bound(n)
            read(inp(ts(14):te(14)),*) Eps% bound(n)
            if(SIMULA==EBM) read(inp(ts(15):te(15)),*) f22 % bound(n)
          end if
          if(SIMULA==K_EPS) then
            read(inp(ts(8):te(8)),*) Kin % bound(n)
            read(inp(ts(9):te(9)),*) Eps % bound(n)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
            read(inp(ts(8):te(8)),*) Kin % bound(n)
            read(inp(ts(9):te(9)),*) Eps % bound(n)
            read(inp(ts(10):te(10)),*) v_2 % bound(n)
            read(inp(ts(11):te(11)),*) f22 % bound(n)
          end if
          if(SIMULA == SPA_ALL) then
            read(inp(ts(8):te(8)),*) VIS % bound(n)
          end if
          if(SIMULA == DES_SPA) then
            read(inp(ts(8):te(8)),*) VIS % bound(n)
          end if
        else  ! HOT .ne. YES
          if(SIMULA==EBM.or.SIMULA==HJ) then
            read(inp(ts(7):te(7)),*)   uu % bound(n)
            read(inp(ts(8):te(8)),*)   vv % bound(n)
            read(inp(ts(9):te(9)),*) ww % bound(n)
            read(inp(ts(10):te(10)),*) uv % bound(n)
            read(inp(ts(11):te(11)),*) uw % bound(n)
            read(inp(ts(12):te(12)),*) vw % bound(n)
            read(inp(ts(13):te(13)),*) Eps% bound(n)
            if(SIMULA==EBM) read(inp(ts(14):te(14)),*) f22 % bound(n)
          end if
          if(SIMULA==K_EPS) then
            read(inp(ts(7):te(7)),*) Kin % bound(n)
            read(inp(ts(8):te(8)),*) Eps % bound(n)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
            read(inp(ts(7):te(7)),*) Kin % bound(n)
            read(inp(ts(8):te(8)),*) Eps % bound(n)
            read(inp(ts(9):te(9)),*) v_2  % bound(n)
            read(inp(ts(10):te(10)),*) f22 % bound(n)
          end if
          if(SIMULA == SPA_ALL) then
            read(inp(ts(7):te(7)),*) VIS % bound(n)
          end if
          if(SIMULA == DES_SPA) then
            read(inp(ts(7):te(7)),*) VIS % bound(n)
          end if
        end if  ! HOT == YES
        name_prof(n)=''
      else   ! typBou .ne. PRESSURE
        if(HOT==YES) then 
          read(inp(ts(6):te(6)),*) T % bound(n)
          if(SIMULA==EBM.or.SIMULA==HJ) then
            read(inp(ts(7):te(7)),*) uu % bound(n)
            read(inp(ts(8):te(8)),*) vv % bound(n)
            read(inp(ts(9):te(9)),*) ww % bound(n)
            read(inp(ts(10):te(10)),*) uv % bound(n)
            read(inp(ts(11):te(11)),*) uw % bound(n)
            read(inp(ts(12):te(12)),*) vw % bound(n)
            read(inp(ts(13):te(13)),*) Eps% bound(n)
            if(SIMULA==EBM) read(inp(ts(14):te(14)),*) f22 % bound(n)
          end if
          if(SIMULA==K_EPS) then
            read(inp(ts(7):te(7)),*) Kin % bound(n)
            read(inp(ts(8):te(8)),*) Eps % bound(n)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
            read(inp(ts(7):te(7)),*) Kin % bound(n)
            read(inp(ts(8):te(8)),*) Eps % bound(n)
            read(inp(ts(9):te(9)),*) v_2 % bound(n)
            read(inp(ts(10):te(10)),*) f22 % bound(n)
          end if
          if(SIMULA == SPA_ALL) then
            read(inp(ts(7):te(7)),*) VIS % bound(n)
          end if
          if(SIMULA == DES_SPA) then
            read(inp(ts(7):te(7)),*) VIS % bound(n)
          end if
        else  ! HOT .ne. YES
          if(SIMULA==EBM.or.SIMULA==HJ) then
            read(inp(ts(6):te(6)),*) uu % bound(n)
            read(inp(ts(7):te(7)),*) vv % bound(n)
            read(inp(ts(8):te(8)),*) ww % bound(n)
            read(inp(ts(9):te(9)),*) uv % bound(n)
            read(inp(ts(10):te(10)),*) uw % bound(n)
            read(inp(ts(11):te(11)),*) vw % bound(n)
            read(inp(ts(12):te(12)),*) Eps% bound(n)
            if(SIMULA==EBM) read(inp(ts(13):te(13)),*) f22 % bound(n)
          end if
          if(SIMULA==K_EPS) then
            read(inp(ts(6):te(6)),*) Kin % bound(n)
            read(inp(ts(7):te(7)),*) Eps % bound(n)
          end if
          if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
            read(inp(ts(6):te(6)),*) Kin % bound(n)
            read(inp(ts(7):te(7)),*) Eps % bound(n)
            read(inp(ts(8):te(8)),*) v_2  % bound(n)
            read(inp(ts(9):te(9)),*) f22 % bound(n)
          end if
          if(SIMULA == SPA_ALL) then
            read(inp(ts(6):te(6)),*) VIS % bound(n)
          end if
          if(SIMULA == DES_SPA) then
            read(inp(ts(6):te(6)),*) VIS % bound(n)
          end if
        end if  ! HOT == YES
        name_prof(n)=''
      end if  ! typBou == PRESSURE
    end if    ! inp .not. file
  end do      

  !------------------------!
  !   Initial conditions   !
  !------------------------!
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) n_initial_cond
  write(*,*) '# Number of initial conditions: ', n_initial_cond
  if(n_initial_cond > grid % n_materials) then
    if(this_proc < 2) write(*,*) 'Warning: there are more initial conditions then materials'
  end if

  do n=1,n_initial_cond
    call ReadC(9,inp,tn,ts,te)
    call To_Upper_Case(inp(ts(2):te(2)))

    ! Initial conditions given in GMV file
    if(inp(ts(2):te(2)) == 'FILE') then
      read(inp(ts(3):te(3)),'(A80)') namIni(n)
      write(*,*) '@Load_Boundary_Conditions: material ', n, '; init. cond. given by file: ', namIni(n)
    else
      namIni(n) = ''

      ! Initial conditions given by constant
      read(inp(ts(2):te(2)),*) U % init(n)
      read(inp(ts(3):te(3)),*) V % init(n)
      read(inp(ts(4):te(4)),*) W % init(n)
 
      if(HOT==YES) then
        read(inp(ts(5):te(5)),*) T % init(n)
        if(SIMULA==EBM.or.SIMULA==HJ) then
          read(inp(ts(6):te(6)),*) uu % init(n)
          read(inp(ts(7):te(7)),*) vv % init(n)
          read(inp(ts(8):te(8)),*) ww % init(n)
          read(inp(ts(9):te(9)),*) uv % init(n)
          read(inp(ts(10):te(10)),*) uw % init(n)
          read(inp(ts(11):te(11)),*) vw % init(n)
          read(inp(ts(12):te(12)),*) Eps% init(n)
          if(SIMULA==EBM) read(inp(ts(13):te(13)),*) f22 % init(n)
        end if
        if(SIMULA==K_EPS) then
          read(inp(ts(6):te(6)),*) Kin % init(n)
          read(inp(ts(7):te(7)),*) Eps % init(n)
        end if
        if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
          read(inp(ts(6):te(6)),*) Kin % init(n)
          read(inp(ts(7):te(7)),*) Eps % init(n)
          read(inp(ts(8):te(8)),*) v_2  % init(n)
          read(inp(ts(9):te(9)),*) f22 % init(n)
        end if
        if(SIMULA == SPA_ALL) then
          read(inp(ts(6):te(6)),*) VIS % init(n)
        end if
        if(SIMULA == DES_SPA) then
          read(inp(ts(6):te(6)),*) VIS % init(n)
        end if
      else ! HOT /= YES
        if(SIMULA==EBM.or.SIMULA==HJ) then
          read(inp(ts(5):te(5)),*) uu % init(n)
          read(inp(ts(6):te(6)),*) vv % init(n)
          read(inp(ts(7):te(7)),*) ww % init(n)
          read(inp(ts(8):te(8)),*) uv % init(n)
          read(inp(ts(9):te(9)),*) uw % init(n)
          read(inp(ts(10):te(10)),*) vw % init(n)
          read(inp(ts(11):te(11)),*) Eps% init(n)
          if(SIMULA==EBM) read(inp(ts(12):te(12)),*) f22 % init(n)
        end if
        if(SIMULA==K_EPS) then
          read(inp(ts(5):te(5)),*) Kin % init(n)
          read(inp(ts(6):te(6)),*) Eps % init(n)
        end if
        if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
          read(inp(ts(5):te(5)),*) Kin % init(n)
          read(inp(ts(6):te(6)),*) Eps % init(n)
          read(inp(ts(7):te(7)),*) v_2  % init(n)
          read(inp(ts(8):te(8)),*) f22 % init(n)
        end if
        if(SIMULA == SPA_ALL) then
          read(inp(ts(5):te(5)),*) VIS % init(n)
        end if
        if(SIMULA == DES_SPA) then
          read(inp(ts(5):te(5)),*) VIS % init(n)
        end if
      end if
    end if
  end do  

  close(9)

  !----------------------------------------------------------------------!
  !   Boundary conditions 2 - distribute them over computational cells   !
  !----------------------------------------------------------------------!
  do n=1,grid % n_boundary_conditions

    ! Boundary condition is given by a single constant
    if(name_prof(n) == '') then 
      do c = -1,-grid % n_bnd_cells,-1
        if(bcmark(c) == n) then
          TypeBC(c) = typBou(n)

          ! If in_out is set to true, set boundary values,
          ! otherwise, just the TypeBC remains set.
          if(in_out) then
            U % n(c) = U % bound(n) 
            V % n(c) = V % bound(n)
            W % n(c) = W % bound(n)
            P % n(c) = P % bound(n) 
            if(HOT == YES) then
              if(TypeBC(c).eq.WALLFL) then
                T % q(c) =  T % bound(n)
              else
                T % n(c) =  T % bound(n)
              endif
            end if  ! for HOT==YES
            if(SIMULA==EBM.or.SIMULA==HJ) then
              uu % n(c) = uu % bound(n)
              vv % n(c) = vv % bound(n)
              ww % n(c) = ww % bound(n)
              uv % n(c) = uv % bound(n)
              uw % n(c) = uw % bound(n)
              vw % n(c) = vw % bound(n)
              Eps % n(c) = Eps % bound(n)
              if(SIMULA==EBM) f22 % n(c)   = f22 % bound(n)
            end if
            if(SIMULA==K_EPS) then
              Kin % n(c) = Kin % bound(n)
              Eps % n(c) = Eps % bound(n)
              Uf(c)        = 0.047
              Ynd(c)       = 30.0
            end if
            if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
              Kin % n(c)   = Kin % bound(n)
              Eps % n(c)   = Eps % bound(n)
              f22 % n(c)   = f22 % bound(n)
              v_2 % n(c)   = v_2 % bound(n)
            end if
            if(SIMULA == SPA_ALL) then
              VIS % n(c)   = VIS % bound(n)
            end if
            if(SIMULA == DES_SPA) then
              VIS % n(c)   = VIS % bound(n)
            end if
          end if
        end if 
      end do

    ! Boundary condition is prescribed in a file 
    else
      open(9, file=name_prof(n))
      if(this_proc < 2) write(*,*) '# Now reading the file:', name_prof(n)
      call ReadC(9,inp,tn,ts,te)
      read(inp(ts(1):te(1)),*) n_points                  ! number of points
      call ReadC(9,inp,tn,ts,te)
      read(inp(ts(1):te(1)),*) dir  ! direction
      call To_Upper_Case(dir)
      if(dir=="XPL" .or. dir=="YPL" .or. dir=="ZPL") then
        do m=1,n_points
          call ReadC(9,inp,tn,ts,te)
          read(inp(ts(1):te(1)),*) x1(m)
          read(inp(ts(2):te(2)),*) x2(m)
          read(inp(ts(3):te(3)),*) U % pro(m)
          read(inp(ts(4):te(4)),*) V % pro(m)
          read(inp(ts(5):te(5)),*) W % pro(m)
          if(SIMULA == EBM) then
            read(inp(ts(6):te(6)),*) uu % pro(m)
            read(inp(ts(7):te(7)),*) vv % pro(m)
            read(inp(ts(8):te(8)),*) ww % pro(m)
            read(inp(ts(9):te(9)),*) uv % pro(m)
            read(inp(ts(10):te(10)),*) uw % pro(m)
            read(inp(ts(11):te(11)),*) vw % pro(m)
            read(inp(ts(12):te(12)),*) f22 % pro(m)
            read(inp(ts(13):te(13)),*) Eps % pro(m)
          end if
        end do  

        ! Set the closest point
        do c = -1,-grid % n_bnd_cells,-1
          if(bcmark(c) == n) then
            TypeBC(c) = typBou(n)
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
              U%n(c) = U % pro(c1)
              V%n(c) = V % pro(c1)
              W%n(c) = W % pro(c1)
              if(HOT==YES) T%n(c) = T%pro(c1)
              if(SIMULA==K_EPS) then
                Kin%n(c) = Kin%pro(c1)
                Eps%n(c) = Eps%pro(c1)
              end if
              if(SIMULA==K_EPS_VV.or.SIMULA==ZETA) then
                Kin%n(c) = Kin%pro(c1)
                Eps%n(c) = Eps%pro(c1)
                v_2%n(c) = v_2%pro(c1)
                f22%n(c) = f22%pro(c1)
              end if
              if(SIMULA == DES_SPA) then
                VIS%n(c) = VIS%pro(c1)
              end if
              if(SIMULA == EBM) then
                uu%n(c) = uu % pro(c1)
                vv%n(c) = vv % pro(c1)
                ww%n(c) = ww % pro(c1)
                uv%n(c) = uv % pro(c1)
                uw%n(c) = uw % pro(c1)
                vw%n(c) = vw % pro(c1)
                f22%n(c) = f22 % pro(c1)
                Eps%n(c) = Eps % pro(c1)
              end if        
            end if    !end if(in_out)
          end if      !end if(bcmark(c) == n)
        end do        !end do c = -1,-grid % n_bnd_cells,-1
      else  ! dir == "XPL" ...
        do m=1,n_points
          call ReadC(9,inp,tn,ts,te)
          read(inp(ts(1):te(1)),*) xyz(m)
          read(inp(ts(2):te(2)),*) U % pro(m)
          read(inp(ts(3):te(3)),*) V % pro(m)
          read(inp(ts(4):te(4)),*) W % pro(m)
          if(HOT==YES) then
            read(inp(ts(5):te(5)),*) T % pro(m)
            if(SIMULA==K_EPS) then
              read(inp(ts(6):te(6)),*) Kin % pro(m)
              read(inp(ts(7):te(7)),*) Eps % pro(m)
            end if
            if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
              read(inp(ts(6):te(6)),*) Kin % pro(m)
              read(inp(ts(7):te(7)),*) Eps % pro(m)
              read(inp(ts(8):te(8)),*) v_2 % pro(m)
              read(inp(ts(9):te(9)),*) f22 % pro(m)
            end if
            if(SIMULA == SPA_ALL) then
              read(inp(ts(6):te(6)),*) VIS % pro(m)
            end if
            if(SIMULA == DES_SPA) then
              read(inp(ts(6):te(6)),*) VIS % pro(m)
            end if
          else
            if(SIMULA==K_EPS) then
              read(inp(ts(5):te(5)),*) Kin % pro(m)
              read(inp(ts(6):te(6)),*) Eps % pro(m)
            end if
            if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
              read(inp(ts(5):te(5)),*) Kin % pro(m)
              read(inp(ts(6):te(6)),*) Eps % pro(m)
              read(inp(ts(7):te(7)),*) v_2 % pro(m)
              read(inp(ts(8):te(8)),*) f22 % pro(m)
            end if
            if(SIMULA == SPA_ALL) then
              read(inp(ts(5):te(5)),*) VIS % pro(m)
            end if
            if(SIMULA == DES_SPA) then
              read(inp(ts(5):te(5)),*) VIS % pro(m)
            end if
            if(SIMULA == EBM) then
              read(inp(ts(5):te(5)),*) uu % pro(m)
              read(inp(ts(6):te(6)),*) vv % pro(m)
              read(inp(ts(7):te(7)),*) ww % pro(m)
              read(inp(ts(8):te(8)),*) uv % pro(m)
              read(inp(ts(9):te(9)),*) uw % pro(m)
              read(inp(ts(10):te(10)),*) vw % pro(m)
              read(inp(ts(11):te(11)),*) f22% pro(m)
              read(inp(ts(12):te(12)),*) Eps% pro(m)
            end if
          end if  
        end do
           
        do c = -1,-grid % n_bnd_cells,-1
          if(bcmark(c) == n) then
            TypeBC(c) = typBou(n)
          
            ! If in_out is set to true, set boundary values,
            ! otherwise, just the TypeBC remains set.
            if(in_out) then
              do m=1,n_points-1
                here = .FALSE. 

                ! Compute the weight factors
                if( (dir == 'X' .or. dir == 'x') .and.                  &
                   grid % xc(c) >= xyz(m) .and. grid % xc(c) <= xyz(m+1) ) then
                  wi = ( xyz(m+1)-grid % xc(c) ) / ( xyz(m+1) - xyz(m) )
                  here = .TRUE.
                else if( (dir == 'Y' .or. dir == 'y') .and.             &
                     grid % yc(c) >= xyz(m) .and. grid % yc(c) <= xyz(m+1) ) then
                  wi = ( xyz(m+1)-grid % yc(c) ) / ( xyz(m+1) - xyz(m) )
                    here = .TRUE.
                else if( (dir == 'Z' .or. dir == 'z') .and.             &
                     grid % zc(c) >= xyz(m) .and. grid % zc(c) <= xyz(m+1) ) then
                  wi = ( xyz(m+1)-grid % zc(c) ) / ( xyz(m+1) - xyz(m) )
                  here = .TRUE.
                else if( (dir == 'RX' .or. dir == 'rx') .and.           &
                       sqrt(grid % yc(c)*grid % yc(c)+grid % zc(c)*grid % zc(c)) >= xyz(m) .and.      &
                     sqrt(grid % yc(c)*grid % yc(c)+grid % zc(c)*grid % zc(c)) <= xyz(m+1) ) then
                  wi = ( xyz(m+1) - sqrt(grid % yc(c)*grid % yc(c)+grid % zc(c)*grid % zc(c)) )     &
                     / ( xyz(m+1) - xyz(m) )
                  here = .TRUE.
                else if( (dir == 'RY' .or. dir == 'ry') .and.           &
                     sqrt(grid % xc(c)*grid % xc(c)+grid % zc(c)*grid % zc(c)) >= xyz(m) .and.      &
                     sqrt(grid % xc(c)*grid % xc(c)+grid % zc(c)*grid % zc(c)) <= xyz(m+1) ) then
                  wi = ( xyz(m+1) - sqrt(grid % xc(c)*grid % xc(c)+grid % zc(c)*grid % zc(c)) )     &
                     / ( xyz(m+1) - xyz(m) )
                  here = .TRUE.
                else if( (dir == 'RZ' .or. dir == 'rz') .and.           &
                     sqrt(grid % xc(c)*grid % xc(c)+grid % yc(c)*grid % yc(c)) <= xyz(m) .and.      &
                     sqrt(grid % xc(c)*grid % xc(c)+grid % yc(c)*grid % yc(c)) >= xyz(m+1) ) then
                    wi = ( xyz(m+1) - sqrt(grid % xc(c)*grid % xc(c)+grid % yc(c)*grid % yc(c)) )     &
                     / ( xyz(m+1) - xyz(m) )
                  here = .TRUE.
                end if

                ! Interpolate the profiles     
                if(here) then
                  U % n(c) = wi*U % pro(m) + (1.-wi)*U % pro(m+1)
                  V % n(c) = wi*V % pro(m) + (1.-wi)*V % pro(m+1)
                  W % n(c) = wi*W % pro(m) + (1.-wi)*W % pro(m+1)
                  if(HOT==YES) &
                    T % n(c) = wi*T % pro(m) + (1.-wi)*T % pro(m+1)
                  if(SIMULA==K_EPS) then
                    Kin % n(c) = wi*Kin % pro(m) + (1.-wi)*Kin % pro(m+1)
                    Eps % n(c) = wi*Eps % pro(m) + (1.-wi)*Eps % pro(m+1)
                  end if
                  if(SIMULA==K_EPS_VV.or.SIMULA == ZETA.or.SIMULA == HYB_ZETA) then
                    Kin % n(c) = wi*Kin % pro(m) + (1.-wi)*Kin % pro(m+1)
                    Eps % n(c) = wi*Eps % pro(m) + (1.-wi)*Eps % pro(m+1)
                    f22 % n(c) = wi*f22 % pro(m) + (1.-wi)*f22 % pro(m+1)
                    v_2 % n(c) = wi*v_2 % pro(m)  + (1.-wi)*v_2% pro(m+1)
                  end if
                  if(SIMULA == SPA_ALL) then
                    VIS % n(c) = wi*VIS % pro(m) + (1.-wi)*VIS % pro(m+1)
                  end if
                  if(SIMULA == DES_SPA) then
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

  !-------------------------------------!
  !   Finally handle the buffer cells   !
  !-------------------------------------!
  do c = -1,-grid % n_bnd_cells,-1
    if(bcmark(c) == BUFFER) TypeBC(c)=BUFFER 
  end do

  end subroutine
