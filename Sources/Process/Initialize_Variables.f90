!==============================================================================!
  subroutine Initialize_Variables(grid)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use par_mod
  use rans_mod
  use Grid_Mod
  use Bulk_Mod
  use Parameters_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------!
  real    :: Stot
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, m, s, n
  integer :: n_wall, n_inflow, n_outflow, n_symmetry, n_heated_wall, n_convect
!==============================================================================!

  Area  = 0.0
  write(*,*) 'grid % n_materials: ', grid % n_materials
  do n = 1, grid % n_materials
    do c = 1, grid % n_cells
      U % n(c)    = U % init(material(c))
      U % o(c)    = U % init(material(c)) 
      U % oo(c)   = U % init(material(c))
      V % n(c)    = V % init(material(c)) 
      V % o(c)    = V % init(material(c))
      V % oo(c)   = V % init(material(c))
      W % n(c)    = W % init(material(c)) 
      W % o(c)    = W % init(material(c))
      W % oo(c)   = W % init(material(c))
      if(HOT==YES) then
        T % n(c)  = T % init(material(c)) 
        T % o(c)  = T % init(material(c)) 
        T % oo(c) = T % init(material(c)) 
      end if 
      if(SIMULA==EBM.or.SIMULA==HJ) then
        uu % n(c)  = uu % init(material(c))
        vv % n(c)  = vv % init(material(c))
        ww % n(c)  = ww % init(material(c))
        uv % n(c)  = uv % init(material(c))
        uw % n(c)  = uw % init(material(c))
        vw % n(c)  = vw % init(material(c))
        Eps % n(c) = Eps % init(material(c))
        if(SIMULA==EBM) f22 % n(c) = f22 % init(material(c))
      end if
      if(SIMULA==K_EPS.or.SIMULA==HYB_PITM) then
        Kin % n(c)  = Kin % init(material(c))
        Kin % o(c)  = Kin % init(material(c))
        Kin % oo(c) = Kin % init(material(c))
        Eps % n(c)  = Eps % init(material(c))
        Eps % o(c)  = Eps % init(material(c))
        Eps % oo(c) = Eps % init(material(c))
        Uf(c)       = 0.047
        Ynd(c)      = 30.0
      end if
      if(SIMULA==K_EPS_VV.or.SIMULA==ZETA.or.SIMULA==HYB_ZETA) then
        Kin % n(c)  = Kin % init(material(c))
        Kin % o(c)  = Kin % init(material(c))
        Kin % oo(c) = Kin % init(material(c))
        Eps % n(c)  = Eps % init(material(c))
        Eps % o(c)  = Eps % init(material(c))
        Eps % oo(c) = Eps % init(material(c))
        f22 % n(c)  = f22 % init(material(c))
        f22 % o(c)  = f22 % init(material(c))
        f22 % oo(c) = f22 % init(material(c))
        v_2 % n(c)  = v_2 % init(material(c))
        v_2 % o(c)  = v_2 % init(material(c))
        v_2 % oo(c) = v_2 % init(material(c))
        Uf(c)       = 0.047
        Ynd(c)      = 30.0
      end if
      if(SIMULA == SPA_ALL .or. SIMULA==DES_SPA) then      
        VIS % n(c)  = VIS % init(material(c))
        VIS % o(c)  = VIS % init(material(c))
        VIS % oo(c) = VIS % init(material(c))
      end if
      if(SIMULA==HYB_ZETA.or.SIMULA==DES_SPA.or.&
      SIMULA==DNS.or.SIMULA==HYB_PITM.or.URANS==YES) then
        U % mean(c) = 0.0
        V % mean(c) = 0.0
        W % mean(c) = 0.0
      end if
    end do 
  end do   !end do n=1,grid % n_materials

  if(TGV == YES) then
    do c = 1, grid % n_cells
      U % n(c)  = -sin(grid % xc(c))*cos(grid % yc(c))
      U % o(c)  = -sin(grid % xc(c))*cos(grid % yc(c))
      U % oo(c) = -sin(grid % xc(c))*cos(grid % yc(c))
      V % n(c)  =  cos(grid % xc(c))*sin(grid % yc(c))
      V % o(c)  =  cos(grid % xc(c))*sin(grid % yc(c))
      V % oo(c) =  cos(grid % xc(c))*sin(grid % yc(c))
      W % n(c)  = 0.0
      W % o(c)  = 0.0
      W % oo(c) = 0.0
      P % n(c)  = 0.25*(cos(2*grid % xc(c)) + cos(2*grid % yc(c)))
    end do
  end if

  !---------------------------------!
  !      Calculate the inflow       !
  !   and initializes the Flux(s)   ! 
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
        Flux(s) = DEnc(material(c1))*( U % n(c2) * grid % sx(s) + &
                                       V % n(c2) * grid % sy(s) + &
                                       W % n(c2) * grid % sz(s) )
                                       
        if(TypeBC(c2)  ==  InFLOW) then
          if(material(c1) == m) bulk(m) % mass_in = bulk(m) % mass_in - Flux(s) 
          Stot  = sqrt(  grid % sx(s)**2  &
                       + grid % sy(s)**2  &
                       + grid % sz(s)**2)
          Area  = Area  + Stot
        endif
        if(TypeBC(c2)  ==  WALL)     n_wall        = n_wall        + 1 
        if(TypeBC(c2)  ==  INFLOW)   n_inflow      = n_inflow      + 1  
        if(TypeBC(c2)  ==  OUTFLOW)  n_outflow     = n_outflow     + 1 
        if(TypeBC(c2)  ==  SYMMETRY) n_symmetry    = n_symmetry    + 1 
        if(TypeBC(c2)  ==  WALLFL)   n_heated_wall = n_heated_wall + 1 
        if(TypeBC(c2)  ==  CONVECT)  n_convect     = n_convect     + 1 
      else
        Flux(s) = 0.0 
      end if
    end do
    call iglsum(n_wall)
    call iglsum(n_inflow)
    call iglsum(n_outflow)
    call iglsum(n_symmetry)
    call iglsum(n_heated_wall)
    call iglsum(n_convect)
    call glosum(bulk(m) % mass_in)
    call glosum(Area)
  end do                  

  !----------------------!
  !   Initializes time   ! 
  !----------------------!
  Time   = 0.0   
  if(this_proc  < 2) then
    write(*,*) '# MassIn=', bulk(1) % mass_in
    write(*,*) '# Average inflow velocity =', bulk(1) % mass_in / Area
    write(*,*) '# number of faces on the wall        : ', n_wall
    write(*,*) '# number of inflow faces             : ', n_inflow
    write(*,*) '# number of outflow faces            : ', n_outflow
    write(*,*) '# number of symetry faces            : ', n_symmetry
    write(*,*) '# number of faces on the heated wall : ', n_heated_wall
    write(*,*) '# number of convective outflow faces : ', n_convect
    write(*,*) '# Variables initialized !'
  end if

  end subroutine
