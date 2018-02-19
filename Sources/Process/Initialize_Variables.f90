!==============================================================================!
  subroutine Initialize_Variables(grid)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Flow_Mod
  use les_mod
  use par_mod
  use rans_mod
  use Grid_Mod
  use Bulk_Mod
  use Control_Mod
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

  area  = 0.0
  print *, 'grid % n_materials: ', grid % n_materials
  do n = 1, grid % n_materials
    do c = 1, grid % n_cells
      u % mean(c) = 0.0
      v % mean(c) = 0.0
      w % mean(c) = 0.0
      u % n(c)    = u % init(material(c))
      u % o(c)    = u % init(material(c)) 
      u % oo(c)   = u % init(material(c))
      v % n(c)    = v % init(material(c)) 
      v % o(c)    = v % init(material(c))
      v % oo(c)   = v % init(material(c))
      w % n(c)    = w % init(material(c)) 
      w % o(c)    = w % init(material(c))
      w % oo(c)   = w % init(material(c))
      if(heat_transfer == YES) then
        t % n(c)  = t % init(material(c)) 
        t % o(c)  = t % init(material(c)) 
        t % oo(c) = t % init(material(c)) 
        Tinf      = t % init(material(c))
      end if 
      if(turbulence_model == REYNOLDS_STRESS_MODEL .or.  &
         turbulence_model == HANJALIC_JAKIRLIC) then
        uu % n(c)  = uu % init(material(c))
        vv % n(c)  = vv % init(material(c))
        ww % n(c)  = ww % init(material(c))
        uv % n(c)  = uv % init(material(c))
        uw % n(c)  = uw % init(material(c))
        vw % n(c)  = vw % init(material(c))
        eps % n(c) = eps % init(material(c))
        if(turbulence_model == REYNOLDS_STRESS_MODEL) then
          f22 % n(c) = f22 % init(material(c))
        end if
      end if
      if(turbulence_model == K_EPS .or.  &
         turbulence_model == HYBRID_PITM) then
        kin % n(c)  = kin % init(material(c))
        kin % o(c)  = kin % init(material(c))
        kin % oo(c) = kin % init(material(c))
        eps % n(c)  = eps % init(material(c))
        eps % o(c)  = eps % init(material(c))
        eps % oo(c) = eps % init(material(c))
        Uf(c)       = 0.047
        Ynd(c)      = 30.0
      end if
      if(turbulence_model == K_EPS_V2      .or.  &
         turbulence_model == K_EPS_ZETA_F  .or.  & 
         turbulence_model == HYBRID_K_EPS_ZETA_F) then
        kin % n(c)  = kin % init(material(c))
        kin % o(c)  = kin % init(material(c))
        kin % oo(c) = kin % init(material(c))
        eps % n(c)  = eps % init(material(c))
        eps % o(c)  = eps % init(material(c))
        eps % oo(c) = eps % init(material(c))
        f22 % n(c)  = f22 % init(material(c))
        f22 % o(c)  = f22 % init(material(c))
        f22 % oo(c) = f22 % init(material(c))
        v2  % n(c)  = v2  % init(material(c))
        v2  % o(c)  = v2  % init(material(c))
        v2  % oo(c) = v2  % init(material(c))
        Uf(c)       = 0.047
        Ynd(c)      = 30.0
      end if
      if(turbulence_model == SPALART_ALLMARAS .or.  &
         turbulence_model == DES_SPALART) then      
        VIS % n(c)  = VIS % init(material(c))
        VIS % o(c)  = VIS % init(material(c))
        VIS % oo(c) = VIS % init(material(c))
      end if
    end do 
  end do   !end do n=1,grid % n_materials

  if(TGV == YES) then
    do c = 1, grid % n_cells
      u % n(c)  = -sin(grid % xc(c))*cos(grid % yc(c))
      u % o(c)  = -sin(grid % xc(c))*cos(grid % yc(c))
      u % oo(c) = -sin(grid % xc(c))*cos(grid % yc(c))
      v % n(c)  =  cos(grid % xc(c))*sin(grid % yc(c))
      v % o(c)  =  cos(grid % xc(c))*sin(grid % yc(c))
      v % oo(c) =  cos(grid % xc(c))*sin(grid % yc(c))
      w % n(c)  = 0.0
      w % o(c)  = 0.0
      w % oo(c) = 0.0
      P % n(c)  = 0.25*(cos(2*grid % xc(c)) + cos(2*grid % yc(c)))
    end do
  end if

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
                                       
        if(Grid_Mod_Bnd_Cond_Type(grid,c2)  ==  InFLOW) then
          if(material(c1) == m) bulk(m) % mass_in = bulk(m) % mass_in - flux(s) 
          Stot  = sqrt(  grid % sx(s)**2  &
                       + grid % sy(s)**2  &
                       + grid % sz(s)**2)
          area  = area  + Stot
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
    call iglsum(n_wall)
    call iglsum(n_inflow)
    call iglsum(n_outflow)
    call iglsum(n_symmetry)
    call iglsum(n_heated_wall)
    call iglsum(n_convect)
    call glosum(bulk(m) % mass_in)
    call glosum(area)
  end do                  

  !----------------------!
  !   Initializes time   ! 
  !----------------------!
  if(this_proc  < 2) then
    print *, '# MassIn=', bulk(1) % mass_in
    print *, '# Average inflow velocity =', bulk(1) % mass_in / area
    print *, '# number of faces on the wall        : ', n_wall
    print *, '# number of inflow faces             : ', n_inflow
    print *, '# number of outflow faces            : ', n_outflow
    print *, '# number of symetry faces            : ', n_symmetry
    print *, '# number of faces on the heated wall : ', n_heated_wall
    print *, '# number of convective outflow faces : ', n_convect
    print *, '# Variables initialized !'
  end if

  end subroutine
