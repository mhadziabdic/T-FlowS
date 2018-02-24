!==============================================================================!
  subroutine Compute_Pressure_Simple(grid)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the S.I.M.p.L.E. method.            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Flow_Mod
  use Comm_Mod
  use Grid_Mod,     only: Grid_Type
  use Info_Mod
  use Numerics_Mod, only: errmax
  use Solvers_Mod,  only: Bicg, Cg, Cgs
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c, c1, c2, niter
  real              :: Us, Vs, Ws, A12, fs
  real              :: p_max, p_min
  real              :: error, tol
  real              :: smdpn
  real              :: dPxi, dPyi, dPzi
  character(len=80) :: precond
  real              :: urf           ! under-relaxation factor                 
!==============================================================================!
!     
!   The form of equations which I am solving:    
!     
!      /               /            
!     |               |             
!     | rho u dS = dt | GRAD pp dS
!     |               |             
!    /               /              
!
!   Dimension of the system under consideration
!   
!     [App] {pp} = {bpp}               [kg/s]
!   
!   Dimensions of certain variables
!
!     APP            [ms]
!     pp             [kg/ms^2]
!     b              [kg/s]
!     Flux           [kg/s]
!   
!------------------------------------------------------------------------------!

  ! Initialize matrix and right hand side
  A % val = 0.0
  b = 0.0 

  !-----------------------------------------!
  !   Initialize the pressure corrections   !
  !-----------------------------------------!
  pp % n = 0.0 

  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    ! Face is inside the domain
    if(c2  > 0 .or. c2  < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER) then

      smdpn = (  grid % sx(s)*grid % sx(s)   &
               + grid % sy(s)*grid % sy(s)   &
               + grid % sz(s)*grid % sz(s) ) &
            / (  grid % sx(s)*grid % dx(s)   & 
               + grid % sy(s)*grid % dy(s)   &
               + grid % sz(s)*grid % dz(s) )  

      ! Interpolate velocity 
      Us = fs * U % n(c1) + (1.0-fs) * U % n(c2)
      Vs = fs * V % n(c1) + (1.0-fs) * V % n(c2)
      Ws = fs * W % n(c1) + (1.0-fs) * W % n(c2)

      ! Calculate coeficients for the system matrix
      if(c2  > 0) then 
        A12 = 0.5 * density * smdpn *           &
           (  grid % vol(c1) / A % sav(c1)       &
            + grid % vol(c2) / A % sav(c2) )  
        A % val(A % pos(1,s))  = -A12
        A % val(A % pos(2,s))  = -A12
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
        A % val(A % dia(c2)) = A % val(A % dia(c2)) +  A12
      else 
        A12 = 0.5 * density * smdpn *           &
             (  grid % vol(c1) / A % sav(c1)     &
              + grid % vol(c2) / A % sav(c2) )  
        A % bou(c2)  = -A12
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
      end if 

      ! Interpolate pressure gradients
      dPxi=.5*( p % x(c1) + p % x(c2) )*grid % dx(s)
      dPyi=.5*( p % y(c1) + p % y(c2) )*grid % dy(s)
      dPzi=.5*( p % z(c1) + p % z(c2) )*grid % dz(s)

      ! Calculate flux through cell face
      flux(s) = density * (  Us*grid % sx(s)       &
                           + Vs*grid % sy(s)       &
                           + Ws*grid % sz(s) )     &
              + A12 * (p % n(c1) - p % n(c2))   &
              + A12 * (dPxi + dPyi + dPzi)                            

      b(c1)=b(c1)-flux(s)
      if(c2  > 0) b(c2)=b(c2)+flux(s)

    ! Side is on the boundary
    else ! (c2 < 0)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) == INFLOW) then 
        Us = U % n(c2)
        Vs = V % n(c2)
        Ws = W % n(c2)
        flux(s) = density * (  Us * grid % sx(s)  &
                             + Vs * grid % sy(s)  &
                             + Ws * grid % sz(s) )
        b(c1) = b(c1)-flux(s)
      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) == OUTFLOW .or.   &
              Grid_Mod_Bnd_Cond_Type(grid,c2) == CONVECT) then 
        Us = U % n(c2)
        Vs = V % n(c2)
        Ws = W % n(c2)
        flux(s) = density * (  Us*grid % sx(s)  &
                             + Vs*grid % sy(s)  &
                             + Ws*grid % sz(s) )
        b(c1) = b(c1)-flux(s)
        smdpn = (  grid % sx(s) * grid % sx(s)   &
                 + grid % sy(s) * grid % sy(s)   &
                 + grid % sz(s) * grid % sz(s) ) &
              / (  grid % sx(s) * grid % dx(s)   &
                 + grid % sy(s) * grid % dy(s)   &
                 + grid % sz(s) * grid % dz(s) )  
        A12 = density * smdpn * grid % vol(c1) / A % sav(c1)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) == PRESSURE) then
        Us = U % n(c1)
        Vs = V % n(c1)
        Ws = W % n(c1)
        flux(s) = density * (  Us * grid % sx(s)  &
                             + Vs * grid % sy(s)  &
                             + Ws * grid % sz(s) )
        b(c1) = b(c1)-flux(s)
        smdpn = ( grid % sx(s) * grid % sx(s)   &
                + grid % sy(s) * grid % sy(s)   &
                + grid % sz(s) * grid % sz(s) ) &
              / ( grid % sx(s) * grid % dx(s)   &
                + grid % sy(s) * grid % dy(s)   &
                + grid % sz(s) * grid % dz(s) )
        A12 = density * smdpn * grid % vol(c1) / A % sav(c1)
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
      else  ! it is SYMMETRY
        flux(s) = 0.0
      end if
    end if

1 end do

  errmax=0.0
  do c=1,grid % n_cells
    errmax=max(errmax, abs(b(c)))
  end do

  ! Don't solve the pressure corection too accurate.
  ! Value 1.e-18 blows the solution.
  ! Value 1.e-12 keeps the solution stable
  call Control_Mod_Tolerance_For_Pressure_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  niter=40
  call bicg(A, pp % n, b, precond, niter, tol, pp % res, error) 
  call Info_Mod_Iter_Fill_At(1, 3, pp % name, niter, pp % res)

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!
  call Control_Mod_Simple_Underrelaxation_For_Pressure(urf)  ! retreive urf
  do c = 1, grid % n_cells
    p % n(c)  =  p % n(c)  +  urf  *  pp % n(c)
  end do

  !------------------------------------!
  !   Normalize the pressure field     !
  !------------------------------------!
  ! p_max  = maxval(p % n(1:grid % n_cells))
  ! p_min  = minval(p % n(1:grid % n_cells))
 
  ! call Comm_Mod_Global_Max_Real(p_max)
  ! call Comm_Mod_Global_Min_Real(p_min)
 
  ! p % n = p % n - 0.5*(p_max+p_min)
 
  call Comm_Mod_Exchange(grid, pp % n)

  end subroutine
