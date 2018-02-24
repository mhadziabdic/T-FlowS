!==============================================================================!
  subroutine Compute_Pressure_Fractional(grid)
!------------------------------------------------------------------------------!
!   Forms and solves pressure equation for the fractional step method.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use Flow_Mod
  use Grid_Mod,     only: Grid_Type  
  use Comm_Mod
  use Info_Mod
  use Solvers_Mod,  only: Bicg, Cg, Cgs
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c, c1, c2, niter
  real              :: p_max, p_min
  real              :: ini_res, tol, mass_err
  real              :: Us, Vs, Ws, fs
  real              :: A12
  character(len=80) :: coupling
  character(len=80) :: precond
  real              :: urf           ! under-relaxation factor                 
!==============================================================================!
!     
!   The form of equations which are being solved:
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
!     pp,            [kg/ms^2]
!     b              [kg/s]
!     Flux           [kg/s]
!   
!------------------------------------------------------------------------------!

  ! Initialize matrix and source term
  A % val = 0.0
  b = 0.0 

  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s=1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    ! Face is inside the domain
    if( c2  > 0 .or. c2  < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) == BUFFER) then 

      ! Extract the "centred" pressure terms from cell velocities
      Us = fs*      (U % n(c1) + p % x(c1)*grid % vol(c1)/A % sav(c1))       &
         + (1.0-fs)*(U % n(c2) + p % x(c2)*grid % vol(c2)/A % sav(c2))

      Vs = fs*      (V % n(c1) + p % y(c1)*grid % vol(c1)/A % sav(c1))       &
         + (1.0-fs)*(V % n(c2) + p % y(c2)*grid % vol(c2)/A % sav(c2))

      Ws = fs*      (W % n(c1) + p % z(c1)*grid % vol(c1)/A % sav(c1))       &
         + (1.0-fs)*(W % n(c2) + p % z(c2)*grid % vol(c2)/A % sav(c2))

      ! Add the "staggered" pressure terms to face velocities
      Us= Us + (p % n(c1)-p % n(c2))*grid % sx(s)*                         &
         ( fs/A % sav(c1) + (1.0-fs)/A % sav(c2) )
      Vs=Vs+(p % n(c1)-p % n(c2))*grid % sy(s)*                            &
         ( fs/A % sav(c1) + (1.0-fs)/A % sav(c2) )
      Ws=Ws+(p % n(c1)-p % n(c2))*grid % sz(s)*                            &
         ( fs/A % sav(c1) + (1.0-fs)/A % sav(c2) )

      ! Now calculate the flux through cell face
      flux(s) = density * ( Us * grid % sx(s) +  &
                            Vs * grid % sy(s) +  &
                            Ws * grid % sz(s) )

      A12=density*(  grid % sx(s)*grid % sx(s)  &
                + grid % sy(s)*grid % sy(s)  &
                + grid % sz(s)*grid % sz(s))
      A12=A12*(fs/A % sav(c1)+(1.-fs)/A % sav(c2))

      if(c2  > 0) then 
        A % val(A % pos(1,s)) = -A12
        A % val(A % pos(2,s)) = -A12
        A % val(A % dia(c1))  = A % val(A % dia(c1)) +  A12
        A % val(A % dia(c2))  = A % val(A % dia(c2)) +  A12
      else
        A % bou(c2) = -A12
        A % val(A % dia(c1)) = A % val(A % dia(c1)) +  A12
      endif

      b(c1)=b(c1)-flux(s)
      if(c2  > 0) b(c2)=b(c2)+flux(s)

    ! Face is on the boundary
    else
      Us = U % n(c2)
      Vs = V % n(c2)
      Ws = W % n(c2)

      flux(s) = density * (  Us * grid % sx(s)  &
                           + Vs * grid % sy(s)  &
                           + Ws * grid % sz(s) )

      b(c1) = b(c1)-flux(s)
    end if

  end do

  !----------------------------------------!
  !   Initialize the pressure correction   !
  !----------------------------------------!
  pp % n = 0.0 

  mass_err = 0.0
  do c = 1, grid % n_cells
    mass_err = mass_err + abs(b(c))
  end do
  call Comm_Mod_Global_Sum_Real(mass_err)                       

  !--------------------------------------------!
  !   Solve the pressure correction equation   !
  !--------------------------------------------!

  ! Don't solve the pressure corection too accurate.
  ! Value 1.e-18 blows the solution.
  ! Value 1.e-12 keeps the solution stable
  call Control_Mod_Tolerance_For_Pressure_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  call Control_Mod_Pressure_Momentum_Coupling(coupling)
  if(coupling == 'PROJECTION') niter = 200
  if(coupling == 'SIMPLE')     niter =  15

  call Cg(A, pp % n, b, precond, niter, tol, ini_res, pp % res)

  call Info_Mod_Iter_Fill_At(1, 3, pp % name, niter, pp % res)   

  !-------------------------------!
  !   Update the pressure field   !
  !-------------------------------!
  urf = 1.0
  if(coupling == 'SIMPLE')  &
    call Control_Mod_Simple_Underrelaxation_For_Pressure(urf)

  p % n  =  p % n  +  urf  *  pp % n

  !----------------------------------!
  !   Normalize the pressure field   !
  !----------------------------------!
  p_max  = maxval(p % n(1:grid % n_cells))
  p_min  = minval(p % n(1:grid % n_cells))

  call Comm_Mod_Global_Max_Real(p_max) 
  call Comm_Mod_Global_Min_Real(p_min) 

  p % n  =  p % n  -  0.5 * (p_max + p_min)

  call Comm_Mod_Exchange(grid, pp % n) 

  end subroutine
