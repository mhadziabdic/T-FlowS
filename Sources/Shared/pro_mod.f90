!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!        for the processor        !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
module pro_mod

  use Var_Mod
  use Matrix_Mod

  implicit none

  ! Right hand side for velocity and pressure equations 
  type(Matrix_Type) :: A  ! system matrix for all variables
  real,allocatable  :: b(:)

  real,allocatable :: VAR1x(:),   VAR1y(:),   VAR1z(:)
  real,allocatable :: VAR2x(:),   VAR2y(:),   VAR2z(:)
  real,allocatable :: VAR3x(:),   VAR3y(:),   VAR3z(:)
  real,allocatable :: VAR4x(:),   VAR4y(:),   VAR4z(:)
  real,allocatable :: VAR5x(:),   VAR5y(:),   VAR5z(:)
  real,allocatable :: VAR6x(:),   VAR6y(:),   VAR6z(:)
  real,allocatable :: VAR7x(:),   VAR7y(:),   VAR7z(:)
  real,allocatable :: VAR8x(:),   VAR8y(:),   VAR8z(:)
  real,allocatable :: VAR9x(:),   VAR9y(:),   VAR9z(:)
  real,allocatable :: VAR10x(:),  VAR10y(:),  VAR10z(:)
  real,allocatable :: VAR11x(:),  VAR11y(:),  VAR11z(:)
  real,allocatable :: VAR12x(:),  VAR12y(:),  VAR12z(:)
  real,allocatable :: PHI1x(:),   PHI1y(:),   PHI1z(:)
  real,allocatable :: PHI2x(:),   PHI2y(:),   PHI2z(:)
  real,allocatable :: PHI3x(:),   PHI3y(:),   PHI3z(:)
  real,allocatable :: PHI4x(:),   PHI4y(:),   PHI4z(:)
  real,allocatable :: PHI5x(:),   PHI5y(:),   PHI5z(:)
  real,allocatable :: PHI6x(:),   PHI6y(:),   PHI6z(:)
  real,allocatable :: PHI7x(:),   PHI7y(:),   PHI7z(:)
  real,allocatable :: PHI8x(:),   PHI8y(:),   PHI8z(:)
  real,allocatable :: PHI9x(:),   PHI9y(:),   PHI9z(:)
  real,allocatable :: PHI10x(:),  PHI10y(:),  PHI10z(:)
  real,allocatable :: PHI11x(:),  PHI11y(:),  PHI11z(:)

  real,allocatable :: PHIx(:),   PHIy(:),   PHIz(:)
  real,allocatable :: PHIside(:)

  ! For advection schemes
  real,allocatable :: phi_max(:), phi_min(:) 

  ! Total dissipation in HJ model
  real,allocatable :: eps_tot(:)

  ! Velocity components
  type(Var_Type) :: u
  type(Var_Type) :: v
  type(Var_Type) :: w

  ! Temperature
  type(Var_Type) :: t
  type(Var_Type) :: tt
  type(Var_Type) :: ut
  type(Var_Type) :: vt
  type(Var_Type) :: wt

  ! Pressure 
  type(Var_Type) :: p  
  type(Var_Type) :: pp

  ! Turbulent viscosity
  real,allocatable :: VISt_sgs(:)
  real,allocatable :: VISt_eff(:)
  real,allocatable :: Ptt(:)

  !-------------------------!
  !   Algorythm parameters  !
  !-------------------------!
  integer, parameter :: AB       =  1
  integer, parameter :: CN       =  2
  integer, parameter :: FI       =  3
  integer, parameter :: LIN      =  4
  integer, parameter :: PAR      =  5
  integer, parameter :: SIMPLE   =  6
  integer, parameter :: FRACT    =  7

  !-----------------------!
  !   Advection schemes   !
  !-----------------------!
  integer, parameter :: CDS       = 20
  integer, parameter :: QUICK     = 21 
  integer, parameter :: LUDS      = 22
  integer, parameter :: MINMOD    = 23
  integer, parameter :: SMART     = 24
  integer, parameter :: AVL_SMART = 25
  integer, parameter :: SUPERBEE  = 26
  integer, parameter :: GAMMA     = 27
  
  !-----------------------!
  !   Turbulence models   !
  !-----------------------!
  integer, parameter :: LES      = 12
  integer, parameter :: DNS      = 13
  integer, parameter :: K_EPS    = 14
  integer, parameter :: K_EPS_VV = 15
  integer, parameter :: SPA_ALL  = 16
  integer, parameter :: DES_SPA  = 17
  integer, parameter :: ZETA     = 18
  integer, parameter :: EBM      = 19
  integer, parameter :: LOW_RE   = 29    
  integer, parameter :: HIGH_RE  = 30
  integer, parameter :: J_L      = 31  
  integer, parameter :: S_L_Y    = 32 
  integer, parameter :: NAG      = 33
  integer, parameter :: RNG      = 34
  integer, parameter :: SMAG     = 35
  integer, parameter :: WALE     = 36
  integer, parameter :: DYN      = 37
  integer, parameter :: MIX      = 38
  integer, parameter :: HYB_ZETA = 39
  integer, parameter :: HYB_PITM = 40
  integer, parameter :: HJ       = 41
  integer, parameter :: ZETAM    = 42
  integer, parameter :: HYB      = 43
  integer, parameter :: PURE     = 44  ! as opposite of hybrid ;)

  ! Mass fluxes throught cell faces
  real,allocatable :: Flux(:) 

  ! Geometrical staff 
  real,allocatable :: Scoef(:)
  real,allocatable :: G(:,:) 
  real,allocatable :: fF(:)   ! weight factors for the fluid phase

  integer,allocatable :: CellFace(:,:)
  integer,allocatable :: WallFace(:)

  logical,allocatable :: IsNearWall(:)

  ! Cells which are bad for calculation of gradients
  logical,allocatable :: BadForG(:)
  integer,allocatable :: NumGood(:),   & 
                         NumNeig(:)

  ! Mass fluxes, bulk velocities and pressure drops for each material
  real,allocatable :: MassIn(:), MasOut(:) 
  real,allocatable :: FLUXx(:),  FLUXy(:),  FLUXz(:)
  real,allocatable :: FLUXoX(:), FLUXoY(:), FLUXoZ(:) 
  real,allocatable :: Ubulk(:),  Vbulk(:),  Wbulk(:)
  real,allocatable :: PdropX(:), PdropY(:), PdropZ(:)

  ! Viscosity, Density, Conductivity
  integer :: StateMat(100)
  integer :: SimulMat(100)
  real    :: VISc, DENc(100), CONc(100), CAPc(100)

  ! angular velocity 
  real    :: omegaX, omegaY, omegaZ, omega

  ! turbulent prandtl number 
  real    :: Prt, Numax


  ! Time step and total time
  real    :: dt, Time

  ! Maarten - Thickness of boundary layer
  real    :: Delta_m

  ! Constants needed for UserProbe2d (cut lines)
  real      :: x_o, y_o, Href
  integer   :: Ncuts
  character :: namCut*80

  ! Integer variable needed for interpolation of
  ! results between different meshes tranfer (LoaIni)
  integer             :: NClast, N_sign, eqn
  integer,allocatable :: near(:)

  ! Residuals                
  real :: errmax, res(100)  

  ! Monitoring planes for each material (domain)
  real,allocatable :: xp(:), yp(:), zp(:)

  !---------------------------!
  !     Solver parameters     !
  !---------------------------!
  real    :: URFC(100), SIMTol, URFC_Tur(100), URFC_Tem(100)
  real    :: TRFC(100)

  ! Under-relaxation parameter for turbulent quantity
  real    :: URFT, Alfa_fin1, Alfa_fin2

  !-----------------------------------!
  !     Area of the cross section     !
  !-----------------------------------!
  real,allocatable :: AreaX(:), AreaY(:), AreaZ(:)           
  real :: Area, Tflux, Qflux, Xmax, Ymax, Zmax, Tref, Tinf           

  !------------------------------!
  !     Algorythm parameters     !
  !------------------------------!
  integer :: INERT,    CONVEC,    CROSS,    DIFFUS 
  integer :: ALGOR
  integer :: SIMULA
  integer :: POSPRO
  integer :: CHANNEL,  TEST,      OTHER,    HOT, HOTini, PIPE, JET, ROT, TGV, URANS 
  integer :: BACKSTEP, AHILL, RB_CONV
  integer :: ROUGH, PER_BC, BUOY 
  integer :: SGDH, GGDH, BS
  integer :: SHAKE(100),    BLEND(100),BLEND_TUR(100), BLEND_TEM(100)
  integer :: SHAKE_PER(100),SHAKE_INT(100)
  integer :: PREC 
  integer :: RES_INI, MODE, BUDG 
  integer :: XHOM,     YHOM,     ZHOM

  integer,parameter :: MAXM=100 
  integer :: Cm(MAXM), Nmon
  real    :: NOM(MAXM), DEN(MAXM), R11(MAXM), U_f(MAXM)

  integer :: Ndt, Ndtt, Nstat, Nini, ini, Ndyn, Nstat2, NewSta, NK, Nbudg 

  !------------------------------------------------------------------!
  ! LineMon:   1:  6 -> Time step 
  ! ~~~~~~~~   7: 18 -> Time                
  !           19: 66 -> U,V,W,P monitoring
  !           67: 78 -> T  monitoring 
  !           79: 90 -> FLUXx
  !           91:102 -> P drop
  !          103:114 -> CFL
  !          115:126 -> Pe
  !          127:138 -> Kin.en.           
  !------------------------------------------------------------------!
  character*138 :: LinMon0 ! everything that goes on the screen
  character*138 :: LinMon1 ! everything that goes on the screen
  character*138 :: LinMon2 ! everything that goes on the screen
  character*138 :: LinMon3 ! everything that goes on the screen
  character*138 :: LinMon4 ! everything that goes on the screen
  character*138 :: LinMon5 ! everything that goes on the screen
  character*138 :: LinMon6 ! everything that goes on the screen
  character*138 :: LinMon7 ! everything that goes on the screen
  character*138 :: LinMon8 ! everything that goes on the screen
  character*138 :: LinMon9 ! everything that goes on the screen

  !------------------------------------------------------------------!
  ! LineRes:   1:  1 -> #
  ! ~~~~~~~~   2:  4 -> ini
  !            5: 16 -> errmax 
  !           17: 28 -> res U
  !           29: 40 -> res V
  !           41: 52 -> res W
  !           53: 64 -> res PP
  !           65: 76 -> res T
  !           77: 80 -> iter U
  !           81: 84 -> iter V
  !           85: 88 -> iter W
  !           89: 92 -> iter P
  !           93: 96 -> iter T
  !------------------------------------------------------------------!
  character*100 :: LineRes              ! everything that goes on the screen
  character :: namIni(128)*80

end module
