!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                 !                                    !
!                                 !   Bojan Niceno                     !
!   Global variable definitions   !   Delft University of Technology   !
!        for the processor        !   Section Heat Transfer            !
!                                 !   niceno@duttwta.wt.tn.tudelft.nl  !
!                                 !                                    !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
module pro_mod

  use allp_mod, only: YES, NO
  use Var_Mod
  use Matrix_Mod

  implicit none

  ! Right hand side for velocity and pressure equations 
  type(Matrix_Type) :: A  ! system matrix for all variables
  real,allocatable  :: b(:)

  ! Used in Dynamic Smgaorinsky model 
  real,allocatable    :: Aval_dif(:)

  ! Correlation points
  real :: R11_1, R11_2, R11_3, R11_4, R11_5
  real :: R11_6, R11_7, R11_8, R11_9, R11_10
  real :: A11_1, A11_2, A11_3, A11_4, A11_5
  real :: A11_6, A11_7, A11_8, A11_9, A11_10


  ! Velocity derivativeses dP/dx .... 
  real,allocatable :: Ux(:), Uy(:), Uz(:)
  real,allocatable :: Vx(:), Vy(:), Vz(:)
  real,allocatable :: Wx(:), Wy(:), Wz(:)

  ! Pressure derivativeses dP/dx .... 
  real,allocatable :: Px(:), Py(:), Pz(:)

  real,allocatable :: Kx(:)

  ! Pressure at the cell faces  
  real,allocatable :: Ps(:)

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

  ! Velocity components
  type(Var_Type) :: U
  type(Var_Type) :: V
  type(Var_Type) :: W

  ! Temperature
  type(Var_Type) :: T
  type(Var_Type) :: tt
  type(Var_Type) :: ut
  type(Var_Type) :: vt
  type(Var_Type) :: wt

  ! Pressure 
  type(Var_Type) :: P  
  type(Var_Type) :: PP

!=====================================================================!
!        Hybrid apriori variables
!=====================================================================!

  ! Velocity components
  type(Var_Type) :: U_r
  type(Var_Type) :: V_r
  type(Var_Type) :: W_r

  ! Pressure
  type(Var_Type) :: P_r
  type(Var_Type) :: PP_r

  ! Turbulent viscosity
  real,allocatable :: VISt_sgs(:)
  real,allocatable :: VISt_eff(:)
  real,allocatable :: Ptt(:)

  ! Reynolds stresses
  type(Var_Type) :: uu_r, vv_r, ww_r, uv_r, uw_r, vw_r

  ! Mass fluxes throught cell faces
  real,allocatable :: Flux_r(:), Alfa_lim(:)

  ! Mass fluxes throught the whole domain
  real,allocatable :: FLUXx_r(:),  FLUXy_r(:),  FLUXz_r(:)
!=====================================================================!
 
  !-------------------------!
  !   Algorythm parameters  !
  !-------------------------!
  integer :: K_EPS
  integer :: K_EPS_VV   
  integer :: HRe  
  integer :: MODE   
  integer :: LRe    
  integer :: SPA_ALL
  integer :: DES_SPA
  integer :: J_L    
  integer :: NAG     
  integer :: S_L_Y   
  integer :: WOLF   
  integer :: ZETA   
  integer :: HYB_ZETA   
  integer :: HYB_PITM   
  integer :: RNG   
  integer :: SMAG
  integer :: DYN 
  integer :: WALE
  integer :: MIX  
  integer :: ZPANS  
  integer :: ZETAM  
  integer :: EBM
  integer :: HYB
  integer :: HJ
  integer :: WF
  integer :: STAN
  integer :: BUOY

  ! Mass fluxes throught cell faces
  real,allocatable :: Flux(:) 

  ! Geometrical staff 
  real,allocatable :: Scoef(:)
  real,allocatable :: G(:,:) 
  real,allocatable :: fF(:)   ! weight factors for the fluid phase

  integer,allocatable :: CellFace(:,:)
  integer,allocatable :: WallFace(:)

  logical,allocatable :: IsNearWall(:)
  logical,allocatable :: IsNearPeri(:)
  logical,allocatable :: IsNearWall_2(:)
  logical,allocatable :: IsNearWall_3(:)
  logical,allocatable :: IsNearInflow(:)
  logical,allocatable :: ConvZone1(:)

  ! Cells which are bad for calculation of gradients
  logical,allocatable :: BadForG(:)
  integer,allocatable :: NumGood(:),   & 
                         NumNeig(:)

  ! Mass fluxes throught the whole domain
  real,allocatable :: MassIn(:), MasOut(:) 
  real,allocatable :: FLUXx(:),  FLUXy(:),  FLUXz(:)
  real,allocatable :: FLUXoX(:), FLUXoY(:), FLUXoZ(:) 
  real,allocatable :: Ubulk(:),  Vbulk(:),  Wbulk(:)

  ! Viscosity, Density, Conductivity
  integer :: StateMat(100)
  integer :: SimulMat(100)
  real    :: VISc, DENc(100), CONc(100), CAPc(100)

  ! Average velocity 
  real    :: Uaver

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
  integer          :: NClast, N_sign, eqn
  integer,allocatable :: near(:)
  integer,allocatable :: near_2(:)
  integer,allocatable :: near_3(:)
  integer,allocatable :: connect(:)
  integer,allocatable :: connect2(:)

  ! Residuals                
  real    :: errmax, res(100)  

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
  integer :: LIN,      PAR,       AB,       CN,       FI
  integer :: ALGOR,    SIMPLE,    FRACT
  integer :: SIMULA,   DNS,       LES 
  integer :: POSPRO,   AVS,       GMV
  integer :: CHANNEL,  TEST,      OTHER,    HOT, HOTini, PIPE, JET, ROT, TGV, BUDG, URANS 
  integer :: BACKSTEP, AHILL, RB_CONV
  integer :: ROUGH, PER_BC 
  integer :: SGDH, GGDH, BS
  integer :: SHAKE(100),    BLEND(100),BLEND_TUR(100), BLEND_TEM(100)
  integer :: SHAKE_PER(100),SHAKE_INT(100)
  integer :: PREC 
  integer :: CDS,      QUICK,    LUDS,     MINMOD,   SMART,    AVL_SMART, &
             SUPERBEE, GAMMA, RES_INI 
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
