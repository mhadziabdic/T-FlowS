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
  use Bulk_Mod
  use Matrix_Mod

  implicit none

  ! Right hand side for velocity and pressure equations 
  type(Matrix_Type) :: A  ! system matrix for all variables
  real,allocatable  :: b(:)

  real,allocatable :: phi_face(:)

  ! For advection schemes
  real,allocatable :: phi_max(:), phi_min(:) 

  ! Velocity components
  type(Var_Type) :: u
  type(Var_Type) :: v
  type(Var_Type) :: w

  ! Temperature
  type(Var_Type) :: t

  ! Pressure 
  type(Var_Type) :: p  
  type(Var_Type) :: pp

  ! Turbulent viscosity
  real,allocatable :: Ptt(:)

  ! Mass fluxes throught cell faces
  real,allocatable :: Flux(:) 

  ! Geometrical staff 
  real,allocatable :: Scoef(:)
  real,allocatable :: G(:,:) 
  real,allocatable :: fF(:)   ! weight factors for the fluid phase

  logical,allocatable :: IsNearWall(:)

  ! Cells which are bad for calculation of gradients
  logical,allocatable :: BadForG(:)
  integer,allocatable :: NumGood(:),   & 
                         NumNeig(:)

  ! Mass fluxes, bulk velocities and pressure drops for each material
  type(Bulk_Type) :: bulk(100)

  ! Viscosity, Density, Conductivity
  integer :: StateMat(100)
  integer :: SimulMat(100)
  real    :: VISc, DENc(100), CONc(100), CAPc(100)

  ! angular velocity 
  real :: omega_x, omega_y, omega_z, omega

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

  !-----------------------------------!
  !     Area of the cross section     !
  !-----------------------------------!
  real :: Area, Tflux, Qflux, Xmax, Ymax, Zmax, Tref, Tinf           

  !------------------------------!
  !     Algorythm parameters     !
  !------------------------------!
  integer :: POSPRO
  integer :: CHANNEL,  TEST,      OTHER,  HOTini, PIPE, JET, ROT, TGV
  integer :: BACKSTEP, AHILL, RB_CONV
  integer :: ROUGH, PER_BC, BUOY 
  integer :: SGDH, GGDH, BS
  integer :: SHAKE(100)
  integer :: SHAKE_PER(100),SHAKE_INT(100)
  integer :: RES_INI, BUDG 
  integer :: XHOM,     YHOM,     ZHOM

  integer,parameter :: MAXM=100 
  integer :: Cm(MAXM), Nmon
  real    :: NOM(MAXM), DEN(MAXM), R11(MAXM), U_f(MAXM)

  end module
