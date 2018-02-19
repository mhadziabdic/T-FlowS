!==============================================================================!
  module Flow_Mod
!------------------------------------------------------------------------------!
!   Module for basic flow field plus temperature.                              !
!   It is a bit of a mumbo-jumbo at this moment, it will furhter have to       !
!   differentiate into numerical and physica parts.                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod
  use Bulk_Mod
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

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
  real,allocatable :: flux(:) 

  ! Variables determining if we are dealing with heat transfer and buoyancy
  ! (both can be either YES or NO)
  integer :: heat_transfer
  integer :: buoyancy

  ! Geometrical staff 
  real,allocatable :: f_coef(:)  ! face coefficient
  real,allocatable :: g(:,:)     ! gradient matrices for each cell
  real,allocatable :: fw(:)      ! weight factors for the fluid phase

  ! Right hand side for velocity and pressure equations 
  type(Matrix_Type) :: A  ! system matrix for all variables
  real, allocatable :: b(:)

  real,allocatable :: phi_face(:)

  ! For advection schemes
  real, allocatable :: phi_max(:), phi_min(:) 

  ! Cells which are bad for calculation of gradients
  logical, allocatable :: BadForG(:)
  integer, allocatable :: NumGood(:),   & 
                          NumNeig(:)

  ! Mass fluxes, bulk velocities and pressure drops for each material
  type(Bulk_Type) :: bulk(100)

  ! Viscosity, Density, Conductivity
  integer :: StateMat(100)
  integer :: SimulMat(100)
  real    :: viscosity, density, conductivity, capacity

  ! angular velocity 
  real :: omega_x, omega_y, omega_z, omega

  ! Constants needed for UserProbe2d (cut lines)
  real      :: x_o, y_o, Href
  integer   :: Ncuts
  character :: namCut*80

  ! Integer variable needed for interpolation of
  ! results between different meshes tranfer (LoaIni)
  integer             :: NClast, N_sign, eqn

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
  integer :: ROUGH, PER_BC
  integer :: SGDH, GGDH, BS
  integer :: SHAKE(100)
  integer :: SHAKE_PER(100),SHAKE_INT(100)
  integer :: RES_INI, BUDG 
  integer :: XHOM,     YHOM,     ZHOM

  integer,parameter :: MAXM=100 
  integer :: Cm(MAXM), Nmon
  real    :: NOM(MAXM), DEN(MAXM), R11(MAXM), U_f(MAXM)

  end module
