!==============================================================================!
  module Parameters_Mod
!------------------------------------------------------------------------------!
!   Module containing all parameters (constants) for the processor.            !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-------------------------!
  !   Algorythm parameters  !
  !-------------------------!

  ! Time integration
  integer, parameter :: AB  = 1
  integer, parameter :: CN  = 2
  integer, parameter :: FI  = 3
  integer, parameter :: LIN = 4
  integer, parameter :: PAR = 5

  ! Velocity=pressure linging
  integer, parameter :: SIMPLE = 6
  integer, parameter :: FRACT  = 7

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

  end module
