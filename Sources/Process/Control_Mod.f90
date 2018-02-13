!==============================================================================!
  module Control_Mod
!------------------------------------------------------------------------------!
  use Tokenizer_Mod
  implicit none
!==============================================================================!

  character(len=80)  :: control_file_name = 'control'
  integer, parameter :: CONTROL_FILE_UNIT = 10

  contains

  !-------------------------!
  !   Basic functionality   !
  !-------------------------!
  include 'Control_Mod/Open_File.f90'
  include 'Control_Mod/Read_Char_Item.f90'
  include 'Control_Mod/Read_Int_Item.f90'
  include 'Control_Mod/Read_Real_Item.f90'
  include 'Control_Mod/Read_Real_Array.f90'
  include 'Control_Mod/Write_File.f90'

  !--------------------!
  !   Input / Output   !
  !--------------------!

    ! Load
    include 'Control_Mod/Problem_Name.f90'

    ! Save
    include 'Control_Mod/Backup_Save_Interval.f90'
    include 'Control_Mod/Results_Save_Interval.f90'

  !-------------!
  !   Physics   !
  !-------------!

    ! Heat transfer
    include 'Control_Mod/Reference_Temperature.f90'

    ! Turbulence
    include 'Control_Mod/Turbulence_Model.f90'
    include 'Control_Mod/Roughness_Coefficient.f90'

    ! Multiphase 
    include 'Control_Mod/Number_Of_Phases.f90'

    ! Species    
    include 'Control_Mod/Number_Of_Species.f90'

    ! User scalars
    include 'Control_Mod/Number_Of_User_Scalars.f90'

  !--------------!
  !   Numerics   !
  !--------------!

    ! Time Stepping
    include 'Control_Mod/Number_Of_Time_Steps.f90'
    include 'Control_Mod/Time_Step.f90'

    ! Discretization
    include 'Control_Mod/Max_Simple_Iterations.f90'
    include 'Control_Mod/Min_Simple_Iterations.f90'
    include 'Control_Mod/Simple_Underrelaxation_Momentum.f90'
    include 'Control_Mod/Simple_Underrelaxation_Pressure.f90'
    include 'Control_Mod/Simple_Underrelaxation_Energy.f90'
    include 'Control_Mod/Simple_Underrelaxation_Turbulence.f90'

    ! Linear solvers
    include 'Control_Mod/Linear_Solver_For_Momentum.f90'
    include 'Control_Mod/Linear_Solver_For_Pressure.f90'
    include 'Control_Mod/Linear_Solver_For_Energy.f90'
    include 'Control_Mod/Linear_Solver_For_Turbulence.f90'

    include 'Control_Mod/Tolerance_Momentum_Solver.f90'
    include 'Control_Mod/Tolerance_Pressure_Solver.f90'
    include 'Control_Mod/Tolerance_Energy_Solver.f90'
    include 'Control_Mod/Tolerance_Turbulence_Solver.f90'

  end module
