!==============================================================================!
  module Control_Mod
!------------------------------------------------------------------------------!
  use Tokenizer_Mod
  implicit none
!==============================================================================!

  character(len=80)  :: control_file_name = 'control'
  integer, parameter :: CONTROL_FILE_UNIT = 10

  contains

  !--------------------!
  !   Input / Output   !
  !--------------------!

    ! Basic functionality
    include 'Control_Mod/Input_Output/Open_File.f90'
    include 'Control_Mod/Input_Output/Read_Char_Item.f90'
    include 'Control_Mod/Input_Output/Read_Int_Item.f90'
    include 'Control_Mod/Input_Output/Read_Real_Item.f90'
    include 'Control_Mod/Input_Output/Read_Real_Array.f90'
    include 'Control_Mod/Input_Output/Write_File.f90'

    ! Load
    include 'Control_Mod/Input_Output/Problem_Name.f90'
    include 'Control_Mod/Input_Output/Load_Restart_Name.f90'
    include 'Control_Mod/Input_Output/Save_Restart_Name.f90'
    include 'Control_Mod/Input_Output/Load_Initial_Solution_Name.f90'
    include 'Control_Mod/Input_Output/Save_Initial_Solution_Name.f90'

    ! Save
    include 'Control_Mod/Input_Output/Backup_Save_Interval.f90'
    include 'Control_Mod/Input_Output/Results_Save_Interval.f90'
    include 'Control_Mod/Input_Output/Monitoring_Points.f90'

  !-------------!
  !   Physics   !
  !-------------!

    ! Heat transfer
    include 'Control_Mod/Physics/Heat_Transfer.f90'
    include 'Control_Mod/Physics/Reference_Temperature.f90'

    ! Turbulence
    include 'Control_Mod/Physics/Turbulence_Model.f90'
    include 'Control_Mod/Physics/Turbulence_Model_Variant.f90'
    include 'Control_Mod/Physics/Roughness_Coefficient.f90'

    ! Other environmental conditions
    include 'Control_Mod/Physics/Angular_Velocity_Vector.f90'
    include 'Control_Mod/Physics/Gravitational_Vector.f90'
    include 'Control_Mod/Physics/Mass_Flow_Rates.f90'
    include 'Control_Mod/Physics/Pressure_Drops.f90'
    include 'Control_Mod/Physics/Point_For_Monitoring_Plane.f90'

    ! Multiphase 
    include 'Control_Mod/Physics/Number_Of_Phases.f90'

    ! Species    
    include 'Control_Mod/Physics/Number_Of_Species.f90'

    ! User scalars
    include 'Control_Mod/User/Number_Of_User_Scalars.f90'

  !--------------!
  !   Numerics   !
  !--------------!

    ! Time Stepping
    include 'Control_Mod/Numerics/Time_Step.f90'
    include 'Control_Mod/Numerics/Number_Of_Time_Steps.f90'
    include 'Control_Mod/Numerics/Starting_Time_Step_For_Statistics.f90'

    ! Discretization
    include 'Control_Mod/Numerics/Pressure_Momentum_Coupling.f90'
    include 'Control_Mod/Numerics/Advection_Scheme_For_Energy.f90'
    include 'Control_Mod/Numerics/Advection_Scheme_For_Momentum.f90'
    include 'Control_Mod/Numerics/Advection_Scheme_For_Turbulence.f90'
    include 'Control_Mod/Numerics/Blending_Coefficient_Energy.f90'
    include 'Control_Mod/Numerics/Blending_Coefficient_Momentum.f90'
    include 'Control_Mod/Numerics/Blending_Coefficient_Turbulence.f90'
    include 'Control_Mod/Numerics/Max_Simple_Iterations.f90'
    include 'Control_Mod/Numerics/Min_Simple_Iterations.f90'
    include 'Control_Mod/Numerics/Simple_Underrelaxation_For_Momentum.f90'
    include 'Control_Mod/Numerics/Simple_Underrelaxation_For_Pressure.f90'
    include 'Control_Mod/Numerics/Simple_Underrelaxation_For_Energy.f90'
    include 'Control_Mod/Numerics/Simple_Underrelaxation_For_Turbulence.f90'
    include 'Control_Mod/Numerics/Time_Integration_For_Inertia.f90'
    include 'Control_Mod/Numerics/Time_Integration_For_Advection.f90'
    include 'Control_Mod/Numerics/Time_Integration_For_Diffusion.f90'
    include 'Control_Mod/Numerics/Time_Integration_For_Cross_Diffusion.f90'

    ! Linear solvers
    include 'Control_Mod/Numerics/Solver_For_Momentum.f90'
    include 'Control_Mod/Numerics/Solver_For_Pressure.f90'
    include 'Control_Mod/Numerics/Solver_For_Energy.f90'
    include 'Control_Mod/Numerics/Solver_For_Turbulence.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Momentum_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Pressure_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Energy_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Turbulence_Solver.f90'
    include 'Control_Mod/Numerics/Tolerance_For_Simple_Algorithm.f90'
    include 'Control_Mod/Numerics/Preconditioner_For_System_Matrix.f90'

  end module