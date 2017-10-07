!==============================================================================!
  program Create_Control_Directory
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: case_name
  character(len=80) :: directory
!------------------------------------------------------------------------------!
! Creates structure of the "Control_Directory.                                 ! 
!==============================================================================!

  case_name = "great_case"

  !-------------------!
  ! Control_Directory !
  !-------------------! 
  call execute_command_line("rm -fr Control_Directory")
  call execute_command_line("mkdir  Control_Directory")

    !------!
    ! Load !
    !------! 
    directory = "Control_Directory/Load/"
    call execute_command_line("mkdir "//trim(directory))

    open (9, file=trim(directory)//'CASE_NAME')
    write(9,*) "#========================="
    write(9,*) "# Specify case name here: "
    write(9,*) "#-------------------------"
    write(9,*) case_name
    write(9,*) "#-------------------------"
    close(9)

    open (9, file=trim(directory)//'INITIAL_SOLUTION_CASE_NAME')
    write(9,*) "#==========================================="
    write(9,*) "# Specify the name of the case with initial "
    write(9,*) "# solution or just type ""skip"" to cancel. "
    write(9,*) "#-------------------------------------------"
    write(9,*) "  SKIP"         
    write(9,*) "#-------------------------------------------"
    close(9)

    open (9, file=trim(directory)//'INTERPOLATED_INITIAL_SOLUTION_CASE_NAME')
    write(9,*) "#===================================================="
    write(9,*) "# Specify the name of the case with initial solution "
    write(9,*) "# from a different grid or type ""skip"" to cancel.  "
    write(9,*) "#----------------------------------------------------"
    write(9,*) "  SKIP"         
    write(9,*) "#----------------------------------------------------"
    close(9)

    !----------!
    ! Numerics !
    !----------! 
    directory = "Control_Directory/Numerics/"
    call execute_command_line("mkdir "//trim(directory))

      !----------------!
      ! Discretization !
      !----------------! 
      directory = "Control_Directory/Numerics/Discretization/"
      call execute_command_line("mkdir "//trim(directory))

      open (9, file=trim(directory)//'PRESSURE_VELOCITY_COUPLING')
      write(9,*) "#================================="
      write(9,*) "# Uncomment one of the following: "
      write(9,*) "#---------------------------------"
      write(9,*) "  SIMPLE"
      write(9,*) "# FRACTIONAL_STEP"
      write(9,*) "#---------------------------------"
      close(9)

      open (9, file=trim(directory)//'ADVECTION_SCHEME_UVW')
      write(9,*) "#================================="
      write(9,*) "# Uncomment one of the following: "
      write(9,*) "#---------------------------------"
      write(9,*) "# UDS" 
      write(9,*) "# CDS" 
      write(9,*) "  MINMOD" 
      write(9,*) "# SMART" 
      write(9,*) "# QUICK" 
      write(9,*) "# LUDS" 
      write(9,*) "# AVL_SMART" 
      write(9,*) "# SUPERBEE" 
      write(9,*) "# BLEND_CDS_UDS"
      write(9,*) "#---------------------------------"
      close(9)

      open (9, file=trim(directory)//'ADVECTION_SCHEME_T')
      write(9,*) "#================================="
      write(9,*) "# Uncomment one of the following: "
      write(9,*) "#---------------------------------"
      write(9,*) "  UDS" 
      write(9,*) "# CDS" 
      write(9,*) "# MINMOD" 
      write(9,*) "# SMART" 
      write(9,*) "# QUICK" 
      write(9,*) "# LUDS" 
      write(9,*) "# AVL_SMART" 
      write(9,*) "# SUPERBEE" 
      write(9,*) "# BLEND_CDS_UDS"
      write(9,*) "#---------------------------------"
      close(9)

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_KIN" )

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_EPS" )

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_V_2" )

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_F22" )

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_VIS" )

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_UU" )

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_VV" )

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_WW" )

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_UV" )

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_UW" )

      call execute_command_line("cp "//trim(directory)//        &
                                       "ADVECTION_SCHEME_T "//  &
                                       trim(directory)//        &
                                       "ADVECTION_SCHEME_VW" )

      open (9, file=trim(directory)//'SIMPLE_MAX_ITERATIONS')
      write(9,*) "#=============================================="
      write(9,*) "# Specify maximum number of SIMPLE iterations: "
      write(9,*) "#----------------------------------------------"
      write(9,*) "  16"
      write(9,*) "#----------------------------------------------"
      close(9)

      open (9, file=trim(directory)//'SIMPLE_MIN_ITERATIONS')
      write(9,*) "#=============================================="
      write(9,*) "# Specify maximum number of SIMPLE iterations: "
      write(9,*) "#----------------------------------------------"
      write(9,*) "   4"
      write(9,*) "#----------------------------------------------"
      close(9)

      open (9, file=trim(directory)//'SIMPLE_UNDERRELAXATION_UVW')
      write(9,*) "#============================================="
      write(9,*) "# Specify underrelaxation factor (0.0 - 1.0): "
      write(9,*) "#---------------------------------------------"
      write(9,*) "  0.7"
      write(9,*) "#---------------------------------------------"
      close(9)

      open (9, file=trim(directory)//'SIMPLE_UNDERRELAXATION_P')
      write(9,*) "#============================================="
      write(9,*) "# Specify underrelaxation factor (0.0 - 1.0): "
      write(9,*) "#---------------------------------------------"
      write(9,*) "  0.3"
      write(9,*) "#---------------------------------------------"
      close(9)

      open (9, file=trim(directory)//'SIMPLE_UNDERRELAXATION_T')
      write(9,*) "#============================================="
      write(9,*) "# Specify underrelaxation factor (0.0 - 1.0): "
      write(9,*) "#---------------------------------------------"
      write(9,*) "  0.6"
      write(9,*) "#---------------------------------------------"
      close(9)

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_KIN" )

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_EPS" )

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_V_2" )

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_F22" )

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_VIS" )

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_UU" )

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_VV" )

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_WW" )

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_UV" )

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_UW" )

      call execute_command_line("cp "//trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_T "//  &
                                       trim(directory)//              &
                                       "SIMPLE_UNDERRELAXATION_VW" )

      !----------------!
      ! Linear_Solvers !
      !----------------! 
      directory = "Control_Directory/Numerics/Linear_Solvers/"
      call execute_command_line("mkdir "//trim(directory))

      open (9, file=trim(directory)//'SOLVER_FOR_UVW')
      write(9,*) "#================================="
      write(9,*) "# Uncomment one of the following: "
      write(9,*) "#---------------------------------"
      write(9,*) "  BiCG"
      write(9,*) "# CG"  
      write(9,*) "# CGS"    
      write(9,*) "#---------------------------------"
      close(9)

      open (9, file=trim(directory)//'SOLVER_FOR_P')
      write(9,*) "#================================="
      write(9,*) "# Uncomment one of the following: "
      write(9,*) "#---------------------------------"
      write(9,*) "# BiCG"
      write(9,*) "# CG"  
      write(9,*) "  CGS"    
      write(9,*) "#---------------------------------"
      close(9)

      open (9, file=trim(directory)//'SOLVER_FOR_T')
      write(9,*) "#================================="
      write(9,*) "# Uncomment one of the following: "
      write(9,*) "#---------------------------------"
      write(9,*) "# BiCG"
      write(9,*) "  CG"  
      write(9,*) "# CGS"    
      write(9,*) "#---------------------------------"
      close(9)

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_KIN" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_EPS" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_V_2" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_F22" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_VIS" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_UU" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_VV" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_WW" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_UV" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_UW" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "SOLVER_FOR_T "//  &
                                       trim(directory)//  &
                                       "SOLVER_FOR_VW" )

      open (9, file=trim(directory)//'TOLERANCE_FOR_UVW')
      write(9,*) "#========================================"
      write(9,*) "# Specify solver tolerance for velocity: "
      write(9,*) "#----------------------------------------"
      write(9,*) "  1.0e-4"
      write(9,*) "#----------------------------------------"
      close(9)

      open (9, file=trim(directory)//'TOLERANCE_FOR_P')
      write(9,*) "#========================================"
      write(9,*) "# Specify solver tolerance for pressure: "
      write(9,*) "#----------------------------------------"
      write(9,*) "  1.0e-6" 
      write(9,*) "#----------------------------------------"
      close(9)

      call execute_command_line("cp "//trim(directory)//       &
                                       "TOLERANCE_FOR_UVW "//  &
                                       trim(directory)//       &
                                       "TOLERANCE_FOR_T" )

      call execute_command_line("cp "//trim(directory)//       &
                                       "TOLERANCE_FOR_UVW "//  &
                                       trim(directory)//       &
                                       "TOLERANCE_FOR_KIN" )

      call execute_command_line("cp "//trim(directory)//       &
                                       "TOLERANCE_FOR_UVW "//  &
                                       trim(directory)//       &
                                       "TOLERANCE_FOR_EPS" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "TOLERANCE_FOR_T "//  &
                                       trim(directory)//  &
                                       "TOLERANCE_FOR_V_2" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "TOLERANCE_FOR_T "//  &
                                       trim(directory)//  &
                                       "TOLERANCE_FOR_F22" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "TOLERANCE_FOR_T "//  &
                                       trim(directory)//  &
                                       "TOLERANCE_FOR_VIS" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "TOLERANCE_FOR_T "//  &
                                       trim(directory)//  &
                                       "TOLERANCE_FOR_UU" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "TOLERANCE_FOR_T "//  &
                                       trim(directory)//  &
                                       "TOLERANCE_FOR_VV" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "TOLERANCE_FOR_T "//  &
                                       trim(directory)//  &
                                       "TOLERANCE_FOR_WW" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "TOLERANCE_FOR_T "//  &
                                       trim(directory)//  &
                                       "TOLERANCE_FOR_UV" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "TOLERANCE_FOR_T "//  &
                                       trim(directory)//  &
                                       "TOLERANCE_FOR_UW" )

      call execute_command_line("cp "//trim(directory)//  &
                                       "TOLERANCE_FOR_T "//  &
                                       trim(directory)//  &
                                       "TOLERANCE_FOR_VW" )

      !---------------!
      ! Time_Stepping !
      !---------------! 
      directory = "Control_Directory/Numerics/Time_Stepping/"
      call execute_command_line("mkdir "//trim(directory))

      open (9, file=trim(directory)//'TIME_STEP')
      write(9,*) "#======================================="
      write(9,*) "# Specify time step (usually < 1.0e-3): "
      write(9,*) "#---------------------------------------"
      write(9,*) "  0.0005"
      write(9,*) "#---------------------------------------"
      close(9)

      open (9, file=trim(directory)//'NUMBER_OF_TIME_STEPS')
      write(9,*) "#==============================="
      write(9,*) "# Specify number of time steps: "
      write(9,*) "#-------------------------------"
      write(9,*) "  100"
      write(9,*) "#-------------------------------"
      close(9)

    !---------!
    ! Physics !
    !---------! 
    directory = "Control_Directory/Physics/"
    call execute_command_line("mkdir "//trim(directory))

      !------------!
      ! Multiphase !
      !------------! 
      directory = "Control_Directory/Physics/Multiphase/"
      call execute_command_line("mkdir "//trim(directory))

      open (9, file=trim(directory)//'NUMBER_OF_PHASES')
      write(9,*) "#==========================="
      write(9,*) "# Specify number of phases: "
      write(9,*) "# (Only 1 is supported now) "
      write(9,*) "#---------------------------"
      write(9,*) "  1"
      write(9,*) "#---------------------------"
      close(9)

      open (9, file=trim(directory)//'PHASE_NAMES')
      write(9,*) "#================================================"
      write(9,*) "# For each phase, a unique name should be given: "
      write(9,*) "#------------------------------------------------"
      write(9,*) "  WATER"
      write(9,*) "  STEAM"
      write(9,*) "  DIESEL"
      write(9,*) "  DUST"  
      write(9,*) "#------------------------------------------------"
      close(9)

      !---------!
      ! Species !
      !---------! 
      directory = "Control_Directory/Physics/Species/"
      call execute_command_line("mkdir "//trim(directory))

      open (9, file=trim(directory)//'NUMBER_OF_SPECIES')
      write(9,*) "#========================================"
      write(9,*) "# Specify number of transported species: "
      write(9,*) "# (None supported at the moment) "
      write(9,*) "#----------------------------------------"
      write(9,*) "  0"
      write(9,*) "#----------------------------------------"
      close(9)

      open (9, file=trim(directory)//'SPECIES_NAMES')
      write(9,*) "#=================================================="
      write(9,*) "# For each species, a unique name should be given: "
      write(9,*) "#--------------------------------------------------"
      write(9,*) "  SALT"
      write(9,*) "  SUGAR"
      write(9,*) "#--------------------------------------------------"
      close(9)

      !------------!
      ! Turbulence !
      !------------! 
      directory = "Control_Directory/Physics/Turbulence/"
      call execute_command_line("mkdir "//trim(directory))

      open (9, file=trim(directory)//'TURBULENCE_MODEL')
      write(9,*) "#============================================================"
      write(9,*) "# Uncomment one of the following lines for turbulence model: "
      write(9,*) "#------------------------------------------------------------"
      write(9,*) "# K_EPS_LRE -> k-eps model with wall integration."
      write(9,*) "  K_EPS_HRE -> k-eps model with wall function."
      write(9,*) "# ZETA      -> zeta-f model with compound wall treatment."
      write(9,*) "# LES_DYN   -> LES with dynamics Smagorinsky SGS model."
      write(9,*) "# LES_SMAG  -> LES with Smagorinsky SGS model."
      write(9,*) "# LES_WALE  -> WALE SGS model."
      write(9,*) "# LES_MIX   -> Mixed SGS model." 
      write(9,*) "# DNS       -> No turbulence model."
      write(9,*) "# EBM       -> RSM with integration up to the wall."  
      write(9,*) "# EBM_HYB   -> RSM to the wall with computed C_mu."
      write(9,*) "# HJ        -> Hanjalic/Jakirlic?"
      write(9,*) "# HJ_HYB    -> Hanjalic/Jakirlic with computed C_mu."
      write(9,*) "# HYB_ZETA  -> Hybrid RANS/LES model with zeta-f model."  
      write(9,*) "# HYB_PITM  -> Hybrid RANS/LES model with k-eps model."  
      write(9,*) "# K_EPS_VV  -> Durbin original model."  
      write(9,*) "# DES_SPA   -> DES model."  
      write(9,*) "# SPA_ALL   -> Spalart-Allmaras model." 
      write(9,*) "#------------------------------------------------------------"
      close(9)

      open (9, file=trim(directory)//'TURBULENCE_MODEL_CONSTANTS_LES_SMAG')
      write(9,*) "#====================================================="
      write(9,*) "# Model constants for LES with Smagorinsky SGS model: "
      write(9,*) "#-----------------------------------------------------"
      write(9,*) "  0.1  -> Cs"                                        
      write(9,*) "#-----------------------------------------------------"
      close(9)

      !------------!
      ! Combustion !
      !------------! 
      directory = "Control_Directory/Physics/Combustion/"
      call execute_command_line("mkdir "//trim(directory))

    !------!
    ! Save !
    !------! 
    directory = "Control_Directory/Save/"
    call execute_command_line("mkdir "//trim(directory))

    open (9, file=trim(directory)//'BACKUP_FILE_NAME')
    write(9,*) "#================================"
    write(9,*) "# Specify backup file name here: "
    write(9,*) "#--------------------------------"
    write(9,*) case_name
    write(9,*) "#--------------------------------"
    close(9)

    open (9, file=trim(directory)//'BACKUP_SAVING_INTERVAL')
    write(9,*) "#================================="
    write(9,*) "# Time-steps between two backups: "
    write(9,*) "#---------------------------------"
    write(9,*) "  10000"
    write(9,*) "#---------------------------------"
    close(9)

    open (9, file=trim(directory)//'RESULT_FILE_NAME')
    write(9,*) "#================================"
    write(9,*) "# Specify backup file name here: "
    write(9,*) "#--------------------------------"
    write(9,*) case_name
    write(9,*) "#--------------------------------"
    close(9)

    open (9, file=trim(directory)//'RESULT_FILE_FORMAT')
    write(9,*) "#================================"
    write(9,*) "# Uncomment desired file format: "
    write(9,*) "#--------------------------------"
    write(9,*) "  GMV -> GMV and VisIt format."    
    write(9,*) "# DAT -> Fluent format."
    write(9,*) "# VTK -> Paraview format."
    write(9,*) "#--------------------------------"
    close(9)

    open (9, file=trim(directory)//'RESULT_SAVING_INTERVAL')
    write(9,*) "#========================================"
    write(9,*) "# Time-steps between two result savings: "
    write(9,*) "#----------------------------------------"
    write(9,*) "  10"
    write(9,*) "#----------------------------------------"
    close(9)

  end program
