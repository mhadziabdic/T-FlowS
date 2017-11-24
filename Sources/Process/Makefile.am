bin_PROGRAMS = Process

Process_SOURCES =	 			\
  Allocate_Memory.f90				\
  Allocate_Variables.f90			\
  Balance_Mass.f90				\
  Update_Boundary_Values.f90			\
  Convective_Outflow.f90			\
  Compute_Gradient_Matrix.f90			\
  CalcMn_Cylind.f90				\
  Compute_Mean.f90				\
  Compute_Pressure_Fractional.f90		\
  Compute_Pressure_Simple.f90			\
  Compute_Shear_And_Vorticity.f90		\
  CalcVort.f90					\
  Compute_Minimum_Maximum.f90			\
  Compute_Fluxes.f90				\
  Compute_Geometry.f90				\
  Compute_Scalar.f90				\
  Advection_Scheme.f90				\
  Correct_Bad.f90				\
  Correct_Velocity.f90				\
  Find_Bad.f90					\
  Pressure_Matrix_Fractional.f90		\
  GraCorNew.f90					\
  GradP3.f90					\
  GradP.f90					\
  GraPhiCor.f90					\
  GraPhi.f90					\
  IniPar.f90					\
  Initialize_Variables.f90			\
  Load_Boundary_Conditions.f90			\
  Load_Cns.f90					\
  Load_Ini.f90					\
  Load_Restart.f90				\
  Load_Restart_Ini.f90				\
  Logo_Process.f90				\
  Main.f90   					\
  NearWallCell.f90				\
  NewUVW.f90					\
  Read_Problem.f90				\
  Save_Dat_Results.f90				\
  Save_Dat_Scalar.f90				\
  Save_Ini.f90					\
  Save_Gmv_Results.f90				\
  Save_Restart.f90				\
  ReaCom.f90					\
  SavParView.f90				\
  Time_And_Length_Scale.f90

SUBDIRS = Communication				\
          Solvers				\
          Turbulence				\
          User					\
          ../Shared

LDADD = Communication/libCommunication.a	\
        Solvers/libSolvers.a			\
        Turbulence/libTurbulence.a		\
        User/libUser.a				\
        ../Shared/libShared.a

INCLUDE = -I@top_srcdir@/
