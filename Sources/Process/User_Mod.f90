!==============================================================================!
  module User_Mod
!------------------------------------------------------------------------------!
!   This is embrio of a future User module, a place where user can
!   define his/her variables and pass them around his functions
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none

  type User_Type

  end type 

  contains

  include 'User_Mod/Allocate.f90'
  include 'User_Mod/Before_Exit.f90'
  include 'User_Mod/Begining_Of_Time_Step.f90'
  include 'User_Mod/Initialize.f90'
  include 'User_Mod/End_Of_Time_Step.f90'
  include 'User_Mod/Force.f90'
  include 'User_Mod/Save_Results.f90'
  include 'User_Mod/Source.f90'

  end module 


