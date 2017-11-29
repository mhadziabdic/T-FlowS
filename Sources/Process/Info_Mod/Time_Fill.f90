!==============================================================================!
  subroutine Info_Mod_Time_Fill(n, sim_time, wall_time)
!------------------------------------------------------------------------------!
!   Fills the info box with information to be written on the screen.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
! use pro_mod    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: n          ! time step          
  real    :: sim_time   ! simulation time
  real    :: wall_time  ! number of seconds of wall-clock time
!-----------------------------------[Locals]-----------------------------------!
  integer :: hours, minutes, seconds
!==============================================================================!

  ! Write time step number
  write(time_info % lines(2)(31:41), '(a11)')    'Time step :'
  write(time_info % lines(2)(43:48),  '(i6)')    n

  ! Write simulation time
  write(time_info % lines(4)(25:41),    '(a17)') 'Simulation time :'
  write(time_info % lines(4)(43:51), '(1pe9.3)') sim_time
  write(time_info % lines(4)(53:55),     '(a3)') '[s]'     

  ! Write wall-clock time
  write(time_info % lines(5)(25:41), '(a17)') 'Wall-clock time :'
  hours   = floor(  wall_time  / 3600.0 )
  minutes = floor( (wall_time - 3600.0 * hours) / 60.0)
  seconds = floor(  wall_time - 3600.0 * hours - 60.0 * minutes )
  write(time_info % lines(5)(43:45), '(i3.3)')  hours
  write(time_info % lines(5)(46:46),   '(a1)')  ':'
  write(time_info % lines(5)(47:48), '(i2.2)')  minutes
  write(time_info % lines(5)(49:49),   '(a1)')  ':'
  write(time_info % lines(5)(50:51), '(i2.2)')  seconds
  write(time_info % lines(5)(53:57),   '(a5)') '[hms]'     

  end subroutine
