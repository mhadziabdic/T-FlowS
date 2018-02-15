!==============================================================================!
  subroutine Control_Mod_Time_Integration_For_Cross_Diffusion(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('TIME_INTEGRATION_FOR_CROSS_DIFFUSION',  &
                                  'fully_implicit',                        &
                                   val,                                    &
                                   verbose)
  call To_Upper_Case(val)

  if( val.ne.'FULLY_IMPLICIT'  .and.  &
      val.ne.'ADAMS_BASHFORTH' .and.  &
      val.ne.'CRANK_NICOLSON') then
    print *, '# Unknown time-integration scheme for cross-diffusion: ',  &
             trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine