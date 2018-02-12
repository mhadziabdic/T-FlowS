!==============================================================================!
  subroutine Control_Mod_Turbulence_Model(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENCE_MODEL', 'k_eps',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if( val .ne. 'K_EPS'    .and.  &
      val .ne. 'K_EPS_VV' .and.  &
      val .ne. 'ZETA'     .and.  &
      val .ne. 'HYB_ZETA' .and.  &
      val .ne. 'LES'      .and.  &
      val .ne. 'DES_SPA'  .and.  &
      val .ne. 'SPA_ALL') then
    print *, '# Unknown turbulence model: ', trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine
