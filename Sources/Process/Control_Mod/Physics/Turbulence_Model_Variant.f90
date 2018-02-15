!==============================================================================!
  subroutine Control_Mod_Turbulence_Model_Variant(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENCE_MODEL_VARIANT', 'none',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if( val .ne. 'NONE'     .and.  &
      val .ne. 'HYBRID'   .and.  &
      val .ne. 'PURE'     .and.  &
      val .ne. 'URANS'    .and.  &
      val .ne. 'LOW_RE'   .and.  &
      val .ne. 'HIGH_RE'  .and.  &
      val .ne. 'WALE'     .and.  &
      val .ne. 'DYNAMIC'  .and.  &
      val .ne. 'SMAGORINSKY') then
    print *, '# Unknown turbulence model variant: ', trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine
