!==============================================================================!
  subroutine Control_Mod_Advection_Scheme_For_Turbulence(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('ADVECTION_SCHEME_FOR_TURBULENCE', 'upwind', &
                                   val, verbose)
  call To_Upper_Case(val)

  if(val .ne. 'YES'       .and.  &
     val .ne. 'NO'        .and.  &
     val .ne. 'UDS'       .and.  &
     val .ne. 'CDS'       .and.  &
     val .ne. 'LUDS'      .and.  &
     val .ne. 'QUICK'     .and.  &
     val .ne. 'SMART'     .and.  &
     val .ne. 'GAMMA'     .and.  &
     val .ne. 'MINMOD'    .and.  &
     val .ne. 'SUPERBEE'  .and.  &
     val .ne. 'AVL_SMART') then

    print *, '# Unknown advection scheme for turbulence: ', trim(val)
    print *, '# Exiting!'
    stop 

  end if

  end subroutine
