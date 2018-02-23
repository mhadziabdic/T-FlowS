!==============================================================================!
  subroutine Cgns_Mod_Get_Arrays_Dimensions(idx, nn_or_nc)
!------------------------------------------------------------------------------!
!   Fetches correct dimensions for arrays in CGNS lib dependent functions      !
!------------------------------------------------------------------------------!
!   Arrays structure in CGNS parallel functions are strictly followings:       !
!   Processor:    |        P_1        |               P_2               | ...  !
!   x,y,z:        |      (1 : NN_1)   |       NN_1 + 1 : NN_1 + NN_2    | ...  !
!------------------------------------------------------------------------------!
!   Cell type:    |      HEXA_8      |     PENTA_6      |       PYRA_5     |...!
!   Connections:  |-p1-|-p2-|...|-pN-|-p1-|-p2-|...|-pN-|-p1-|-p2-|...|-pN-|...!
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use par_mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: idx      !(out)
  integer :: nn_or_nc !(in )
!-----------------------------------[Locals]-----------------------------------!
  integer :: tmp(1:n_proc)
!------------------------------------------------------------------------------!

  ! single-processor case
  if(this_proc .eq. 0) then
    idx = 1
    return
  end if

  ! multi-processor case
  tmp = 0
  call wait
  tmp(this_proc) = nn_or_nc

  call IglSumArray(tmp, n_proc)

  idx = 1

  if (this_proc > 1) idx = sum(tmp(1:this_proc-1)) + 1

  end subroutine
