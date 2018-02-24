!==============================================================================!
  subroutine Info_Mod_Iter_Fill_At(r, c, name_var, n_iter, res)
!------------------------------------------------------------------------------!
!   Inserts infromation at specified position in the information box.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: r         ! row
  integer          :: c         ! column
  character(len=*) :: name_var
  integer          :: n_iter
  real             :: res       
!==============================================================================!

  if (this_proc < 2) then

    ! Normal variables
    if(n_iter >= 0) then
      write(iter_info % lines(r)((c-1)*L_BOX+ 3 :  &
                                 (c-1)*L_BOX+ 5),  '(a3)')  name_var
      write(iter_info % lines(r)((c-1)*L_BOX+ 6 :  &
                                 (c-1)*L_BOX+ 6),  '(a1)')  ':'         
      write(iter_info % lines(r)((c-1)*L_BOX+ 7 :  &
                                 (c-1)*L_BOX+ 9),  '(i3)')  n_iter      
    endif

    ! Residual 
    write(iter_info % lines(r)((c-1)*L_BOX+11 :  &
                               (c-1)*L_BOX+19),  '(1pe9.3)') res

  end if

  end subroutine
