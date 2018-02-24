!==============================================================================!
  subroutine Info_Mod_Bulk_Fill(courant, peclet, fx, fy, fz, px, py, pz)
!------------------------------------------------------------------------------!
!   Fills the info box with information to be written on the screen.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: courant, peclet, fx, fy, fz, px, py, pz
!==============================================================================!

  if (this_proc < 2) then

    ! Courant and Peclet numbers
    write(bulk_info % lines(1)( 5:27),    '(a23)') 'Maximum Courant number:'
    write(bulk_info % lines(1)(29:37), '(1pe9.3)') courant

    write(bulk_info % lines(1)(45:66),    '(a22)') 'Maximum Peclet number:'
    write(bulk_info % lines(1)(68:76), '(1pe9.3)') peclet 

    write(bulk_info % lines(2)( 5:12),     '(a8)') 'Flux x :'
    write(bulk_info % lines(2)(14:22), '(1pe9.2)') fx
    write(bulk_info % lines(2)(32:39),     '(a8)') 'Flux y :'
    write(bulk_info % lines(2)(41:49), '(1pe9.2)') fy
    write(bulk_info % lines(2)(59:66),     '(a8)') 'Flux z :'
    write(bulk_info % lines(2)(68:76), '(1pe9.2)') fz

    write(bulk_info % lines(3)( 5:12),     '(a8)') 'Pdrop x:'
    write(bulk_info % lines(3)(14:22), '(1pe9.2)') px
    write(bulk_info % lines(3)(32:39),     '(a8)') 'Pdrop y:'
    write(bulk_info % lines(3)(41:49), '(1pe9.2)') py
    write(bulk_info % lines(3)(59:66),     '(a8)') 'Pdrop z:'
    write(bulk_info % lines(3)(68:76), '(1pe9.2)') pz

  end if

  end subroutine
