!==============================================================================!
  subroutine Source_Zeta_K_Eps_Zeta_F(grid, n_step)
!------------------------------------------------------------------------------!
!   Calculate source terms in equation for v2                                  !
!   Term which is negative is put on left hand side in diagonal of             !
!   martix of coefficient                                                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: n_step
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c
  integer           :: s, c1, c2,j     
!==============================================================================!
!   In transport equation for v2 two source terms exist which have form:       !
!                                                                              !
!     /                     /                                                  !
!    |                     |                                                   !
!    | f22 * kin * dV  -   | (vi2 * eps / kin) * dV                            !
!    |                     |                                                   !
!   /                      /                                                   !
!                                                                              !
!   First term can appire as pozitiv and as negative as well so depend of      !
!   sign of term , it is placed on left or right hand side.  Second, negative  !
!   source term is added to main diagonal left hand side coefficient matrix    !
!   in order to increase stability of solver                                   !
!------------------------------------------------------------------------------!      

  if(turbulence_model == K_EPS_ZETA_F .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F) then

    ! Positive source term 
    ! The first option in treating the source is making computation very 
    ! sensitive to initial condition while the second one can lead to 
    ! instabilities for some cases such as flow around cylinder. That is why we 
    ! choose this particular way to the add source term. 
    do c = 1, grid % n_cells
      if(n_step > 500) then
        b(c) = b(c) + f22 % n(c) * grid % vol(c)
      else
        b(c) = b(c) + max(0.0,f22 % n(c)*grid % vol(c))
        A % val(A % dia(c)) = A % val(A % dia(c))               &
                            + max(0.0, -f22 % n(c) * grid % vol(c)  &
                            / (zeta % n(c) + TINY))    
      end if      
      A % val(A % dia(c)) =  A % val(A % dia(c))  &
                          + grid % vol(c) * p_kin(c)     &
                          / (kin % n(c)+TINY) 
    end do

  end if

  end subroutine
