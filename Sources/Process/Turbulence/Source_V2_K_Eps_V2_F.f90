!==============================================================================!
  subroutine Source_V2_K_Eps_V2_F(grid, Nstep)
!------------------------------------------------------------------------------!
!   Calculate source terms in equation for v2                                  !
!   Term which is negative is put on left hand side in diagonal of             !
!   martix of coefficient                                                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: Nstep
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c
  integer           :: s, c1, c2,j     
  character(len=80) :: turbulence_model
  character(len=80) :: turbulence_model_variant
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

  call Control_Mod_Turbulence_Model(turbulence_model)

  if(turbulence_model == 'ZETA' .or.  &
     turbulence_model=='HYB_ZETA') then

    ! Positive source term 
    ! The first option in treating the source is making computation very 
    ! sensitive to initial condition while the second one can lead to 
    ! instabilities for some cases such as flow around cylinder. That is why we 
    ! choose this particular way to the add source term. 
    do c = 1, grid % n_cells
      if(Nstep > 500) then
        b(c) = b(c) + f22%n(c)*grid % vol(c)
      else
        b(c) = b(c) + max(0.0,f22 % n(c)*grid % vol(c))
        A % val(A % dia(c)) = A % val(A % dia(c))               &
                            + max(0.0, -f22 % n(c) * grid % vol(c)  &
                            / (v_2 % n(c) + TINY))    
      end if      
      A % val(A % dia(c)) =  A % val(A % dia(c))  &
                          + grid % vol(c) * p_kin(c)     &
                          / (kin % n(c)+TINY) 
    end do
  else if(turbulence_model == 'K_EPS_VV') then
    do c = 1, grid % n_cells
      b(c) = b(c) + max(0.0, f22 % n(c) * kin % n(c) * grid % vol(c))
      A % val(A % dia(c)) = A % val(A % dia(c))                            &
                          + max(0.0, -f22 % n(c) * kin % n(c) * grid % vol(c)  &
                          / (v_2 % n(c) + TINY))
    end do
    do c = 1, grid % n_cells  
      A % val(A % dia(c)) = A % val(A % dia(c))  &
                          + grid % vol(c) * eps % n(c) / (kin%n(c)+TINY) 
    end do
  end if

  end subroutine
