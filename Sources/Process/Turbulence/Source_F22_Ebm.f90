!==============================================================================!
  subroutine SourceF22_EBM(grid)
!------------------------------------------------------------------------------!
!   Calculate source terms in eliptic relaxation equation                      !
!   and imposing  boundary condition forf22                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Moules]-----------------------------------!
  use all_mod
  use pro_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer s, c, c1, c2, j
  real    Sor11,  f22hg
  real    A0
!==============================================================================!
!                                                                              !
!   The form of source terms are :                                             !
!                                                                              !
!     /                                                                        !
!    |                                                                         !
!    | f22hg*dV ; f22hg - f22hg homogenious is placed in a source              !
!    |                     coefficients b(c)                                   !
!   /                                                                          !
!      f22hg = (1.0 - Cv_1)*(vi2(c)/Kin(c) - 2.0/3.0)/Tsc(c)     &             !
!              + 2.0*Cv_2*Pk(c)/(3.0*Kin(c))                                   !
!                                                                              !
!     /                                                                        !
!    |                                                                         !
!    | f22*dV ;this term is placed in a diagonal of coefficient matrice        !
!    |                                                                         !
!   /                                                                          !
!                                                                              !
!                                                                              !
!  Dimensions of certain variables                                             !
!                                                                              !
!     Tsc            [s]                                                       !
!     Kin            [m^2/s^2]                                                 !
!     Eps            [m^3/s^2]                                                 !
!     vi2            [m^2/s^2]                                                 !
!     f22            [-]                                                       !
!     Lsc            [m]                                                       !
!                                                                              !
!------------------------------------------------------------------------------!

  call Time_And_Length_Scale(grid)

  ! Source term f22hg
  do c = 1, grid % n_cells
    f22hg = 1.0
    Sor11 = grid % vol(c)/Lsc(c)**2
    A % val(A % dia(c)) = A % val(A % dia(c)) + Sor11     
    b(c) = b(c) + f22hg*grid % vol(c)/Lsc(c)**2
  end do

  ! Source term
  do s = 1, grid % n_faces
    c1=grid % faces_c(1,s)
    c2=grid % faces_c(2,s)
    if(c2 < 0 .and. TypeBC(c2) /= BUFFER ) then
      if(TypeBC(c2)==WALL .or. TypeBC(c2)==WALLFL) then

          f22 % n(c2) = 0.0

      end if   ! end if of BC=wall
    end if    ! end if of c2<0
  end do

  end subroutine
