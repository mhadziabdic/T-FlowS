!==============================================================================!
  subroutine CalcVort(grid)
!------------------------------------------------------------------------------!
!  Computes the magnitude of the vorticity                                     !
!------------------------------------------------------------------------------!
!  Vort = sqrt( 2 * Sij * Sij )                                                !
!  Sij = 1/2 ( dUi/dXj - dUj/dXi )                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i 
  real    :: Sjk
!==============================================================================!

  call Exchange(grid, U % n)
  call Exchange(grid, V % n)
  call Exchange(grid, W % n)

  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  do i = 1, 3        

    if(i == 1) then
      call GraPhi(grid, W % n, 2, PHIy, .TRUE.)  ! dW/dy
      call GraPhi(grid, V % n, 3, PHIz, .TRUE.)  ! dV/dz
    end if
    if(i == 2) then
      call GraPhi(grid, W % n, 1, PHIx, .TRUE.)  ! dW/dx
      call GraPhi(grid, U % n, 3, PHIz, .TRUE.)  ! dU/dz
    end if
    if(i == 3) then
      call GraPhi(grid, V % n, 1, PHIx, .TRUE.)  ! dV/dx
      call GraPhi(grid, U % n, 2, PHIy, .TRUE.)  ! dU/dy
    end if

    do c = 1, grid % n_cells
      if(i == 1) then
        Sjk = 0.5*(PHIy(c)-PHIz(c)) ! Syz=Szy:  .5(dV/dz+dW/dy)
        Vort(c) = 2.0*Sjk*Sjk   
      end if
      if(i == 2) then
        Sjk = 0.5*(PHIx(c)-PHIz(c)) ! Sxz=Szx:  .5(dU/dz+dW/dx)
        Vort(c) = Vort(c) + 2.0*Sjk*Sjk  
      end if
      if(i == 3) then
        Sjk = 0.5*(PHIx(c)-PHIy(c)) ! Sxy=Syx:  .5(dV/dx+dU/dy)
        Vort(c) = Vort(c) + 2.0*Sjk*Sjk 
      end if
    end do 

  end do  ! i

  do c = 1, grid % n_cells
    Vort(c) = sqrt(abs(2.0 * Vort(c)))
  end do

  end subroutine
