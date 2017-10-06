!======================================================================!
  subroutine CalcVort()
!----------------------------------------------------------------------!
!  Computes the magnitude of the vorticity                             !
!----------------------------------------------------------------------!
!  Vort = sqrt( 2 * Sij * Sij )                                       !
!  Sij = 1/2 ( dUi/dXj - dUj/dXi )                                     !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use pro_mod
  use les_mod
  use rans_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer :: c, i 
  real    :: Sii, Sjk
!======================================================================!

  call Exchng(U % n)
  call Exchng(V % n)
  call Exchng(W % n)

!===================!
!                   !
!     SGS terms     !
!                   !
!===================!
  do i=1,3        

    if(i == 1) then
      call GraPhi(W % n,2,PHIy, .TRUE.)  ! dW/dy
      call GraPhi(V % n,3,PHIz, .TRUE.)  ! dV/dz
    end if
    if(i == 2) then
      call GraPhi(W % n,1,PHIx, .TRUE.)  ! dW/dx
      call GraPhi(U % n,3,PHIz, .TRUE.)  ! dU/dz
    end if
    if(i == 3) then
      call GraPhi(V % n,1,PHIx, .TRUE.)  ! dV/dx
      call GraPhi(U % n,2,PHIy, .TRUE.)  ! dU/dy
    end if

    do c=1,NC
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

  do c = 1, NC
    Vort(c) = sqrt(abs(2.0 * Vort(c)))
  end do
  end subroutine CalcVort
