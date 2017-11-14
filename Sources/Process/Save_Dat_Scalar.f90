!==============================================================================!
  subroutine Save_Dat_Scalar(grid, namSc, idSc, phi)
!------------------------------------------------------------------------------!
!   Writes: NAME.dat                                                           !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use pro_mod
  use par_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!---------------------------------[Arguments]----------------------------------!
  integer   :: idSc
  character :: namSc*(*)
  real      :: phi(-grid % n_boundary_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer   :: N, s, c, c1, c2, Nfac(10), NtotFac
!==============================================================================!

  !------------!
  !   Inside   !
  !------------!
  write(9,'(A4,A,A2)') '(0 "', namSc, '")'
  write(9,'(A6,I3,A9,I8,2X,I9,A2)')  '(300 (', idSc, ' 1 1 0 0 ',  &
                                       1, grid % n_cells, ')(' 
  do c=1,grid % n_cells
    write(9,'(F14.6)') phi(c)
  end do  
  write(9,'(A2)') '))'

  !---------------------!
  !   On the boundary   !
  !---------------------!
  NtotFac = 0
  do n=1,10   ! browse through boundary condition types
    Nfac(n) = 0
    do s = 1, grid % n_faces   ! count the faces with boundary condition "n"
      c2 = SideC(2,s)
      if(c2 < 0) then
        if(BCmark(c2) == n) Nfac(n)=Nfac(n)+1
      end if
    end do    ! sides

    if(Nfac(n) > 0) then
      write(9,'(A4,A,A17,I3,A3)') '(0 "', namSc, ' on the boundary ', n, ' ")'
      write(9,'(A6,I3,I4,A7,I7,I7,A4)') '(300 (', idSc, 100+n, ' 1 0 0 ', &
                                        NtotFac+1, NtotFac+Nfac(n), ')('
      do s = 1, grid % n_faces !@@@ +NSsh
        c1 = SideC(1,s)
        c2 = SideC(2,s)
        if(c2 < 0) then
          if(BCmark(c2) == n) then
            if(TypeBC(c2) == SYMMETRY) then
              write(9,'(F14.6)') phi(c1)
            else
              write(9,'(F14.6)') phi(c2)
            end if
          end if
        end if
      end do
      write(9,'(A2)') '))'
    end if

    ! Prepare for next boundary
    NtotFac = NtotFac+Nfac(n)

  end do   ! n -> boundary condition types

  end subroutine
