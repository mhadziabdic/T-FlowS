!======================================================================!
  subroutine Find_Bad(grid)
!----------------------------------------------------------------------!
! Searches for cells which are "bad" for calculation of pressure       !
! gradients.                                                           !
!                                                                      !
! Practically, these are the tetrahedronal cells with two faces on the !
! boundary and two in the domain.                                      ! 
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use Flow_Mod
  use Comm_Mod
  use Grid_Mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  type(Grid_Type) :: grid
!-------------------------------[Locals]-------------------------------!
  integer :: s, c, c1, c2, n_bad
!======================================================================!

  BadForG = .false. 
  NumGood = 0

  n_bad = 0

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      NumGood(c1) = NumGood(c1) + 1 
    end if  
  end do

  do c = 1, grid % n_cells
    if(NumGood(c)==2) then
      BadForG(c) = .true.
      n_bad = n_bad + 1
    end if
  end do 

  NumGood = 0
  
  call Comm_Mod_Global_Sum_Int(n_bad)

  if(this_proc < 2) print *, '# There are ', n_bad, &
                          ' bad cells for gradients.'

  end subroutine
