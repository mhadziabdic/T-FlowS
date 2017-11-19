!==============================================================================!
  subroutine Save_Gmv_Faces(grid, sub, NNsub)
!------------------------------------------------------------------------------!
! Writes .faces.gmv file.                                                      !
!----------------------------------[Modules]-----------------------------------!
  use all_mod, only: bcmark
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, NNsub
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c1, c2, n, s
  character(len=80) :: name_out
!==============================================================================!

  !-----------------------------------------!
  !                                         !
  !   Create boundary condition .gmv file   !
  !                                         !
  !-----------------------------------------!
  if(sub /= 0) return

  call Name_File(sub, name_out, '.faces.gmv', len_trim('.faces.gmv'))
  open(9, file=name_out)
  write(6, *) '# Now creating the file:', trim(name_out)

  !-----------!
  !   Start   !
  !-----------!
  write(9,'(A14)') 'gmvinput ascii'

  !-----------!
  !   Nodes   !
  !-----------!
  write(9,*) 'nodes', NNsub

  do n = 1, grid % n_nodes
    write(9, '(1PE14.7)') grid % xn(n)
  end do
  do n = 1, grid % n_nodes
    write(9, '(1PE14.7)') grid % yn(n)
  end do
  do n = 1, grid % n_nodes
    write(9, '(1PE14.7)') grid % zn(n)
  end do

  !-----------!
  !   Cells   !
  !-----------!

  ! Count the cell faces on the periodic boundaries
  write(9,*) 'cells', grid % n_faces
  do s = 1, grid % n_faces
    if(grid % faces_n_nodes(s) == 4) then
      write(9,*) 'quad 4'
      write(9,'(4I9)')                             &
        grid % faces_n(1,s), grid % faces_n(2,s),  &
        grid % faces_n(3,s), grid % faces_n(4,s)
    else if(grid % faces_n_nodes(s) == 3) then
      write(9,*) 'tri 3'
      write(9,'(3I9)')                             &
        grid % faces_n(1,s), grid % faces_n(2,s),  &
        grid % faces_n(3,s)
    else
      write(*,*) '# Unsupported cell type ',       &
                 grid % faces_n_nodes(s), ' nodes.'
      write(*,*) '# Exiting'
      stop 
    end if
  end do  

  !---------------!
  !   Materials   !
  !---------------!
  write(9,*) 'materials', grid % n_boundary_conditions + 1, 0
  do n = 1, grid % n_boundary_conditions
    write(9,*) grid % boundary_conditions(n) % name
  end do        
  write(9,*) 'DEFAULT_INSIDE'

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
   
    ! If boundary 
    if( c2 < 0 ) then 
      write(9,*) bcmark(c2) 

    ! If inside 
    else 
      write(9,*) grid % n_boundary_conditions + 1 

    end if
  end do

  write(9,'(A6)') 'endgmv'            !  end the GMV file
  close(9)

  end subroutine
