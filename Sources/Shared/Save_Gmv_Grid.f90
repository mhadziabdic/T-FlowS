!==============================================================================!
  subroutine Save_Gmv_Grid(grid, sub, NNsub, NCsub)
!------------------------------------------------------------------------------!
! Writes: NAME.gmv, NAME.faces.gmv, NAME.shadow.gmv                            !
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod
  use par_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: sub, NNsub, NCsub, NmaterBC
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c,  c1,  c2,  n, s
  character(len=80) :: name_out
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .gmv file   !
  !                      !
  !----------------------!
  call Name_File(sub, name_out, '.gmv', len_trim('.gmv'))
  open(9, file=name_out)
  write(6, *) '# Now creating the file:', name_out

  !-----------!
  !   Start   !
  !-----------!
  write(9,'(A14)') 'gmvinput ascii'

  !-----------!
  !   Nodes   !
  !-----------!
  write(9,*) 'nodes', NNsub

  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') grid % xn(n)
  end do
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') grid % yn(n)
  end do
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') grid % zn(n)
  end do

  !-----------!
  !   Cells   !
  !-----------!
  write(9,*) 'cells', NCsub
  do c=1,NC
    if(NewC(c) /= 0) then
      if(grid % cells_n_nodes(c) == 8) then
        write(9,*) 'hex 8'
        write(9,'(8I9)')                                          &
          NewN(grid % cells_n(1,c)), NewN(grid % cells_n(2,c)),   &
          NewN(grid % cells_n(4,c)), NewN(grid % cells_n(3,c)),   &
          NewN(grid % cells_n(5,c)), NewN(grid % cells_n(6,c)),   &
          NewN(grid % cells_n(8,c)), NewN(grid % cells_n(7,c))
      else if(grid % cells_n_nodes(c) == 6) then
        write(9,*) 'prism 6'
        write(9,'(6I9)')                                          &
          NewN(grid % cells_n(1,c)), NewN(grid % cells_n(2,c)),   &
          NewN(grid % cells_n(3,c)), NewN(grid % cells_n(4,c)),   &
          NewN(grid % cells_n(5,c)), NewN(grid % cells_n(6,c))
      else if(grid % cells_n_nodes(c) == 4) then
        write(9,*) 'tet 4'
        write(9,'(4I9)')                                          &
          NewN(grid % cells_n(1,c)), NewN(grid % cells_n(2,c)),   &
          NewN(grid % cells_n(3,c)), NewN(grid % cells_n(4,c))
      else if(grid % cells_n_nodes(c) == 5) then
        write(9,*) 'pyramid 5'
        write(9,'(5I9)')                                          &
          NewN(grid % cells_n(5,c)), NewN(grid % cells_n(1,c)),   &
          NewN(grid % cells_n(2,c)), NewN(grid % cells_n(4,c)),   &
          NewN(grid % cells_n(3,c))
      else
        write(*,*) '# Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), ' nodes.'
        write(*,*) '# Exiting'
        stop 
      end if 
    end if
  end do  

  !---------------!
  !   Materials   !
  !---------------!
  write(9,'(A10,2I5)') 'materials', Nmat, 0
  do n=1,1024
    if(Mater(n)) write(9,*) n 
  end do        
  do c=1,NC
    if(NewC(c) /= 0) then
      write(9,*) material(c)
    end if
  end do        

  write(9,'(A6)') 'endgmv'            !  end the GMV file
  close(9)

  !-----------------------------------------!
  !                                         !
  !   Create boundary condition .gmv file   !
  !                                         !
  !-----------------------------------------!
  if(sub /= 0) return

  call Name_File(sub, name_out, '.faces.gmv', len_trim('.faces.gmv'))
  open(9, file=name_out)
  write(6, *) '# Now creating the file:', name_out

  !-----------!
  !   Start   !
  !-----------!
  write(9,'(A14)') 'gmvinput ascii'

  !-----------!
  !   Nodes   !
  !-----------!
  write(9,*) 'nodes', NNsub

  do n=1,NN
    write(9, '(1PE14.7)') grid % xn(n)
  end do
  do n=1,NN
    write(9, '(1PE14.7)') grid % yn(n)
  end do
  do n=1,NN
    write(9, '(1PE14.7)') grid % zn(n)
  end do

  !-----------!
  !   Cells   !
  !-----------!

  ! Count the cell faces on the periodic boundaries
  write(9,*) 'cells', NS
  do s=1,NS
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
                 grid % cells_n_nodes(c), ' nodes.'
      write(*,*) '# Exiting'
      stop 
    end if
  end do  

  !---------------!
  !   Materials   !
  !---------------!
  NmaterBC = 0
  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    if( c2 < 0 ) then 
      NmaterBC=max(NmaterBC,bcmark(c2))
    end if
  end do  

  write(9,*) 'materials', NmaterBC + 1, 0
  do n = 1, NmaterBC + 1
    write(9,*) n 
  end do        

  do s=1,NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)
    if( c2 < 0 ) then 
      write(9,*) bcmark(c2) 
    else 
      write(9,*) NmaterBC+1 
    end if
  end do

  write(9,'(A6)') 'endgmv'            !  end the GMV file
  close(9)

  !-----------------------------!
  !                             !
  !   Create shadow .gmv file   !
  !                             !
  !-----------------------------!
  if(sub /= 0) return

  call Name_File(sub, name_out, '.shadow.gmv', len_trim('.shadow.gmv'))
  open(9, file=name_out)
  write(6, *) '# Now creating the file:', name_out

  do s=NS+1,NS+NSsh
    write(9,*) grid % faces_n_nodes(s) 
    if(grid % faces_n_nodes(s)==3) then
      write(9,*) grid % faces_n(1,s),  &
                 grid % faces_n(2,s),  &
                 grid % faces_n(3,s),  &
                 SideC(1,s), SideC(2,s) 
    else if(grid % faces_n_nodes(s)==4) then
      write(9,*) grid % faces_n(1,s),  &
                 grid % faces_n(2,s),  &
                 grid % faces_n(3,s),  &
                 grid % faces_n(4,s),  &
                 SideC(1,s), SideC(2,s) 
    end if
  end do  

  close(9)

  end subroutine Save_Gmv_Grid
