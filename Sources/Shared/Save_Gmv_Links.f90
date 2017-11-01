!==============================================================================!
  subroutine Save_Gmv_Links(sub, NNsub, NCsub, NSsub, NBCsub, NBFsub) 
!------------------------------------------------------------------------------!
!   Creates the file "name.ln.gmv" to check the cell connections.              !
!                                                                              !
!   Links between the computational cells have been introduced as aditional    !
!   cells of general type. Cell centers are introduced as aditional nodes      !
!   with numbers NN+1 to NN+NC. Material of this links is different than       !
!   from the cells, so that they can be visualised  more easily in GMV.        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod 
  use gen_mod 
  use par_mod 
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: sub, NNsub, NCsub, NSsub, NBCsub, NBFsub
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n, c, c1, c2, s 
  integer           :: nf_sub_non_per, nf_sub_per
  character(len=80) :: name_out
!==============================================================================!

  !----------------------!
  !                      !
  !   Create .gmv file   !
  !                      !
  !----------------------!
  name_out = name         

  call Name_File(sub, name_out, '.ln.gmv', len_trim('.ln.gmv'))
  open(9, file=name_out)
  write(6, *) '# Now creating the file:', name_out

  !-----------!
  !   Nodes   !
  !-----------!
  write(9,'(A14)') 'gmvinput ascii'          !  start of GMV file
  write(9,*) 'nodes', NNsub + NCsub + NBCsub + NBFsub

  ! X
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') grid % xn(n)
  end do
  do c=1,NC
    if(NewC(c)  > 0) write(9, '(1PE14.7)') xc(c)
  end do 
  do c=-1,-NBC,-1
    if(NewC(c) /= 0) write(9, '(1PE14.7)') xc(c)
  end do 
  do c=1,NBFsub
    write(9, '(1PE14.7)') xc(BuReIn(c))
  end do

  ! Y
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') grid % yn(n)
  end do
  do c=1,NC
    if(NewC(c)  > 0) write(9, '(1PE14.7)') yc(c)
  end do 
  do c=-1,-NBC,-1
    if(NewC(c) /= 0) write(9, '(1PE14.7)') yc(c)
  end do 
  do c=1,NBFsub
    write(9, '(1PE14.7)') yc(BuReIn(c))
  end do

  ! Z
  do n=1,NN
    if(NewN(n) /= 0) write(9, '(1PE14.7)') grid % zn(n)
  end do
  do c=1,NC
    if(NewC(c)  > 0) write(9, '(1PE14.7)') zc(c)
  end do 
  do c=-1,-NBC,-1
    if(NewC(c) /= 0) write(9, '(1PE14.7)') zc(c)
  end do 
  do c=1,NBFsub
    write(9, '(1PE14.7)') zc(BuReIn(c))
  end do

  !-----------!
  !   Cells   !
  !-----------!

  ! Regular (ordinary) cells
  write(9,*) 'cells', NCsub + NSsub + NBFsub ! + NBFsub
  do c=1,NC
    if(NewC(c)  > 0) then
      if(grid % cells_n_nodes(c) == 8) then
        write(9,*) 'hex 8'
        write(9,*)                                                         &
              NewN(grid % cells_n(1,c)), NewN(grid % cells_n(2,c)),  &
              NewN(grid % cells_n(4,c)), NewN(grid % cells_n(3,c)),  &
              NewN(grid % cells_n(5,c)), NewN(grid % cells_n(6,c)),  &
              NewN(grid % cells_n(8,c)), NewN(grid % cells_n(7,c))
      else if(grid % cells_n_nodes(c) == 6) then
        write(9,*) 'prism 6'
        write(9,*)                                                         &
              NewN(grid % cells_n(1,c)), NewN(grid % cells_n(2,c)),  &
              NewN(grid % cells_n(3,c)), NewN(grid % cells_n(4,c)),  &
              NewN(grid % cells_n(5,c)), NewN(grid % cells_n(6,c))
      else if(grid % cells_n_nodes(c) == 4) then
        write(9,*) 'tet 4'
        write(9,*)                                                         &
              NewN(grid % cells_n(1,c)), NewN(grid % cells_n(2,c)),  &
              NewN(grid % cells_n(3,c)), NewN(grid % cells_n(4,c))
      else if(grid % cells_n_nodes(c) == 5) then
        write(9,*) 'pyramid 5'
        write(9,*)                                                         &
              NewN(grid % cells_n(5,c)), NewN(grid % cells_n(1,c)),  &
              NewN(grid % cells_n(2,c)), NewN(grid % cells_n(4,c)),  &
              NewN(grid % cells_n(3,c))      
      else
        write(*,*) '# Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), 'nodes. '
        write(*,*) '# Exiting !'
        stop
      end if
    end if 
  end do  

  ! Physical links; non-periodic
  nf_sub_non_per = 0 
  do s=1, NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)

    if( (NewS(s) > 0) .and. (NewS(s) <= NSsub) ) then

      if( (Sx(s) * (xc(c2)-xc(c1) )+               &
           Sy(s) * (yc(c2)-yc(c1) )+               &
           Sz(s) * (zc(c2)-zc(c1) ))  > 0.0 ) then 

        nf_sub_non_per = nf_sub_non_per + 1

        c1 = NewC(SideC(1,s))
        c2 = NewC(SideC(2,s))
        if( c2  > 0 ) then
          write(9,*) 'general 1'
          write(9,*) '  2'
          write(9,*) NNsub+c1, NNsub+c2
        else
          write(9,*) 'general 1'
          write(9,*) '  2'
          write(9,*) NNsub+c1, NNsub+NCsub-c2
        end if
      end if

    end if
  end do  

  ! Physical links; periodic
  nf_sub_per    = 0 
  do s=1, NS
    c1 = SideC(1,s)
    c2 = SideC(2,s)

    if( (NewS(s) > 0) .and. (NewS(s) <= NSsub) ) then

      if( (Sx(s) * (xc(c2)-xc(c1) )+               &
           Sy(s) * (yc(c2)-yc(c1) )+               &
           Sz(s) * (zc(c2)-zc(c1) ))  < 0.0 ) then 

        nf_sub_per = nf_sub_per + 1

        c1 = NewC(SideC(1,s))
        c2 = NewC(SideC(2,s))
        if( c2  > 0 ) then
          write(9,*) 'general 1'
          write(9,*) '  2'
          write(9,*) NNsub+c1, NNsub+c2
        else
          write(9,*) 'general 1'
          write(9,*) '  2'
          write(9,*) NNsub+c1, NNsub+NCsub-c2
        end if
      end if

    end if
  end do  

  write(*,*) '# Non-periodic links:', nf_sub_non_per
  write(*,*) '# Periodic links    :', nf_sub_per
  write(*,*) '# Number of sides   :', NSsub

  ! Interprocessor links
  do c=1,NBFsub
    c1 = BuSeIn(c) 
    write(9,*) 'general 1'
    write(9,*) '  2'
    write(9,*) NNsub+c1, NNsub+NCsub+NBCsub+c
  end do  

  !---------------!
  !   Materials   !
  !---------------!
  write(9,*) 'material  4  0'
  write(9,*) 'cells'
  write(9,*) 'links'
  write(9,*) 'periodic'
  write(9,*) 'buffers'

  do c=1, NCsub
    write(9,*) 1
  end do
  do s=1, nf_sub_non_per
    write(9,*) 2
  end do
  do s=1, nf_sub_per
    write(9,*) 3
  end do
  do s=1,NBFsub
    write(9,*) 4
  end do

  write(9,'(A6)') 'endgmv'            !  end the GMV file
  close(9)

  end subroutine Save_Gmv_Links
