!==============================================================================!
  subroutine Number
!------------------------------------------------------------------------------!
!   Number the cells in each subdomain for subsequent separate saving.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
  use gen_mod 
  use par_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer              :: b, c, n, s, c1, c2, sub, subo, ln
  integer              :: n_nodes_sub, n_cells_sub, n_faces_sub,  &
                          n_b_cells_sub,n_buff_sub,NCSsub,NCFsub
  character(len=80)    :: name_out
  integer,allocatable  :: side_cell(:,:)
!==============================================================================!
!   Each subdomain needs two buffers: a send buffer and a receive buffer.
!   A receive buffer will be stored as aditional boundary cells for each
!   subdomain. So each subdomain will have NBC physical boundary faces and 
!   NBBC-NBC buffer bounndary cells. It is handy to do it that way, because 
!   most of the algorythms can remain the same as they are now.  They won't 
!   even "know" that they use values from other processors.  On the other 
!   hand, a sending buffer has to be allocated in a new separate array called 
!   simply buffer(). An additional array is needed to keep track of all the 
!   indexes. That one is called buffind().  buffind() has stored cell numbers 
!   from it's own subdomain so that later they can be copied with (well, 
!   something like that):
!
!   do i=1,BUFFSIZ
!     buffer(i) = U(buffind(i))
!   end do
!------------------------------------------------------------------------------!

  !-------------------------------!
  !                               !
  !   Browse through subdomains   !
  !                               !
  !-------------------------------!
  do sub=1,n_sub

    call Name_File(sub, name_out, '.buf', len_trim('.buf'))
    open(9, file=name_out)
    write(*, *) '# Now creating the file:', name_out

    write(9,'(A20)') '%%%%%%%%%%%%%%%%%%%%'
    write(9,'(A20)') '%                  %'
    write(9,'(A20)') '%  Buffer indexes  %'
    write(9,'(A20)') '%                  %'
    write(9,'(A20)') '%%%%%%%%%%%%%%%%%%%%'

    ! Cells
    n_cells_sub = 0     ! number of cells in subdomain
    do c=1,NC
      NewC(c)=0
    end do
    do c=1,NC
      if(proces(c) == sub) then
        n_cells_sub=n_cells_sub+1
        NewC(c)=n_cells_sub
      end if
    end do

    ! Nodes
    n_nodes_sub = 0     ! number of cells in subdomain
    do n=1,NN
      NewN(n)=0
    end do
    do c=1,NC
      if(proces(c) == sub) then
        do ln=1,grid % cells(c) % n_nodes
          NewN(grid % cells(c) % n(ln))=-1
        end do
      end if
    end do
    do n=1,NN
      if(NewN(n) == -1) then
        n_nodes_sub=n_nodes_sub+1
        NewN(n)=n_nodes_sub
      end if
    end do

    ! Faces & real boundary cells
    n_faces_sub   = 0  ! number of sides in subdomain
    n_b_cells_sub = 0  ! number of real boundary cells in subdomain
    NCSsub = 0
    do s=1,NS
      NewS(s)=0
    end do
    do c=-NBC,-1
      NewC(c)=0
    end do
    do s=1,NS
      c1=SideC(1,s)  
      c2=SideC(2,s) 
      if(c2  > 0) then
        if( (proces(c1) == sub).and.(proces(c2) == sub) ) then
          n_faces_sub=n_faces_sub+1
          NewS(s)=n_faces_sub
        end if
      else ! c2 < 0
        if( proces(c1) == sub ) then
          n_faces_sub =n_faces_sub+1
          NewS(s)=n_faces_sub     ! new number for the side
          n_b_cells_sub=n_b_cells_sub+1
          NewC(c2)=-n_b_cells_sub ! new number for the boundary cell
        end if
      end if 
    end do

    do s=1,n_copy
      c1=CopyS(1,s)
      c2=CopyS(2,s)
      if( (proces(c1) == sub).and.(proces(c2) == sub) ) then
        NCSsub=NCSsub+1
      end if
    end do

    write(*,*) 'Now saving subdomain ', sub, ' with:'
    write(*,*) n_cells_sub, ' cells'
    write(*,*) n_nodes_sub, ' nodes' 
    write(*,*) n_faces_sub, ' sides' 
    write(*,*) n_b_cells_sub, ' physical boundary cells' 

    !--------------------!
    !   Create buffers   !
    !--------------------!
    n_buff_sub = 0
    NCFsub = 0
    write(9,'(A33)') '#--- Number of physical b. cells:'
    write(9,'(I8)')  n_b_cells_sub   
    do subo=1,n_sub
      if(subo /= sub) then
        NBBs(subo)=n_buff_sub+1

        ! Faces inside the domain
        do s=1,NS
          c1=SideC(1,s)  
          c2=SideC(2,s) 
          if(c2  > 0) then
            if( (proces(c1) == sub).and.(proces(c2) == subo) ) then
              n_buff_sub = n_buff_sub+1
              BuSeIn(n_buff_sub)=NewC(c1) ! Buffer Send Index 
              BuReIn(n_buff_sub)=c2       ! important for coordinate
              BufPos(n_buff_sub)=-n_b_cells_sub-n_buff_sub

              NewS(s)=n_faces_sub+n_buff_sub
            end if
            if( (proces(c2) == sub).and.(proces(c1) == subo) ) then
              n_buff_sub = n_buff_sub+1
              BuSeIn(n_buff_sub)=NewC(c2) ! Buffer Send Index
              BuReIn(n_buff_sub)=c1       ! important for coordinate
              BufPos(n_buff_sub)=-n_b_cells_sub-n_buff_sub

              NewS(s)=n_faces_sub+n_buff_sub
            end if
          end if  ! c2 > 0
        end do    ! through sides

        ! Faces on the "copy" boundary
        do s=1,n_copy
          c1=CopyS(1,s)  
          c2=CopyS(2,s) 
          if( (proces(c1) == sub).and.(proces(c2) == subo) ) then
            n_buff_sub = n_buff_sub+1
            NCFsub = NCFsub+1
            BuSeIn(n_buff_sub)=NewC(c1) ! Buffer Send Index 
            BuReIn(n_buff_sub)=c2 
            BufPos(n_buff_sub)= - (-n_b_cells_sub-n_buff_sub) ! watch the sign
          end if
          if( (proces(c2) == sub).and.(proces(c1) == subo) ) then
            n_buff_sub = n_buff_sub+1
            NCFsub = NCFsub+1
            BuSeIn(n_buff_sub)=NewC(c2) ! Buffer Send Index
            BuReIn(n_buff_sub)=c1 
            BufPos(n_buff_sub)= - (-n_b_cells_sub-n_buff_sub) ! watch the sign
          end if
        end do    ! through sides
        NBBe(subo)=n_buff_sub

        ! Write to buffer file
        write(9,'(A33)') '#===============================#' 
        write(9,'(A33)') '#   Conections with subdomain:  #' 
        write(9,'(A33)') '#===============================#' 
        write(9,'(I8)')  subo 
        write(9,'(A33)') '#--- Number of local connections:'
        write(9,'(I8)')  NBBe(subo)-NBBs(subo)+1 
        write(9,'(A40)') '#--- Local number in a buffer and index:'
        do b=NBBs(subo),NBBe(subo)
          write(9,'(2I8)') b-NBBs(subo)+1, BuSeIn(b) 
        end do
      end if 

    end do ! for subo

    call Save_Gmv_Grid(sub, n_nodes_sub, n_cells_sub)
    call Save_Cns_Geo(sub, n_cells_sub, n_faces_sub, n_b_cells_sub,  &
                      n_buff_sub,NCFsub)
    call Save_Gmv_Links(sub, n_nodes_sub, n_cells_sub, n_faces_sub,   &
                        n_b_cells_sub, n_buff_sub)

    write(*,*) 'Test:'
    write(*,*) 'n_nodes_sub   =', n_nodes_sub
    write(*,*) 'n_cells_sub   =', n_cells_sub
    write(*,*) 'n_faces_sub   =', n_faces_sub
    write(*,*) 'n_b_cells_sub =', n_b_cells_sub

    write(*,*) '====================================' 
    write(*,*) 'Subdomain   ', sub
    write(*,*) 'Buffer size ', n_buff_sub
    do subo=1,n_sub
      if(subo /= sub) then
        write(*,*) 'Connections with ', subo ,' : ',                &
          NBBe(subo)-NBBs(subo)+1,                                  &
          n_b_cells_sub+NBBs(subo),                                 &
          n_b_cells_sub+NBBe(subo) 
      end if 
    end do ! for subo
    write(*,*) '------------------------------------' 

  end do   ! through subdomains

  close(9)

  !--------------------------------------------------!
  !                                                  !
  !   Save the entire domain with renumbered cells   !
  !                                                  !
  !--------------------------------------------------!
  do n=1,NN
    NewN(n)=n
  end do
  do c=1,NC
    NewC(c)=0
  end do
  do s=1,NS
    NewS(s)=0
  end do

  n_cells_sub = 0     ! number of cells renumbered
  do sub=1,n_sub
    do c=1,NC
      if(proces(c) == sub) then
        n_cells_sub=n_cells_sub+1
        NewC(c)=n_cells_sub
      end if
    end do
  end do

  n_faces_sub = 0     ! number of sides renumbered
  do sub=1,n_sub
    do s=1,NS
      c1 = SideC(1,s)
      c2 = SideC(2,s)
      if(proces(c1) == sub) then
        n_faces_sub=n_faces_sub+1
        NewS(s)=n_faces_sub
      end if
    end do
  end do
  write(*,*) 'Number of sides: ', NS, n_faces_sub

  call Sort_Cells_By_Index(grid % cells, NewC(1), NC)
  call Sort_Int_By_Index(proces(1),  NewC(1),NC)
  call Sort_Int_By_Index(material(1),NewC(1),NC)

  call Sort_Int_By_Index(SideN(1,0), NewS(1),NS)
  call Sort_Int_By_Index(SideN(1,1), NewS(1),NS)
  call Sort_Int_By_Index(SideN(1,2), NewS(1),NS)
  call Sort_Int_By_Index(SideN(1,3), NewS(1),NS)
  call Sort_Int_By_Index(SideN(1,4), NewS(1),NS)
  call RNSort(Dx(1), NewS(1), NS)  ! this is important
  call RNSort(Dy(1), NewS(1), NS)  ! for plotting the
  call RNSort(Dz(1), NewS(1), NS)  ! grid with EpsPar()
  allocate(side_cell(NS,2))
  do s=1,NS
    side_cell(s,1) = SideC(1,s)
    side_cell(s,2) = SideC(2,s)
  end do
  call Sort_Int_By_Index(side_cell(1,1), NewS(1),NS)
  call Sort_Int_By_Index(side_cell(1,2), NewS(1),NS)
  do s=1,NS
    SideC(1,s) = side_cell(s,1)
    SideC(2,s) = side_cell(s,2)
  end do
  deallocate(side_cell)

  call Count_Materials

  call Save_Gmv_Grid(0, NN, NC)

  call Save_Cas(0, NN, NC, NS+NSsh)

  call Save_Eps_Decomposed

  end subroutine Number
