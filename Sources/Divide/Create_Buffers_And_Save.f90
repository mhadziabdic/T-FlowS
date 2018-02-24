!==============================================================================!
  subroutine Create_Buffers_And_Save(grid)
!------------------------------------------------------------------------------!
!   Number the cells in each subdomain for subsequent separate saving.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod
  use all_mod
  use gen_mod 
  use div_mod
  use Comm_Mod, only: nbb_s, nbb_e
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: b, c, n, s, c1, c2, sub, subo, ln
  integer              :: n_nodes_sub, n_cells_sub, n_faces_sub,  &
                          n_bnd_cells_sub, n_buff_sub, NCSsub, NCFsub
  character(len=80)    :: name_out
  integer,allocatable  :: side_cell(:,:)
!==============================================================================!
!   Each subdomain needs two buffers: a send buffer and a receive buffer.      !
!   A receive buffer will be stored as aditional boundary cells for each       !
!   subdomain. So each subdomain will have grid % n_bnd_cells physical         !
!   boundary faces and nbb_C-grid % n_bnd_cells buffer bounndary cells.         !
!   It is handy to do it that way, because most of the algorythms can remain   !
!   the same as they are now.  They won't even "know" that they use values     !
!   from other processors.  On the other hand, a sending buffer has to be      !
!   allocated in a new separate array called simply buffer(). An additional    !
!   array is needed to keep track of all the indexes. That one is called       !
!   buffind().  buffind() has stored cell numbers from it's own subdomain      !
!   so that later they can be copied with (well, something like that):         !
!                                                                              !
!   do i=1,BUFFSIZ                                                             !
!     buffer(i) = U(buffind(i))                                                !
!   end do                                                                     !
!------------------------------------------------------------------------------!

  allocate (nbb_s(0:n_sub))
  allocate (nbb_e(0:n_sub))

  !-------------------------------!
  !                               !
  !   Browse through subdomains   !
  !                               !
  !-------------------------------!
  do sub = 1, n_sub

    call Name_File(sub, name_out, '.buf')
    open(9, file=name_out)
    print *, '# Creating the file: ', trim(name_out)

    write(9,'(A20)') '#------------------#'
    write(9,'(A20)') '#                  #'
    write(9,'(A20)') '#  Buffer indexes  #'
    write(9,'(A20)') '#                  #'
    write(9,'(A20)') '#------------------#'

    ! Cells
    n_cells_sub = 0     ! number of cells in subdomain
    do c = 1, grid % n_cells
      NewC(c)=0
    end do
    do c = 1, grid % n_cells
      if(proces(c) == sub) then
        n_cells_sub=n_cells_sub+1
        NewC(c)=n_cells_sub
      end if
    end do

    ! Nodes
    n_nodes_sub = 0     ! number of cells in subdomain
    do n = 1, grid % n_nodes
      NewN(n)=0
    end do
    do c = 1, grid % n_cells
      if(proces(c) == sub) then
        do ln = 1, grid % cells_n_nodes(c)
          NewN(grid % cells_n(ln,c))=-1
        end do
      end if
    end do
    do n = 1, grid % n_nodes
      if(NewN(n) == -1) then
        n_nodes_sub=n_nodes_sub+1
        NewN(n)=n_nodes_sub
      end if
    end do

    ! Faces & real boundary cells
    n_faces_sub     = 0  ! number of sides in subdomain
    n_bnd_cells_sub = 0  ! number of real boundary cells in subdomain
    NCSsub = 0
    do s = 1, grid % n_faces
      NewS(s)=0
    end do
    do c=-grid % n_bnd_cells,-1
      NewC(c)=0
    end do
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)  
      c2 = grid % faces_c(2,s) 
      if(c2  > 0) then
        if( (proces(c1) == sub).and.(proces(c2) == sub) ) then
          n_faces_sub=n_faces_sub+1
          NewS(s)=n_faces_sub
        end if
      else ! c2 < 0
        if( proces(c1) == sub ) then
          n_faces_sub =n_faces_sub+1
          NewS(s)=n_faces_sub        ! new number for the side
          n_bnd_cells_sub=n_bnd_cells_sub+1
          NewC(c2)=-n_bnd_cells_sub  ! new number for the boundary cell
        end if
      end if 
    end do

    do s = 1, grid % n_copy
      c1=grid % bnd_cond % copy_s(1,s)
      c2=grid % bnd_cond % copy_s(2,s)
      if( (proces(c1) == sub).and.(proces(c2) == sub) ) then
        NCSsub=NCSsub+1
      end if
    end do

    print '(a,i5,a)', ' # Saving subdomain ', sub, ' with:'
    print '(a,i9,a)', ' # ', n_cells_sub, ' cells'
    print '(a,i9,a)', ' # ', n_nodes_sub, ' nodes' 
    print '(a,i9,a)', ' # ', n_faces_sub, ' sides' 
    print '(a,i9,a)', ' # ', n_bnd_cells_sub, ' physical boundary cells' 

    !--------------------!
    !   Create buffers   !
    !--------------------!
    n_buff_sub = 0
    NCFsub = 0
    write(9,'(A30)') '# Number of physical boundary cells:'
    write(9,'(I8)')  n_bnd_cells_sub   
    do subo=1,n_sub
      if(subo /= sub) then
        nbb_s(subo)=n_buff_sub+1

        ! Faces inside the domain
        do s = 1, grid % n_faces
          c1 = grid % faces_c(1,s)  
          c2 = grid % faces_c(2,s) 
          if(c2  > 0) then
            if( (proces(c1) == sub).and.(proces(c2) == subo) ) then
              n_buff_sub = n_buff_sub+1
              BuSeIn(n_buff_sub)=NewC(c1)  ! buffer send index 
              BuReIn(n_buff_sub)=c2        ! important for coordinate
              BufPos(n_buff_sub)=-n_bnd_cells_sub-n_buff_sub

              NewS(s)=n_faces_sub+n_buff_sub
            end if
            if( (proces(c2) == sub).and.(proces(c1) == subo) ) then
              n_buff_sub = n_buff_sub+1
              BuSeIn(n_buff_sub)=NewC(c2)  ! buffer send index
              BuReIn(n_buff_sub)=c1        ! important for coordinate
              BufPos(n_buff_sub)=-n_bnd_cells_sub-n_buff_sub

              NewS(s)=n_faces_sub+n_buff_sub
            end if
          end if  ! c2 > 0
        end do    ! through sides

        ! Faces on the "copy" boundary
        do s = 1, grid % n_copy
          c1=grid % bnd_cond % copy_s(1,s)  
          c2=grid % bnd_cond % copy_s(2,s) 
          if( (proces(c1) == sub).and.(proces(c2) == subo) ) then
            n_buff_sub = n_buff_sub+1
            NCFsub = NCFsub+1
            BuSeIn(n_buff_sub)=NewC(c1) ! buffer send index 
            BuReIn(n_buff_sub)=c2 
            BufPos(n_buff_sub)= - (-n_bnd_cells_sub-n_buff_sub) ! watch the sign
          end if
          if( (proces(c2) == sub).and.(proces(c1) == subo) ) then
            n_buff_sub = n_buff_sub+1
            NCFsub = NCFsub+1
            BuSeIn(n_buff_sub)=NewC(c2) ! buffer send index
            BuReIn(n_buff_sub)=c1 
            BufPos(n_buff_sub)= - (-n_bnd_cells_sub-n_buff_sub) ! watch the sign
          end if
        end do    ! through sides
        nbb_e(subo)=n_buff_sub

        ! Write to buffer file
        write(9,'(A33)') '#-------------------------------#' 
        write(9,'(A33)') '#   Conections with subdomain:  #' 
        write(9,'(A33)') '#-------------------------------#' 
        write(9,'(I8)')  subo 
        write(9,'(A30)') '# Number of local connections:'
        write(9,'(I8)')  nbb_e(subo)-nbb_s(subo)+1 
        write(9,'(A37)') '# Local number in a buffer and index:'
        do b=nbb_s(subo),nbb_e(subo)
          write(9,'(2I8)') b-nbb_s(subo)+1, BuSeIn(b) 
        end do
      end if 

    end do ! for subo

    call Save_Cns_Geo(grid,              &
                      sub,               &
                      n_nodes_sub,       &
                      n_cells_sub,       &
                      n_faces_sub,       &
                      n_bnd_cells_sub,   &
                      n_buff_sub,        &
                      NCFsub)

    call Save_Vtu_Cells(grid,         &
                        sub,          &
                        n_nodes_sub,  &
                        n_cells_sub)

    call Save_Vtu_Links(grid,             &
                        sub,              &
                        n_nodes_sub,      &
                        n_cells_sub,      &
                        n_faces_sub,      &
                        n_bnd_cells_sub,  &
                        n_buff_sub)

    call Save_Shadows(grid,         &
                      sub,          &
                      n_cells_sub)

    print *, '# Test:'
    print *, '# n_nodes_sub   =', n_nodes_sub
    print *, '# n_cells_sub   =', n_cells_sub
    print *, '# n_faces_sub   =', n_faces_sub
    print *, '# n_bnd_cells_sub =', n_bnd_cells_sub

    print *, '#=====================================' 
    print *, '# Subdomain   ', sub
    print *, '# Buffer size ', n_buff_sub
    do subo=1,n_sub
      if(subo /= sub) then
        print '(a,i9,a,3i9)', ' # Connections with ', subo ,' : ',  &
          nbb_e(subo)-nbb_s(subo)+1,                      &
          n_bnd_cells_sub+nbb_s(subo),                   &
          n_bnd_cells_sub+nbb_e(subo) 
      end if 
    end do ! for subo
    print *, '#-----------------------------------' 

  end do   ! through subdomains

  close(9)

  !--------------------------------------------------!
  !                                                  !
  !   Save the entire domain with renumbered cells   !
  !                                                  !
  !--------------------------------------------------!
  do n = 1, grid % n_nodes
    NewN(n)=n
  end do
  do c = 1, grid % n_cells
    NewC(c)=0
  end do
  do s = 1, grid % n_faces
    NewS(s)=0
  end do

  n_cells_sub = 0     ! number of cells renumbered
  do sub = 1, n_sub
    do c = 1, grid % n_cells
      if(proces(c) == sub) then
        n_cells_sub=n_cells_sub+1
        NewC(c)=n_cells_sub
      end if
    end do
  end do

  n_faces_sub = 0     ! number of sides renumbered
  do sub = 1, n_sub
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(proces(c1) == sub) then
        n_faces_sub=n_faces_sub+1
        NewS(s)=n_faces_sub
      end if
    end do
  end do
  print *, '# Number of faces: ', grid % n_faces, n_faces_sub

  ! It is not sorting nodes ... is it good?  I doubt
  call Grid_Mod_Sort_Cells_By_Index(grid, NewC(1), grid % n_cells)
  call Grid_Mod_Sort_Faces_By_Index(grid, NewS(1), grid % n_faces)

  call Sort_Int_By_Index(proces(1),  NewC(1), grid % n_cells)
  call Sort_Int_By_Index(material(1),NewC(1), grid % n_cells)

  ! This is important for plotting the grid with EpsPar()
  call Sort_Real_By_Index(grid % dx(1), NewS(1), grid % n_faces)
  call Sort_Real_By_Index(grid % dy(1), NewS(1), grid % n_faces)
  call Sort_Real_By_Index(grid % dz(1), NewS(1), grid % n_faces)

  allocate(side_cell(grid % n_faces, 2))
  do s = 1, grid % n_faces
    side_cell(s,1) = grid % faces_c(1,s)
    side_cell(s,2) = grid % faces_c(2,s)
  end do
  call Sort_Int_By_Index(side_cell(1,1), NewS(1), grid % n_faces)
  call Sort_Int_By_Index(side_cell(1,2), NewS(1), grid % n_faces)
  do s = 1, grid % n_faces
    grid % faces_c(1,s) = side_cell(s,1)
    grid % faces_c(2,s) = side_cell(s,2)
  end do
  deallocate(side_cell)

  call Save_Vtu_Cells(grid, 0, grid % n_nodes, grid % n_cells)
  call Save_Vtu_Faces(grid)

  end subroutine
