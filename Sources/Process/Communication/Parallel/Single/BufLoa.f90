!======================================================================!
  subroutine BufLoa 
!----------------------------------------------------------------------!
! Reads: NAME.buf                                                      !
! ~~~~~~                                                               !
!------------------------------[Modules]-------------------------------!
  use all_mod
  use par_mod
!----------------------------------------------------------------------!
  implicit none
!-------------------------------[Locals]-------------------------------!
  integer           :: c, dummy 
  integer           :: sub, subo, NNsub,NCsub,NSsub,NBCsub,NBFsub
  character(len=80) :: name_in
!======================================================================!
!  Each subdomain needs two buffers: a send buffer and a receive buffer.
!  A receive buffer will be stored as aditional boundary cells for each
!  subdomain. So each subdomain will have NBC physical boundary faces
!  and NBBC-NBC buffer bounndary cells. It is handy to do it that way,
!  because most of the algorythms can remain the same as they are now.
!  They won't even "know" that they use values from other processors.
!  On the other hand, a sending buffer has to be allocated in a new 
!  separate array called simply buffer(). An additional array is needed 
!  to keep track of all the indexes. That one is called BufInd().
!  BufInd() has stored cell numbers from it's own subdomain so that
!  later they can be copied with (well, something like that):
!  do i=1,BUFFSIZ
!    buffer(i) = U(BufInd(i))
!  end do
!----------------------------------------------------------------------!

  if(n_proc == 0) return

  call Name_File(this_proc, name_in, '.buf')
  open(9, file=name_in)
  if(this_proc < 2) print *, '# Now reading the file:', name_in

  allocate (NBBs(0:n_proc))
  allocate (NBBe(0:n_proc))

!///// number of physical boundary cells
  call ReadC(9,inp,tn,ts,te)
  read(inp,*) NBCsub

!///// initialize 
  do sub=0,n_proc
    NBBs(sub) = -(NBCsub) 
    NBBe(sub) = -(NBCsub)
  end do

!///// fill the indexes and the buffers
  do sub=1,n_proc
    if(sub  /=  this_proc) then

!----- connections with subdomain          
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) subo 

!----- number of local connections with subdomain sub 
      call ReadC(9,inp,tn,ts,te)
      read(inp,*) NBBe(sub)

      NBBs(sub) = NBBe(sub-1) - 1  
      NBBe(sub) = NBBs(sub) - NBBe(sub) + 1

      do c=NBBs(sub),NBBe(sub),-1
        call ReadC(9,inp,tn,ts,te)
        read(inp,*) dummy, BufInd(c) 
      end do 
    else
      NBBs(sub) = NBBe(sub-1)-1  ! just to become "sloppy" 
      NBBe(sub) = NBBe(sub-1)    ! this_proc will be needed for next 
    end if
  end do   ! through subdomains

  close(9)

!///// correct the "sloppy" indexes
  do sub=1,n_proc
    if(NBBe(sub)  > NBBs(sub)) then  
      NBBs(sub) = -1 
      NBBe(sub) = 0 
    end if
  end do 

  call wait

!->>>  print *, 'PE',this_proc, '#===================#' 
!->>>  print *, 'PE',this_proc, '# Check connections #' 
!->>>  print *, 'PE',this_proc, '#-------------------#' 
!->>>  do sub=1,n_proc
!->>>    write(*,'(A2,I2,3I7)') 'PE',this_proc, sub, NBBs(sub), NBBe(sub)
!->>>  end do   ! through subdomains

  end subroutine
