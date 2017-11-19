!======================================================================!
  subroutine Exchange(grid, phi) 
!----------------------------------------------------------------------!
!   Exchanges the values between the processors.                       !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use par_mod
  use pro_mod
  use Grid_Mod
!----------------------------------------------------------------------!
  implicit none
!------------------------------[Include]-------------------------------!
  include 'mpif.h'
!-----------------------------[Arguments]------------------------------!
  type(Grid_Type) :: grid
  real            :: phi(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  integer :: c1, c2, sub, rtag, stag, length, error
  integer :: status(MPI_STATUS_SIZE)
!======================================================================!

!///// fill the buffers with new values
  do sub=1,n_proc
    if( NBBe(sub)  <=  NBBs(sub) ) then  
      do c2=NBBs(sub),NBBe(sub),-1
        c1 = BufInd(c2)
        phi(c2) = phi(c1)
      end do
    end if
  end do   

!///// exchange the values
  do sub=1,n_proc
    if( NBBe(sub)  <=  NBBs(sub) ) then  

      length = NBBs(sub) - NBBe(sub) + 1
      stag = (n_proc)*this_proc + sub    ! tag for sending
      rtag = (n_proc)*sub + this_proc    ! tag for receivinging

!===============================================
      call MPI_SENDRECV_REPLACE & 
!-------------------------------------+---------
             (phi(NBBe(sub)),   & ! buffer  
              length,           & ! length   
              MPI_REAL,         & ! datatype  
!-------------------------------------+---------
              (sub-1),          & ! dest,      
              stag,             & ! sendtag,    
!-------------------------------------+---------
              (sub-1),          & ! source,      
              rtag,             & ! recvtag,      
!-------------------------------------+---------
              MPI_COMM_WORLD,   &
              status,           &
              error) 
!===============================================

    end if  !  NBBe(sub)  /=  NBBs(sub)
  end do

  end subroutine Exchng
