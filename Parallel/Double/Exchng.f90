!======================================================================!
  subroutine Exchng(PHI) 
!----------------------------------------------------------------------!
!   Exchanges the values between the processors.                       !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
  use par_mod
  use pro_mod
!----------------------------------------------------------------------!
  implicit none
!------------------------------[Include]-------------------------------!
  include 'mpif.h'
!-----------------------------[Parameters]-----------------------------!
  real    :: PHI(-NbC:NC)
!-------------------------------[Locals]-------------------------------!
  integer :: c1, c2, sub, rtag, stag, length, error
  integer :: status(MPI_STATUS_SIZE)
!======================================================================!

!///// fill the buffers with new values
  do sub=1,NPro
    if( NBBe(sub)  <=  NBBs(sub) ) then  
      do c2=NBBs(sub),NBBe(sub),-1
        c1 = BufInd(c2)
        PHI(c2) = PHI(c1)
      end do
    end if
  end do   

!///// exchange the values
  do sub=1,NPro
    if( NBBe(sub)  <=  NBBs(sub) ) then  

      length = NBBs(sub) - NBBe(sub) + 1
      stag = MAXPRO*this + sub    ! tag for sending
      rtag = MAXPRO*sub + this    ! tag for receivinging

!===============================================
      call MPI_SENDRECV_REPLACE & 
!-------------------------------------+---------
             (PHI(NBBe(sub)),   & ! buffer  
              length,           & ! length   
              MPI_DOUBLE_PRECISION,       & ! datatype  
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
