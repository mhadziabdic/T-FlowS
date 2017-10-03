!======================================================================!
  subroutine ToUppr(string)
!----------------------------------------------------------------------!
!   Transforms string to uppercase.                                    !
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Parameters]-----------------------------!
  character*(*) :: string  
!-------------------------------[Locals]-------------------------------!
  integer :: i, value
!======================================================================!

  do i=1,len_trim(string)
    value = ichar(string(i:i))
    if (value  >=  97 .and. value  <=  122) then 
      string(i:i) = char(value-32) 
    end if
  end do

  end subroutine ToUppr
