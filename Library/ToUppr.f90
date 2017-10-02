!======================================================================!
  SUBROUTINE ToUppr(string)
!----------------------------------------------------------------------!
!   Transforms string to uppercase.                                    !
!----------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------[Parameters]-----------------------------!
  CHARACTER*(*) :: string  
!-------------------------------[Locals]-------------------------------!
  INTEGER :: i, value
!======================================================================!

  do i=1,len_trim(string)
    value = ichar(string(i:i))
    if (value  >=  97 .and. value  <=  122) then 
      string(i:i) = char(value-32) 
    end if
  end do

  END SUBROUTINE ToUppr
