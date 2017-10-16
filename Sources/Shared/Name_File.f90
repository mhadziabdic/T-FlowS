!======================================================================!
  subroutine Name_File(sub, namOut, ext, lext)
!----------------------------------------------------------------------!
!   Creates the file name depending on the subdomain and file type.    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use all_mod
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  integer       :: sub, lext
  character*(*) :: ext
  character*(*) :: namOut
!-------------------------------[Locals]-------------------------------!
  integer   :: c
  character :: numb*4
!======================================================================!

  namOut = name

  if(sub == 0) then
    namOut(len_trim(name)+1:len_trim(name)+1+lext-1) = ext(1:lext) 
  else
    write(numb,'(I4)') sub
    write(namOut(len_trim(name)+1:len_trim(name)+5),'(A5)') '-0000' 
    do c=1,4
      if( numb(c:c) >= '0' .and. numb(c:c) <= '9' )                 &
        namOut(len_trim(name)+1+c:len_trim(name)+1+c) = numb(c:c)
    end do
    namOut(len_trim(name)+6:len_trim(name)+6+lext-1) = ext(1:lext) 
  end if 

  end subroutine Name_File
