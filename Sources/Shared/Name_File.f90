!==============================================================================!
  subroutine Name_File(sub, name_out, ext, lext)
!------------------------------------------------------------------------------!
!   Creates the file name depending on the subdomain and file type.            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use all_mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer       :: sub, lext
  character*(*) :: ext
  character*(*) :: name_out
!-----------------------------------[Locals]-----------------------------------!
  integer   :: c
  character :: numb*4
!==============================================================================!

  name_out = name

  if(sub == 0) then
    name_out(len_trim(name)+1:len_trim(name)+1+lext-1) = ext(1:lext) 
  else
    write(numb,'(I4)') sub
    write(name_out(len_trim(name)+1:len_trim(name)+5),'(A5)') '-0000' 
    do c=1,4
      if( numb(c:c) >= '0' .and. numb(c:c) <= '9' )                 &
        name_out(len_trim(name)+1+c:len_trim(name)+1+c) = numb(c:c)
    end do
    name_out(len_trim(name)+6:len_trim(name)+6+lext-1) = ext(1:lext) 
  end if 

  end subroutine Name_File
