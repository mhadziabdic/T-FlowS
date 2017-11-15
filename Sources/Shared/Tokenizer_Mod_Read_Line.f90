!==============================================================================!
  subroutine Tokenizer_Mod_Read_Line(un) 
!------------------------------------------------------------------------------!
!  Reads a line from a file (unit 9) and discards if it is comment.            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use allp_mod, only: CMN_FILE
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: un  ! unit
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!
!   A comment is each line which begins with "!", "#", "$", "%" or "*".        !
!   Input line must not exceed the lenght of 300 characters.                   !
!------------------------------------------------------------------------------!

  !-----------------------------------!
  !  Read the whole line into string  !
  !-----------------------------------!
1 read(un,'(A300)') token % string

  !------------------------------------------!
  !  If you are reading from command file    !
  !  (T-FlowS.cmn), increase the line count  !
  !------------------------------------------!
  if( un == CMN_FILE ) then
    cmn_line_count = cmn_line_count + 1
  end if

  !--------------------!
  !  Skip empty lines  !
  !--------------------!
  if( token % string  ==  '' ) goto 1 ! see: man ascii

  !----------------------!
  !  Skip comment lines  !
  !----------------------!
  if( token % string(1:1) == '!' .or.               &
      token % string(1:1) == '#' .or.               &
      token % string(1:1) == '%' ) goto 1

  !--------------------------------------!
  !  Parse tokens. This is somehow cool  !
  !--------------------------------------!
  token % n = 0
  if(token % string(1:1) >= '!') then
    token % n = 1
    token % s(1)=1
  end if
  do i=1,298
    if( token % string(i:  i  ) <  '!' .and.  &
        token % string(i+1:i+1) >= '!') then
      token % n = token % n + 1
      token % s(token % n) = i+1
    end if
    if( token % string(i  :i  ) >= '!' .and.  &
        token % string(i+1:i+1) <  '!') then
      token % e(token % n) = i
    end if
  end do

  end subroutine
