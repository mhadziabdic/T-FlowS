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
  !  Read the whole line into whole  !
  !-----------------------------------!
1 read(un,'(A300)') line % whole

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
  if( line % whole  ==  '' ) goto 1 ! see: man ascii

  !----------------------!
  !  Skip comment lines  !
  !----------------------!
  if( trim(line % whole(1:1)) == '!' .or.               &
      trim(line % whole(1:1)) == '#' .or.               &
      trim(line % whole(1:1)) == '%' ) goto 1

  !--------------------------------------!
  !  Parse tokens. This is somehow cool  !
  !--------------------------------------!
  line % n_tokens = 0
  if(line % whole(1:1) >= '!') then
    line % n_tokens = 1
    line % s(1)=1
  end if
  do i=1,298
    if( line % whole(i:  i  ) <  '!' .and.  &
        line % whole(i+1:i+1) >= '!') then
      line % n_tokens = line % n_tokens + 1
      line % s(line % n_tokens) = i+1
    end if
    if( line % whole(i  :i  ) >= '!' .and.  &
        line % whole(i+1:i+1) <  '!') then
      line % e(line % n_tokens) = i
    end if
  end do

  ! Chop them up
  do i = 1, line % n_tokens
    line % tokens(i) = line % whole(line % s(i) : line % e(i))
  end do

  end subroutine
