!======================================================================!
  subroutine ReadC(un, string, tn, ts, te) 
!----------------------------------------------------------------------!
!  Reads a line from a file (unit 9) and discards if it is comment.    !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use allp_mod, only: CMN_FILE          
  use all_mod,  only: cmn_line_count
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  character :: string*300
  integer   :: un, tn, ts(300), te(300)
!-------------------------------[Locals]-------------------------------!
  integer :: i
!======================================================================!
!  A comment is each line which begins with one of the following       !
!  characters: ! # $ % / & [ ] { } : @ < > * (Quite a lot, huh ?)      !
!  Input line must not exceed the lenght of 300 characters.            !
!  Note: not very carefully checked, but so far it works.              !
!----------------------------------------------------------------------!

  !-----------------------------------!
  !  Read the whole line into string  !
  !-----------------------------------!
1 read(un,'(A300)') string

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
  if( string  ==  '' ) goto 1 ! see: man ascii

  !----------------------!
  !  Skip comment lines  !
  !----------------------!
  if( string(1:1) == '!' .or.                                       &
      string(1:1) == '#' .or.                                       &
      string(1:1) == '$' .or.                                       &
      string(1:1) == '%' .or.                                       &
      string(1:1) == '/' .or.                                       &
      string(1:1) == '&' .or.                                       &
      string(1:1) == '[' .or.                                       &
      string(1:1) == ']' .or.                                       &
      string(1:1) == '{' .or.                                       &
      string(1:1) == '}' .or.                                       &
      string(1:1) == ':' .or.                                       &
      string(1:1) == ';' .or.                                       &
      string(1:1) == '@' .or.                                       &
      string(1:1) == '<' .or.                                       &
      string(1:1) == '>' .or.                                       &
      string(1:1) == '*' ) goto 1

  !--------------------------------------!
  !  Parse tokens. This is somehow cool  !
  !--------------------------------------!
  tn = 0
  if(string(1:1) >= '!') then
    tn = 1
    ts(1)=1
  end if
  do i=1,298
    if( string(i:i)  < '!' .and. string(i+1:i+1) >= '!') then
      tn = tn + 1
      ts(tn) = i+1
    end if
    if( string(i:i) >= '!' .and. string(i+1:i+1)  < '!') then
      te(tn) = i
    end if
  end do

  end subroutine ReadC
