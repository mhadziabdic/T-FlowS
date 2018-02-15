!==============================================================================!
  subroutine Control_Mod_Read_Real_Item(keyword, def, val, verbose)
!------------------------------------------------------------------------------!
!   Working horse function to read integer value (argument "val") behind a     !
!   keyword (argument "keyword") in control file.  If not found, a default     !
!   vaue specified in argument "def" is used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: keyword
  real              :: def      ! default value
  real              :: val      ! spefified value, if found
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: found, reached_end
!==============================================================================!

  rewind(CONTROL_FILE_UNIT)

  ! Set default number of time steps
  val = def

  ! Browse through command file to see if user specificed it
  found = .false.
  do
    call Tokenizer_Mod_Read_Line(CONTROL_FILE_UNIT, reached_end)
    if(reached_end) goto 1
    if(line % tokens(1) == trim(keyword)) then
      read(line % tokens(2), *) val
      found = .true.
    end if
  end do

1 if( .not. found) then
    if(present(verbose)) then
      if(verbose) then
        print '(a,a,a)', '# Couldn''t find keyword: ', keyword, '.'
        print '(a,1pe12.5)', '# Using the default value of ', def
      end if
    end if
  end if 

  end subroutine