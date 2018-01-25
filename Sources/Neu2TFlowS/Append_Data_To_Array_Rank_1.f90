subroutine Append_Data_To_Array_Rank_1(data_ar_r1, data_ar_to_append_r1)
!------------------------------------------------------------------------------!
!   Appends one real*8 array to another                                        !
!------------------------------------------------------------------------------!
!   Requires explicit interface                                                !
!------------------------------------------------------------------------------!
  implicit none  
!---------------------------------[Arguments]----------------------------------!
  real*8, allocatable  :: data_ar_r1(:)
  real*8, allocatable  :: data_ar_to_append_r1(:)
  real*8, allocatable  :: tmp_buffer_real(:)
!-----------------------------------[Locals]-----------------------------------!
  integer              :: x1, x2, y1, y2
!------------------------------------------------------------------------------!

  y1 = lbound(data_ar_to_append_r1, dim=1) ! data_ar_to_append_r1
  y2 = ubound(data_ar_to_append_r1, dim=1) ! data_ar_to_append_r1


  if (.not.allocated(data_ar_r1)) then
    allocate(data_ar_r1(y1:y2))

    data_ar_r1(y1:y2) = data_ar_to_append_r1(y1:y2)
  else
    x1 = lbound(data_ar_r1, dim=1) ! data_ar_r1
    x2 = ubound(data_ar_r1, dim=1) ! data_ar_r1

    allocate(tmp_buffer_real(x1:x2))
    tmp_buffer_real(x1:x2) = data_ar_r1(x1:x2)

    deallocate(data_ar_r1)
    allocate(data_ar_r1(x1:x2 + y2 - y1 + 1))
    
    data_ar_r1(x1:x2) = tmp_buffer_real(x1:x2)
    deallocate(tmp_buffer_real)

    data_ar_r1(x2+1:x2 + y2 - y1 + 1) = data_ar_to_append_r1(y1:y2)
  end if

  deallocate(data_ar_to_append_r1)

end subroutine Append_Data_To_Array_Rank_1