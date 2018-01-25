subroutine Append_Data_To_Array_Rank_2(data_ar_r2, data_ar_to_append_r2)
!------------------------------------------------------------------------------!
!   Appends one real*8 array to another                                        !
!------------------------------------------------------------------------------!
!   Requires explicit interface                                                !
!------------------------------------------------------------------------------!
  implicit none  
!---------------------------------[Arguments]----------------------------------!
  integer, allocatable :: data_ar_r2(:,:)
  integer, allocatable :: data_ar_to_append_r2(:,:)
  integer, allocatable :: tmp_buffer_int(:,:)
!-----------------------------------[Locals]-----------------------------------!
  integer              :: x1, x2, y1, y2, d1, d2
!------------------------------------------------------------------------------!

  y1 = lbound(data_ar_to_append_r2, dim=2) ! data_ar_to_append_r2
  y2 = ubound(data_ar_to_append_r2, dim=2) ! data_ar_to_append_r2
  d1 = lbound(data_ar_to_append_r2, dim=1) ! data_ar_r2
  d2 = lbound(data_ar_to_append_r2, dim=1) ! data_ar_r2

  if (.not.allocated(data_ar_r2)) then
    allocate(data_ar_r2(d1:d2, y1:y2))

    data_ar_r2(d1:d2, y1:y2) = data_ar_to_append_r2(d1:d2, y1:y2)
  else
    x1 = lbound(data_ar_r2, dim=2) ! data_ar_r2
    x2 = ubound(data_ar_r2, dim=2) ! data_ar_r2

    allocate(tmp_buffer_int(d1:d2, x1:x2))
    tmp_buffer_int(d1:d2, x1:x2) = data_ar_r2(d1:d2, x1:x2)

    deallocate(data_ar_r2)
    allocate(data_ar_r2(d1:d2, x1:x2 + y2 - y1 + 1))
    
    data_ar_r2(d1:d2, x1:x2) = tmp_buffer_int(d1:d2, x1:x2)
    deallocate(tmp_buffer_int)

    data_ar_r2(d1:d2, x2+1:x2 + y2 - y1 + 1) = data_ar_to_append_r2(d1:d2, y1:y2)
  end if

  deallocate(data_ar_to_append_r2)

end subroutine Append_Data_To_Array_Rank_2