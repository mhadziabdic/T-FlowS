!==============================================================================!
  subroutine Work_Mod_Allocate_Integer_Faces(grid, n)
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: grid
  integer                 :: n    ! number of arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nf
!==============================================================================!

  ! Get number of cells and boundary cells
  nf = grid % n_faces

  ! Store the pointer to the grid
  pnt_grid => grid

  ! Allocate requested memory
  allocate(i_face_01(nf));   i_face_01 = 0;   if(n ==  1) return
  allocate(i_face_02(nf));   i_face_02 = 0;   if(n ==  2) return
  allocate(i_face_03(nf));   i_face_03 = 0;   if(n ==  3) return
  allocate(i_face_04(nf));   i_face_04 = 0;   if(n ==  4) return
  allocate(i_face_05(nf));   i_face_05 = 0;   if(n ==  5) return
  allocate(i_face_06(nf));   i_face_06 = 0;   if(n ==  6) return
  allocate(i_face_07(nf));   i_face_07 = 0;   if(n ==  7) return
  allocate(i_face_08(nf));   i_face_08 = 0;   if(n ==  8) return
  allocate(i_face_09(nf));   i_face_09 = 0;   if(n ==  9) return
  allocate(i_face_10(nf));   i_face_10 = 0;   if(n == 10) return
  allocate(i_face_11(nf));   i_face_11 = 0;   if(n == 11) return
  allocate(i_face_12(nf));   i_face_12 = 0;   if(n == 12) return
  allocate(i_face_13(nf));   i_face_13 = 0;   if(n == 13) return
  allocate(i_face_14(nf));   i_face_14 = 0;   if(n == 14) return
  allocate(i_face_15(nf));   i_face_15 = 0;   if(n == 15) return
  allocate(i_face_16(nf));   i_face_16 = 0;   if(n == 16) return
  allocate(i_face_17(nf));   i_face_17 = 0;   if(n == 17) return
  allocate(i_face_18(nf));   i_face_18 = 0;   if(n == 18) return
  allocate(i_face_19(nf));   i_face_19 = 0;   if(n == 19) return
  allocate(i_face_20(nf));   i_face_20 = 0;   if(n == 20) return
  allocate(i_face_21(nf));   i_face_21 = 0;   if(n == 21) return
  allocate(i_face_22(nf));   i_face_22 = 0;   if(n == 22) return
  allocate(i_face_23(nf));   i_face_23 = 0;   if(n == 23) return
  allocate(i_face_24(nf));   i_face_24 = 0;   if(n == 24) return
  allocate(i_face_25(nf));   i_face_25 = 0;   if(n == 25) return
  allocate(i_face_26(nf));   i_face_26 = 0;   if(n == 26) return
  allocate(i_face_27(nf));   i_face_27 = 0;   if(n == 27) return
  allocate(i_face_28(nf));   i_face_28 = 0;   if(n == 28) return
  allocate(i_face_29(nf));   i_face_29 = 0;   if(n == 29) return
  allocate(i_face_30(nf));   i_face_30 = 0;   if(n == 30) return

  end subroutine
