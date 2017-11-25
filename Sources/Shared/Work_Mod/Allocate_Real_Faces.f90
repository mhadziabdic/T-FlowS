!==============================================================================!
  subroutine Work_Mod_Allocate_Real_Faces(grid, n)
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
  allocate(r_face_01(nf));   r_face_01 = 0.0;   if(n ==  1) return
  allocate(r_face_02(nf));   r_face_02 = 0.0;   if(n ==  2) return
  allocate(r_face_03(nf));   r_face_03 = 0.0;   if(n ==  3) return
  allocate(r_face_04(nf));   r_face_04 = 0.0;   if(n ==  4) return
  allocate(r_face_05(nf));   r_face_05 = 0.0;   if(n ==  5) return
  allocate(r_face_06(nf));   r_face_06 = 0.0;   if(n ==  6) return
  allocate(r_face_07(nf));   r_face_07 = 0.0;   if(n ==  7) return
  allocate(r_face_08(nf));   r_face_08 = 0.0;   if(n ==  8) return
  allocate(r_face_09(nf));   r_face_09 = 0.0;   if(n ==  9) return
  allocate(r_face_10(nf));   r_face_10 = 0.0;   if(n == 10) return
  allocate(r_face_11(nf));   r_face_11 = 0.0;   if(n == 11) return
  allocate(r_face_12(nf));   r_face_12 = 0.0;   if(n == 12) return
  allocate(r_face_13(nf));   r_face_13 = 0.0;   if(n == 13) return
  allocate(r_face_14(nf));   r_face_14 = 0.0;   if(n == 14) return
  allocate(r_face_15(nf));   r_face_15 = 0.0;   if(n == 15) return
  allocate(r_face_16(nf));   r_face_16 = 0.0;   if(n == 16) return
  allocate(r_face_17(nf));   r_face_17 = 0.0;   if(n == 17) return
  allocate(r_face_18(nf));   r_face_18 = 0.0;   if(n == 18) return
  allocate(r_face_19(nf));   r_face_19 = 0.0;   if(n == 19) return
  allocate(r_face_20(nf));   r_face_20 = 0.0;   if(n == 20) return
  allocate(r_face_21(nf));   r_face_21 = 0.0;   if(n == 21) return
  allocate(r_face_22(nf));   r_face_22 = 0.0;   if(n == 22) return
  allocate(r_face_23(nf));   r_face_23 = 0.0;   if(n == 23) return
  allocate(r_face_24(nf));   r_face_24 = 0.0;   if(n == 24) return
  allocate(r_face_25(nf));   r_face_25 = 0.0;   if(n == 25) return
  allocate(r_face_26(nf));   r_face_26 = 0.0;   if(n == 26) return
  allocate(r_face_27(nf));   r_face_27 = 0.0;   if(n == 27) return
  allocate(r_face_28(nf));   r_face_28 = 0.0;   if(n == 28) return
  allocate(r_face_29(nf));   r_face_29 = 0.0;   if(n == 29) return
  allocate(r_face_30(nf));   r_face_30 = 0.0;   if(n == 30) return

  end subroutine
