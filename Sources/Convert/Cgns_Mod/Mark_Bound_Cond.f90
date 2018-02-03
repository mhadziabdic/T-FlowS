!==============================================================================!
  subroutine Cgns_Mod_Mark_Bound_Cond
!------------------------------------------------------------------------------!
!   Mark nodes in x,y,z arrays from tria_connections & quad_connections        !
!   with sect_id                                                               !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

!  ! Color quad nodes with b.c.
!  if ( bc_id(sect_id) .ne. 0) then
!    do c = 1, n_nodes
!      do i = 1, last_quad
!        do j = 1, 4
!          if (c .eq. quad_connections(j,i)) then
!            bc_mark(c) = sect_id
!          end if
!        end do
!      end do
!    end do
!  end if

!  ! Color tria nodes with b.c.
!  if ( bc_id(sect_id) .ne. 0) then
!    do c = 1, n_nodes
!      do i = 1, last_tria
!        do j = 1, 3
!          if (c .eq. quad_connections(j,i)) then
!            bc_mark(c) = sect_id
!          end if
!        end do
!      end do
!    end do
!  end if

  end subroutine
