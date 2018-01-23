! * ------------------------------------------------------------------------- *
! * CGNS - CFD General Notation System (http://www.cgns.org)                  *
! * CGNS/MLL - Mid-Level Library header file                                  *
! * Please see cgnsconfig.h file for this local installation configuration    *
! * ------------------------------------------------------------------------- *
!
! * ------------------------------------------------------------------------- *
!
!  This software is provided 'as-is', without any express or implied warranty.
!  In no event will the authors be held liable for any damages arising from
!  the use of this software.
!
!  Permission is granted to anyone to use this software for any purpose,
!  including commercial applications, and to alter it and redistribute it
!  freely, subject to the following restrictions:
!
!  1. The origin of this software must not be misrepresented; you must not
!     claim that you wrote the original software. If you use this software
!     in a product, an acknowledgment in the product documentation would be
!     appreciated but is not required.
!
!  2. Altered source versions must be plainly marked as such, and must not
!     be misrepresented as being the original software.
!
!  3. This notice may not be removed or altered from any source distribution.
!
! * ------------------------------------------------------------------------- *
!
!
!     file open modes
!
      integer, parameter :: CGIO_MODE_READ   = 0
      integer, parameter :: CGIO_MODE_WRITE  = 1
      integer, parameter :: CGIO_MODE_MODIFY = 2
!
!     database file types
!
      integer, parameter :: CGIO_FILE_NONE = 0
      integer, parameter :: CGIO_FILE_ADF  = 1
      integer, parameter :: CGIO_FILE_HDF5 = 2
      integer, parameter :: CGIO_FILE_ADF2 = 3
!
!     dimension limits
!
      integer, parameter :: CGIO_MAX_DATATYPE_LENGTH = 2
      integer, parameter :: CGIO_MAX_DIMENSIONS      = 12
      integer, parameter :: CGIO_MAX_NAME_LENGTH     = 32
      integer, parameter :: CGIO_MAX_LABEL_LENGTH    = 32
      integer, parameter :: CGIO_MAX_VERSION_LENGTH  = 32
      integer, parameter :: CGIO_MAX_ERROR_LENGTH    = 80
      integer, parameter :: CGIO_MAX_LINK_DEPTH      = 100
      integer, parameter :: CGIO_MAX_FILE_LENGTH     = 1024
      integer, parameter :: CGIO_MAX_LINK_LENGTH     = 4096

