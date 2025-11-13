!> QMeCha Open release
!>
!> All rights reserved @2025 Matteo Barborini
!> Released under CC-BY-NC-ND license
!>
!> @Author       : Matteo Barborini
!> @Version      : 0.1
!> @Release date : 24.10.2025
!> @Repository   : github.com/QMeCha/QMeCha_open
!>
module fortran_kinds_v
   implicit none
   intrinsic selected_int_kind, selected_real_kind
   integer, public, parameter :: sp    = selected_real_kind(6, 37)
   integer, public, parameter :: dp    = selected_real_kind(15, 307)
   integer, public, parameter :: qp    = selected_real_kind(33, 4931)
   integer, public, parameter :: int8  = selected_int_kind(2)
   integer, public, parameter :: int16 = selected_int_kind(4)
   integer, public, parameter :: int32 = selected_int_kind(9)
   integer, public, parameter :: int64 = selected_int_kind(18)
   integer, public, parameter :: ascii = kind('x')
   integer, public, parameter :: stdin  = 5
   integer, public, parameter :: stdout = 6
   integer, public, parameter :: stderr = 0
end module fortran_kinds_v
