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
function is_numeric(string)
   use fortran_kinds_v, only: dp, int32
   implicit none
   character(len=*), intent(in) :: string
   real(dp)                     :: x
   integer(int32)               :: e
   logical                      :: is_numeric
   read(string,*,iostat = e) x
   is_numeric = e == 0
end function is_numeric
