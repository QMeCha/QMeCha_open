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
pure function regularization_f( d, eps ) result(f)
   use fortran_kinds_v, only: dp
   implicit none
   real(dp), intent(in) :: d, eps
   real(dp)             :: r, f
   if (eps.le.0.0_dp) then
      f=1.0_dp
   else
      r = abs(d/eps)
      if (r.lt.1.0_dp) then
         f=7.0_dp*r**6 - 15.0_dp * r**4 + 9.0_dp * r**2
      else
         f=1.0_dp
      endif
   endif
end function regularization_f
