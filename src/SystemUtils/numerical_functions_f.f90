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
pure function erf_1(x) result(f)
   use fortran_kinds_v, only: dp
   implicit none
   real(dp), intent(in) :: x
   real(dp), parameter :: a1 = 0.278393_dp
   real(dp), parameter :: a2 = 0.230389_dp
   real(dp), parameter :: a3 = 0.000972_dp
   real(dp), parameter :: a4 = 0.078108_dp
   real(dp) :: x_sgn, x_abs
   real(dp) :: f
   x_sgn = sign(1.0_dp,x)
   x_abs = abs(x)
   if (x_abs.le.1.0d-10) then
      f=0.0_dp
   else if (x_abs.gt.10.0_dp) then
      f=1.0_dp
   else
      f = 1.0_dp + a1*x_abs + a2*x_abs**2+a3*x_abs**3+a4*x_abs**4
      f = 1.0_dp - 1.0_dp / f**4
   endif
   f= x_sgn * f
end function erf_1
