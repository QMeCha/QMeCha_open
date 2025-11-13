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
function lu_dec_det(n, A, A_LU, ipiv) result(det_A)
   use fortran_kinds_v, only: dp, int32
   implicit none
   integer(int32)                                :: n
   real(dp), dimension(:,:),       intent(in)    :: A
   real(dp), dimension(n,n),   intent(inout) :: A_LU
   integer(int32), dimension(n), intent(inout) :: ipiv
   real(dp)                                      :: det_A
   integer(int32)                                :: info
   integer(int32)                                :: i1
   external dgetrf, dlacpy
   if(n.eq.2) then
      det_A = A(1,1) * A(2,2) - A(1,2) * A (2,1)
   else
      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      call dlacpy('A',n, n, A, n, A_LU, n)
      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call dgetrf(n, n, A_LU, n, ipiv, info)
      if (info /= 0) stop '(lu_dec_det) Matrix is numerically singular'
      det_A = 1.0_dp
      do i1  = 1_int32, n
         if (ipiv(i1).eq.i1) then
            det_A =   det_A * A_LU(i1,i1)
         else
            det_A = - det_A * A_LU(i1,i1)
         endif
      enddo
      if (info /= 0) stop '(lu_dec_det) Computation of matrix determinant failed'
   endif
end function lu_dec_det
