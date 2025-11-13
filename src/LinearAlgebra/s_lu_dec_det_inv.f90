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
subroutine lu_dec_det_inv( n, A, det_A, inv_A )
   use fortran_kinds_v, only: dp, int32
   implicit none
   integer(int32),                     intent(in)    :: n
   real(dp),       dimension(n,n), intent(in)    :: A
   real(dp),                           intent(inout) :: det_A
   real(dp),       dimension(n,n), intent(inout) :: inv_A
   real(dp),       dimension(n) :: work   ! work array for LAPACK
   integer(int32), dimension(n) :: ipiv   ! pivot indices
   integer(int32) :: info, i1
   external dgetrf, dgetri, dlacpy
   if (n.eq.1_int32) then
      det_A = A(1,1)
      inv_A(1,1) = 1.0_dp / A(1,1)
      return
   elseif ( n.eq.2_int32) then
      det_A = A(1,1) * A(2,2) - A(1,2) * A(2,1)
      inv_A(1,2) = -A(1,2) ; inv_A(2,1) = -A(2,1)
      inv_A(2,2) =  A(1,1) ; inv_A(1,1) =  A(2,2)
      inv_A = inv_A / det_A
      return
   endif
   inv_A = A
   call dgetrf(n, n, inv_A, n, ipiv, info)
   if (info /= 0) stop '(lu_dec_det_inv) Matrix is numerically singular'
   det_A = 1.0_dp
   do i1  = 1_int32, n
      if (ipiv(i1).eq.i1) then
         det_A =   det_A * inv_A(i1,i1)
      else
         det_A = - det_A * inv_A(i1,i1)
      endif
   enddo
   call dgetri(n, inv_A, n, ipiv, work, n, info)
   if (info /= 0) stop '(lu_dec_det_inv) Matrix inversion failed'
end subroutine lu_dec_det_inv
subroutine lu_dec_det_inv_t( n, work, ipiv, A, det_A, inv_A )
   use fortran_kinds_v, only: dp, int32
   implicit none
   integer(int32),                     intent(in)    :: n
   real(dp),       dimension(n),     intent(inout) :: work   ! work array for LAPACK
   integer(int32), dimension(n),     intent(inout) :: ipiv   ! pivot indices
   real(dp),       dimension(n,n), intent(in)    :: A
   real(dp),                           intent(inout) :: det_A
   real(dp),       dimension(n,n), intent(inout) :: inv_A
   integer(int32) :: info, i1
   external dgetrf, dgetri, dlacpy
   if (n.eq.1_int32) then
      det_A = A(1,1)
      inv_A(1,1) = 1.0_dp / A(1,1)
      return
   elseif ( n.eq.2_int32) then
      det_A = A(1,1) * A(2,2) - A(1,2) * A(2,1)
      inv_A(1,2) = -A(1,2) ; inv_A(2,1) = -A(2,1)
      inv_A(2,2) =  A(1,1) ; inv_A(1,1) =  A(2,2)
      inv_A = inv_A / det_A
      return
   endif
   inv_A = A
   call dgetrf(n, n, inv_A, n, ipiv, info)
   if (info /= 0) stop '(lu_dec_det_inv) Matrix is numerically singular'
   det_A = 1.0_dp
   do i1  = 1_int32, n
      if (ipiv(i1).eq.i1) then
         det_A =   det_A * inv_A(i1,i1)
      else
         det_A = - det_A * inv_A(i1,i1)
      endif
   enddo
   call dgetri(n, inv_A, n, ipiv, work, n, info)
   if (info /= 0) stop '(lu_dec_det_inv) Matrix inversion failed'
end subroutine lu_dec_det_inv_t
subroutine lu_dec_inv( n, work, ipiv, A, inv_A )
   use fortran_kinds_v, only: dp, int32
   implicit none
   integer(int32),                     intent(in)    :: n
   real(dp),       dimension(n),     intent(inout) :: work   ! work array for LAPACK
   integer(int32), dimension(n),     intent(inout) :: ipiv   ! pivot indices
   real(dp),       dimension(n,n), intent(inout) :: A
   real(dp),       dimension(n,n), intent(inout) :: inv_A
   integer(int32) :: info
   external dgetrf, dgetri
   ipiv = 0_int32; work = 0.0_dp; inv_A = A
   if (n.eq.1_int32) then
      inv_A(1,1) = 1.0_dp / inv_A(1,1)
      return
   endif
! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
   call dgetrf(n, n, inv_A, n, ipiv, info)
   if (info /= 0) stop '(lu_dec_inv) Matrix is numerically singular'
   ! DGETRI computes the inverse of a matrix using the LU factorization
   ! computed by DGETRF.
   call dgetri(n, inv_A, n, ipiv, work, n, info)
   if (info /= 0) stop '(lu_dec_inv) Matrix inversion failed'
end subroutine lu_dec_inv
