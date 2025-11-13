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
subroutine sm_ss_upd_inv(n, i_p, vec_P_new, vec_v_g, vec_v, v_n, inv_P)
   use fortran_kinds_v, only: dp, int32
   implicit none
   integer(int32),               intent(in)    :: n
   integer(int32),               intent(in)    :: i_p
   real(dp), dimension(n),     intent(in)    :: vec_P_new
   real(dp), dimension(n),     intent(in)    :: vec_v_g
   real(dp), dimension(n),     intent(inout) :: vec_v
   real(dp),                     intent(in)    :: v_n
   real(dp), dimension(n,n), intent(inout) :: inv_P
   external                                    :: dgemv, dger
   call dgemv( 'T', n, n, -1.0_dp, inv_P, n, vec_P_new, 1_int32, &
   & 0.0_dp, vec_v, 1_int32)
   vec_v(i_p) = vec_v(i_p) - 1.0_dp
   call dger( n, n, -1.0_dp/v_n, vec_v, 1_int32, vec_v_g, 1_int32, inv_P, n )
   call dger( n, n, +1.0_dp/v_n, vec_v_g, 1_int32, vec_v, 1_int32, inv_P, n )
end subroutine sm_ss_upd_inv
