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
subroutine sm_upd_inv(dim, d, i_n, vec_m, vec_Minv_n, vec_Mm, g, inv_M)
   use fortran_kinds_v, only: dp, int32
   implicit none
   character(1),                      intent(in)    :: dim
   integer(int32),                    intent(in)    :: d
   integer(int32),                    intent(in)    :: i_n
   real(dp), dimension(d),          intent(in)    :: vec_m
   real(dp), dimension(d),          intent(in)    :: vec_Minv_n
   real(dp), dimension(d),          intent(inout) :: vec_Mm
   real(dp),                          intent(in)    :: g
   real(dp), dimension(d,d),      intent(inout) :: inv_M
   external                                :: dgemv, dger
   if(dim.eq.'C') then
      call dgemv('N', d, d, 1.0_dp, inv_M, d, vec_m,&
      & 1_int32, 0.0_dp, vec_Mm, 1_int32)
      vec_Mm(i_n) = vec_Mm(i_n) - 1.0_dp
      call dger( d, d, -1.0_dp/g, vec_Mm, 1_int32, vec_Minv_n, 1_int32, &
      &  inv_M, d )
   elseif(dim.eq.'R') then
      call dgemv('T', d, d, 1.0_dp, inv_M, d, vec_m,&
      & 1_int32, 0.0_dp, vec_Mm, 1_int32)
      vec_Mm(i_n) = vec_Mm(i_n) - 1.0_dp
      call dger( d, d, -1.0_dp / g, vec_Minv_n, 1_int32, vec_Mm, 1_int32, &
      & inv_M, d )
   endif
end subroutine sm_upd_inv
subroutine sm_upd_inv_t( dir, n, vec_dm, vec_Mn, g, work, inv_M )
   use fortran_kinds_v, only: dp, int32
   implicit none
   character(1),                 intent(in)    :: dir
   integer(int32),               intent(in)    :: n
   real(dp), dimension(n),     intent(in)    :: vec_dm
   real(dp), dimension(n),     intent(in)    :: vec_Mn
   real(dp), dimension(n),     intent(inout) :: work
   real(dp),                     intent(in)    :: g
   real(dp), dimension(n,n), intent(inout) :: inv_M
   external                                :: dgemv, dger
   if(dir.eq.'C') then
      call dgemv('N', n, n, 1.0_dp, inv_M, n, vec_dm,&
      & 1_int32, 0.0_dp, work, 1_int32)
      call dger( n, n, -1.0_dp/g, work, 1_int32, vec_Mn, 1_int32, &
      &  inv_M, n )
   elseif(dir.eq.'R') then
      call dgemv('T', n, n, 1.0_dp, inv_M, n, vec_dm,&
      & 1_int32, 0.0_dp, work, 1_int32)
      call dger( n, n, - 1.0_dp / g, vec_Mn, 1_int32, work, 1_int32, &
      & inv_M, n )
   endif
end subroutine sm_upd_inv_t
