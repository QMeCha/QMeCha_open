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
subroutine norm_sqr_mat(dir, dim, A, norm_A )
   use fortran_kinds_v, only: int32, dp
   implicit none
   character(1),                     intent(in)    :: dir
   integer(int32),                   intent(in)    :: dim
   real(dp), dimension(dim,dim), intent(in)    :: A
   real(dp),                         intent(inout) :: norm_A
   real(dp)       :: sum_dir
   integer(int32) :: i1
   norm_A = 1.0_dp
   if(dir.eq.'C') then
      do i1 = 1_int32, dim
         sum_dir = sum(A(1:dim,i1)**2)
         if ( sum_dir.lt.norm_A) norm_A = sum_dir
      enddo
   elseif(dir.eq.'R') then
      do i1 = 1_int32, dim
         sum_dir = sum(A(i1,1:dim)**2)
         if ( sum_dir.lt.norm_A) norm_A = sum_dir
      enddo
   elseif(dir.eq.'B') then
      do i1 = 1_int32, dim
         sum_dir = sum(A(1:dim,i1)**2)
         if ( sum_dir.lt.norm_A) norm_A = sum_dir
         sum_dir = sum(A(i1,1:dim)**2)
         if ( sum_dir.lt.norm_A) norm_A = sum_dir
      enddo
   else
      print *, "Error in norm_sqr_mat function"
   endif
end subroutine norm_sqr_mat
