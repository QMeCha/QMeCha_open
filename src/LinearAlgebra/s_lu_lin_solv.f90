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
subroutine lu_lin_solv( n_tot_par, n_par, mat, eps, vec_M, f, da_vec)
   use fortran_kinds_v, only:dp, int32
   implicit none
   integer(int32),                           intent(in)    :: n_tot_par
   integer(int32),                           intent(in)    :: n_par
   real(dp), dimension(n_tot_par,n_tot_par),     intent(inout) :: mat
   real(dp),                                 intent(in)    :: eps
   real(dp), dimension(n_tot_par),             intent(inout) :: vec_M
   real(dp), dimension(n_tot_par),             intent(inout) :: f
   real(dp), dimension(n_tot_par),             intent(out)   :: da_vec
   integer(int32)                                :: ierr
   integer(int32)                                :: i1, i2
   da_vec(1:n_par) = f(1:n_par)
   if (n_par.lt.n_tot_par) then
      f(n_par+1:n_tot_par) = 0.0_dp
      da_vec(n_par+1:n_tot_par) = 0.0_dp
   endif
   do i1 = 1_int32, n_par
      vec_M(i1) = 1.0_dp / sqrt( mat(i1,i1) + eps )
      da_vec(i1) = da_vec(i1) * vec_M(i1)
   enddo
   do i1 = 1_int32, n_par
      do i2 = i1+1_int32, n_par
         mat(i2,i1) = mat(i2,i1) * vec_M(i1) * vec_M(i2)
         mat(i1,i2) = mat(i2,i1)
      enddo ; enddo
   do i1 = 1_int32, n_tot_par
      mat(i1,i1) = 1.0_dp
   enddo
   call dposv( 'L', n_tot_par, 1_int32, mat(1:n_tot_par,1:n_tot_par), n_tot_par, da_vec(1:n_tot_par), n_tot_par, ierr )
   do i1 = 1_int32, n_par
      da_vec(i1) = da_vec(i1) * vec_M(i1)
      mat(i1,i1) = 1.0_dp / ( vec_M(i1)**2 ) - eps
      do i2 = i1+1_int32, n_par
         mat(i2,i1) = mat(i2,i1) / ( vec_M(i1) * vec_M(i2) )
         mat(i1,i2) = mat(i2,i1)
      enddo ; enddo
   do i1 = n_par+1, n_tot_par
      mat(i1,i1) = 0.0_dp
   enddo
   if (n_par.lt.n_tot_par) da_vec(n_par+1:n_tot_par) = 0.0_dp
end subroutine lu_lin_solv
subroutine lu_lin_solv_ext( n_tot_par, n_par, mat, eps, vec_M, f, da_vec)
   use fortran_kinds_v, only:dp, int32
   implicit none
   integer(int32),                           intent(in)    :: n_tot_par
   integer(int32),                           intent(in)    :: n_par
   real(dp), dimension(n_tot_par,n_tot_par), intent(inout) :: mat
   real(dp),                                 intent(in)    :: eps
   real(dp), dimension(n_tot_par),           intent(inout) :: vec_M
   real(dp), dimension(n_tot_par),           intent(inout) :: f
   real(dp), dimension(n_tot_par),           intent(out)   :: da_vec
   integer(int32) :: ierr
   integer(int32) :: i1, i2
   da_vec(1:n_par) = f(1:n_par)
   if (n_par.lt.n_tot_par) then
      f(n_par+1:n_tot_par) = 0.0_dp
      da_vec(n_par+1:n_tot_par) = 0.0_dp
   endif
   do i1 = 1_int32, n_par
      vec_M(i1) = 1.0_dp / sqrt( mat(i1,i1) )
      da_vec(i1) = da_vec(i1) * vec_M(i1)
   enddo
   do i1 = 1_int32, n_par
      mat(i1,i1) = 1.0_dp + eps
      do i2 = i1+1_int32, n_par
         mat(i2,i1) = mat(i2,i1) * vec_M(i1) * vec_M(i2)
         mat(i1,i2) = mat(i2,i1)
      enddo ; enddo
   do i1 = n_par + 1_int32, n_tot_par
      mat(i1,i1) = 1.0_dp  + eps
   enddo
   call dposv( 'L', n_tot_par, 1_int32, mat(1:n_tot_par,1:n_tot_par), n_tot_par, da_vec(1:n_tot_par), n_tot_par, ierr )
   do i1 = 1_int32, n_par
      da_vec(i1) = da_vec(i1) * vec_M(i1)
      mat(i1,i1) = 1.0_dp /  vec_M(i1)**2
      do i2 = i1+1_int32, n_par
         mat(i2,i1) = mat(i2,i1) / ( vec_M(i1) * vec_M(i2) )
         mat(i1,i2) = mat(i2,i1)
      enddo ; enddo
   do i1 = n_par+1, n_tot_par
      mat(i1,i1) = 0.0_dp
   enddo
   if (n_par.lt.n_tot_par) da_vec(n_par+1:n_tot_par) = 0.0_dp
end subroutine lu_lin_solv_ext
