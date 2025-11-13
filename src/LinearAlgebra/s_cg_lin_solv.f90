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
subroutine cg_lin_solv( cg_acc_ini, n_tot_par, n_par, mat, eps, vec_M, f, da_vec)
   use fortran_kinds_v, only: dp, int32
   use openmp_mpi_m
   use write_lines_m
   use timings_m,       only: time_compute, time_compute_difference, t_lnalg_step, t_lnalg_tot
   implicit none
   real(dp),                                         intent(inout) :: cg_acc_ini
   integer(int32),                               intent(in)    :: n_tot_par
   integer(int32),                               intent(in)    :: n_par
   real(dp), dimension(n_tot_par,n_tot_par), intent(inout) :: mat
   real(dp),                                     intent(in)    :: eps
   real(dp), dimension(n_tot_par),             intent(inout) :: vec_M
   real(dp), dimension(n_tot_par),             intent(inout) :: f
   real(dp), dimension(n_tot_par),             intent(out)   :: da_vec
   integer(int32)                                :: ierr
   integer(int32)                                :: i1, i2
   real(dp), dimension(n_tot_par) :: r, p, Ap, z
   real(dp) :: r2, alpha, beta, rz, cg_acc
   call time_compute(t_lnalg_step)
   cg_acc = minval(abs(f(1:n_par)))*0.1_dp
   if (cg_acc.gt.cg_acc_ini) then
      cg_acc = cg_acc_ini
   endif
   call write_variable_line(stdout,0,mpi_rank,2,"CG solver: accuracy set to",cg_acc)
   do i1 = 1_int32, n_par
      vec_M(i1) = 1.0_dp / sqrt( mat(i1,i1) + eps )
      f(i1) = f(i1) * vec_M(i1)
   enddo
   if (n_par.lt.n_tot_par) f(n_par+1:n_tot_par) = 0.0_dp
   do i1 = 1_int32, n_tot_par
      mat(i1,i1) = 1.0_dp
      do i2 = i1+1_int32, n_tot_par
         mat(i2,i1) = mat(i2,i1) * vec_M(i1) * vec_M(i2)
         mat(i1,i2) = mat(i2,i1)
      enddo ; enddo
   r(:) = 0.0_dp
   p(:) = 0.0_dp
   da_vec(:) = f(:)
   r(:) = f(:)
   call dsymv ( 'L', n_par, -1.0_dp, mat, n_tot_par, da_vec, 1_int32, 1.0_dp, r, 1_int32 )
   p(:) = r(:)
   r2 = dot_product( r(1:n_par), r(1:n_par) )
   i1 = 0_int32
   do
      call dsymv ( 'L', n_tot_par, 1.0_dp, mat, n_tot_par, p, 1_int32, 0.0_dp, Ap, 1_int32 )
      alpha = r2 / dot_product ( p(1:n_par), Ap(1:n_par))
      da_vec(:) = da_vec(:) + alpha * p(1:n_par)
      beta = 1.0_dp / r2
      r(1:n_par) = r(1:n_par) - alpha * Ap(1:n_par)
      r2 = dot_product( r(1:n_par), r(1:n_par) )
      if ( sqrt (r2) .lt. cg_acc ) exit
      beta = beta * r2
      p(1:n_par) = r(1:n_par) + beta * p(1:n_par)
      i1 = i1 + 1_int32
   end do
   call write_variable_line(stdout,0,mpi_rank,2,"CG solver: number of steps to convergence",i1)
   call write_variable_line(stdout,0,mpi_rank,2,"CG solver: residual error",sqrt(r2))
   do i1 = 1_int32, n_tot_par
      da_vec(i1) = da_vec(i1) * vec_M(i1)
      mat(i1,i1) = 1.0_dp / ( vec_M(i1)**2 ) - eps
      do i2 = i1+1_int32, n_tot_par
         mat(i2,i1) = mat(i2,i1) / ( vec_M(i1) * vec_M(i2) )
         mat(i1,i2) = mat(i2,i1)
      enddo ; enddo
   if (n_par.lt.n_tot_par) da_vec(n_par+1:n_tot_par) = 0.0_dp
   call time_compute_difference(t_lnalg_step,time_cpu_tot=t_lnalg_tot)
   call write_variable_line(stdout,0,mpi_rank,2,"CG solver: execution time",t_lnalg_step,units='sec.' )
end subroutine cg_lin_solv
subroutine cg_lin_solv_mpi(cg_acc_ini, n_tot_par, n_par, n_bin, n_wlk, A, eps, vec_M, f, da_vec)
   use fortran_kinds_v, only: dp, int32, stdout
   use openmp_mpi_m
   use write_lines_m
   use timings_m,       only: time_compute, time_compute_difference, t_lnalg_step, t_lnalg_tot
   implicit none
   real(dp),                                         intent(inout) :: cg_acc_ini
   integer(int32),                                   intent(in)    :: n_tot_par
   integer(int32),                                   intent(in)    :: n_par
   integer(int32),                                   intent(in)    :: n_bin
   integer(int32),                                   intent(in)    :: n_wlk
   real(dp), dimension(n_tot_par,n_bin,n_wlk), intent(in)    :: A
   real(dp),                                         intent(in)    :: eps
   real(dp), dimension(n_tot_par),                 intent(inout) :: vec_M
   real(dp), dimension(n_tot_par),                 intent(inout) :: f
   real(dp), dimension(n_tot_par),                 intent(inout) :: da_vec
   real(dp), dimension(n_bin,n_wlk) :: OTx_pw
   real(dp), dimension(n_tot_par,n_wlk) :: r_pw
   real(dp), dimension(n_par) :: r, Ap, z
   real(dp), dimension(n_tot_par) :: p
   real(dp) :: alpha, beta, r2, rz, rs, cg_acc 
   integer(int32)  :: i1, iw
   call time_compute(t_lnalg_step)
   da_vec(1:n_par) = f(1:n_par)
   if (n_tot_par.gt.n_par) then
      da_vec(n_par+1:n_tot_par) = 0.0_dp
      f(n_par+1:n_tot_par) = 0.0_dp
   endif
   cg_acc = minval(abs(f(1:n_par)))*0.1_dp
   if (cg_acc.gt.cg_acc_ini) then
      cg_acc = cg_acc_ini
   endif
   call write_variable_line(stdout,0,mpi_rank,2,"CG solver: accuracy set to",cg_acc)
   r(:) = 0.0_dp
   p(:) = 0.0_dp
!$omp parallel default(shared) private(iw) reduction(+:r)
!$omp do
   do iw = 1_int32, n_wlk
      call dgemv( 'T', n_tot_par, n_bin, 1.0_dp, A(1:n_tot_par,:,iw), n_tot_par, da_vec, 1_int32, 0.0_dp, OTx_pw(:,iw), 1_int32)
      call dgemv( 'N', n_tot_par, n_bin, 1.0_dp , A(1:n_tot_par,:,iw), n_tot_par, OTx_pw(:,iw), 1_int32,   0.0_dp, r_pw(:,iw), 1_int32)
      r(:) = r(:) - r_pw(1:n_par,iw)
   enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
   call mpi_allreduce(MPI_IN_PLACE, r(:) ,     n_par, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
   vec_M(1:n_par) = 1.0_dp / ( vec_M(1:n_par) + eps )
   r(1:n_par) = r(1:n_par) + f(1:n_par) - eps * da_vec(1:n_par)
   z(1:n_par) = vec_M(1:n_par) * r(1:n_par)
   p(1:n_par) = z(1:n_par)
   r2 = dot_product( r(1:n_par), r(1:n_par) )
   rz = dot_product( r(1:n_par), z(1:n_par) )
   i1 = 0_int32
   do
      Ap(:) = 0.0_dp
!$omp parallel default(shared) private(iw) reduction(+:Ap)
!$omp do
      do iw = 1_int32, n_wlk
         call dgemv( 'T', n_tot_par, n_bin, 1.0_dp , A(1:n_tot_par,:,iw), n_tot_par, p, 1_int32, 0.0_dp, OTx_pw(:,iw), 1_int32)
         call dgemv( 'N', n_tot_par, n_bin, 1.0_dp , A(1:n_tot_par,:,iw), n_tot_par, OTx_pw(:,iw), 1_int32,   0.0_dp, r_pw(:,iw), 1_int32)
         Ap(1:n_par) = Ap(1:n_par) + r_pw(1:n_par,iw)
      enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE, Ap(1:n_par),  n_par, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      Ap(1:n_par) = Ap(1:n_par) + eps * p(1:n_par)
      alpha = rz / dot_product ( p(1:n_par), Ap)
      beta = 1.0_dp / rz
      da_vec(1:n_par) = da_vec(1:n_par) + alpha * p(1:n_par)
      r(1:n_par) = r(1:n_par) - alpha * Ap(1:n_par)
      r2 = dot_product( r(1:n_par), r(1:n_par) )
      rs=sqrt(r2)
      if ( (rs.lt. cg_acc)) exit
      z(1:n_par) = r(1:n_par) * vec_M(1:n_par)
      rz = dot_product( r(1:n_par), z(1:n_par) )
      p(1:n_par) = z(1:n_par) + rz *  beta * p(1:n_par)
      i1 = i1 + 1_int32
   enddo
   vec_M(1:n_par) = 1.0_dp / vec_M(1:n_par) - eps
   call write_variable_line(stdout,0,mpi_rank,2,"CG solver: number of steps to convergence",i1)
   call write_variable_line(stdout,0,mpi_rank,2,"CG solver: residual error",rs)
   call time_compute_difference(t_lnalg_step,time_cpu_tot=t_lnalg_tot)
   call write_variable_line(stdout,0,mpi_rank,2,"CG solver: execution time",t_lnalg_step,units='sec.' )
end subroutine cg_lin_solv_mpi
