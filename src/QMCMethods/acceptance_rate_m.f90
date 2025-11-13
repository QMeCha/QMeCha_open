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
module acceptance_rate_m
   use fortran_kinds_v,       only: int32, dp, stdout
   use write_lines_m
   use openmp_mpi_m
   use quantum_monte_carlo_v, only: n_wlk_max, n_wlk, dt_e, dt_p, dt_d, var_dt, arate
   use molecular_system_v,    only: n_po, n_el
   use qdo_system_v,          only: n_qdo
   implicit none
   integer(int32), save, public, allocatable, dimension(:) :: n_acc_mv_pw_e
   integer(int32), save, public, allocatable, dimension(:) :: n_tot_mv_pw_e
   real(dp),       save, public                            :: n_tot_mv_pn_e
   real(dp),       save, public                            :: n_acc_mv_pn_e
   real(dp),       save, public                            :: n_tot_mv_e
   real(dp),       save, public                            :: n_acc_mv_e
   real(dp),       save, public                            :: acc_rate_e
   integer(int32), save, public, allocatable, dimension(:) :: n_acc_mv_pw_p
   integer(int32), save, public, allocatable, dimension(:) :: n_tot_mv_pw_p
   real(dp),       save, public                            :: n_tot_mv_pn_p
   real(dp),       save, public                            :: n_acc_mv_pn_p
   real(dp),       save, public                            :: n_tot_mv_p
   real(dp),       save, public                            :: n_acc_mv_p
   real(dp),       save, public                            :: acc_rate_p
   integer(int32), save, public, allocatable, dimension(:) :: n_acc_mv_pw_d
   integer(int32), save, public, allocatable, dimension(:) :: n_tot_mv_pw_d
   real(dp),       save, public                            :: n_tot_mv_pn_d
   real(dp),       save, public                            :: n_acc_mv_pn_d
   real(dp),       save, public                            :: n_tot_mv_d
   real(dp),       save, public                            :: n_acc_mv_d
   real(dp),       save, public                            :: acc_rate_d
   public  :: ini_acc_rates, upd_var_dt, cmp_acc_rates, write_acc_rates, rst_acc_rates
contains
   subroutine ini_acc_rates()
      if (n_el.gt.0_int32) then
         allocate( n_acc_mv_pw_e(1:n_wlk_max) ) ; n_acc_mv_pw_e = 0_int32
         allocate( n_tot_mv_pw_e(1:n_wlk_max) ) ; n_tot_mv_pw_e = 0_int32
      endif
      if ( n_po.gt.0_int32 ) then
         allocate( n_acc_mv_pw_p(1:n_wlk_max) ) ; n_acc_mv_pw_p = 0_int32
         allocate( n_tot_mv_pw_p(1:n_wlk_max) ) ; n_tot_mv_pw_p = 0_int32
      endif
      if (n_qdo.gt.0_int32) then
         allocate( n_acc_mv_pw_d(1:n_wlk_max) ) ; n_acc_mv_pw_d = 0_int32
         allocate( n_tot_mv_pw_d(1:n_wlk_max) ) ; n_tot_mv_pw_d = 0_int32
      endif
      n_tot_mv_e = 0.0_dp ; n_acc_mv_e = 0.0_dp
      n_tot_mv_p = 0.0_dp ; n_acc_mv_p = 0.0_dp
      n_tot_mv_d = 0.0_dp ; n_acc_mv_d = 0.0_dp
   end subroutine ini_acc_rates
   subroutine rst_acc_rates( )
      if(n_el.gt.0_int32) then
         n_tot_mv_pw_e = 0_int32 ; n_acc_mv_pw_e = 0_int32
      endif
      if(n_po.gt.0_int32) then
         n_tot_mv_pw_p = 0_int32 ; n_acc_mv_pw_p = 0_int32
      endif
      if (n_qdo.gt.0_int32) then
         n_tot_mv_pw_d = 0_int32 ; n_acc_mv_pw_d = 0_int32
      endif
      n_tot_mv_e = 0.0_dp ; n_acc_mv_e = 0.0_dp
      n_tot_mv_p = 0.0_dp ; n_acc_mv_p = 0.0_dp
      n_tot_mv_d = 0.0_dp ; n_acc_mv_d = 0.0_dp
   end subroutine rst_acc_rates
   subroutine cmp_acc_rates
      if(n_el.gt.0_int32) then
         n_tot_mv_pn_e = dble(sum(n_tot_mv_pw_e(1:n_wlk)))
         n_acc_mv_pn_e = dble(sum(n_acc_mv_pw_e(1:n_wlk)))
         n_tot_mv_pw_e = 0_int32
         n_acc_mv_pw_e = 0_int32
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_allreduce( MPI_IN_PLACE, n_acc_mv_pn_e, 1, MPI_DOUBLE_PRECISION,&
         & MPI_SUM, MPI_COMM_WORLD, mpierr )
         call mpi_allreduce( MPI_IN_PLACE, n_tot_mv_pn_e, 1, MPI_DOUBLE_PRECISION,&
         & MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
         n_tot_mv_e = n_tot_mv_e + n_tot_mv_pn_e
         n_acc_mv_e = n_acc_mv_e + n_acc_mv_pn_e
         acc_rate_e = n_acc_mv_e / n_tot_mv_e
      endif
      if(n_po.gt.0_int32) then
         n_tot_mv_pn_p = dble(sum(n_tot_mv_pw_p(1:n_wlk)))
         n_acc_mv_pn_p = dble(sum(n_acc_mv_pw_p(1:n_wlk)))
         n_tot_mv_pw_p = 0_int32
         n_acc_mv_pw_p = 0_int32
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_allreduce( MPI_IN_PLACE, n_acc_mv_pn_p, 1, MPI_DOUBLE_PRECISION,&
         & MPI_SUM, MPI_COMM_WORLD, mpierr )
         call mpi_allreduce( MPI_IN_PLACE, n_tot_mv_pn_p, 1, MPI_DOUBLE_PRECISION,&
         & MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
         n_tot_mv_p = n_tot_mv_p + n_tot_mv_pn_p
         n_acc_mv_p = n_acc_mv_p + n_acc_mv_pn_p
         acc_rate_p = n_acc_mv_p / n_tot_mv_p
      endif
      if (n_qdo.gt.0_int32) then
         n_tot_mv_pn_d = dble(sum(n_tot_mv_pw_d(1:n_wlk)))
         n_acc_mv_pn_d = dble(sum(n_acc_mv_pw_d(1:n_wlk)))
         n_tot_mv_pw_d = 0_int32
         n_acc_mv_pw_d = 0_int32
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_allreduce( MPI_IN_PLACE, n_acc_mv_pn_d, 1, MPI_DOUBLE_PRECISION,&
         & MPI_SUM, MPI_COMM_WORLD, mpierr )
         call mpi_allreduce( MPI_IN_PLACE, n_tot_mv_pn_d, 1, MPI_DOUBLE_PRECISION,&
         & MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
         n_tot_mv_d = n_tot_mv_d + n_tot_mv_pn_d
         n_acc_mv_d = n_acc_mv_d + n_acc_mv_pn_d
         acc_rate_d = n_acc_mv_d / n_tot_mv_d
      endif 
   end subroutine cmp_acc_rates
   subroutine write_acc_rates( )
      if(n_el.gt.0_int32) then
         call write_variable_line(stdout,0,mpi_rank,2,"Acceptance rate of Monte Carlo moves for elec.", acc_rate_e * 100.0_dp,units="%" )
      endif
      if(n_po.gt.0_int32) then
         call write_variable_line(stdout,0,mpi_rank,2,"Acceptance rate of Monte Carlo moves for posi.", acc_rate_p * 100.0_dp,units="%" )
      endif
      if (n_qdo.gt.0_int32) then
         call write_variable_line(stdout,0,mpi_rank,2,"Acceptance rate of Monte Carlo moves for drud.", acc_rate_d * 100.0_dp,units="%" )
      endif
   end subroutine write_acc_rates
   subroutine upd_var_dt( prnt )
      logical, optional :: prnt
      logical           :: prnt_in
      if (.not.present(prnt)) then
         prnt_in = .true.
      else
         prnt_in = prnt
      endif
      if ( var_dt ) then
         if(n_el.gt.0_int32) then
            dt_e = dt_e * (erf(10.0_dp * (acc_rate_e - arate) * arate )+1.0_dp)
            if (prnt_in) call write_variable_line(stdout,0,mpi_rank,2,"Changing Time step for next sampling for elec.", dt_e )
         endif
         if(n_po.gt.0_int32) then
            dt_p = dt_p * (erf(10.0_dp * (acc_rate_p - arate) * arate )+1.0_dp)
            if (prnt_in) call write_variable_line(stdout,0,mpi_rank,2,"Changing Time step for next sampling for posi.", dt_p )
         endif
         if (n_qdo.gt.0_int32) then
            dt_d = dt_d * (erf(10.0_dp * (acc_rate_d - arate) * arate )+1.0_dp)
            if (prnt_in) call write_variable_line(stdout,0,mpi_rank,2,"Changing Time step for next sampling for drud.", dt_d )
         endif
         call rst_acc_rates()
      endif 
   end subroutine upd_var_dt
end module acceptance_rate_m
