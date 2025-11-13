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
module wavefunction_optimization_v
   use fortran_kinds_v, only: int32, dp
   use openmp_mpi_m,    only: n_mpi_tasks
   use quantum_monte_carlo_v, only: n_wlk, n_bin, bin_l, n_tot_wlk
   implicit none
   real(dp),       save, public :: da
   real(dp),       save, public :: eps
   real(dp),       save, public :: eps_force_reg
   real(dp),       save, public :: cg_acc
   integer(int32), save, public :: n_opt_step
   integer(int32), save, public :: n_crs_step
   integer(int32), save, public :: n_avg_step
   real(dp), save, public :: n_tot_sampling, n_red_sampling
   real(dp), save, public :: n_tot_bins, n_red_bins
   integer(int32), save, public :: n_sampling_task
   character(3),   save, public :: opt_mthd
   character(3),   save, public :: lin_solv_mthd
   logical,  save, public       :: stop_opt
   logical,  save, public       :: corr_samp
   logical,  save, public       :: eps_auto_tune
   real(dp), save, public, allocatable, dimension(:,:)   :: e_all
   real(dp), save, public, allocatable, dimension(:,:,:) :: O_all
   !real(dp), save, public, allocatable, dimension(:,:,:) :: Oe_tot_measures
   real(dp), save, public, allocatable, dimension(:,:,:) :: f_b
   real(dp), save, public, allocatable, dimension(:)     :: O
   real(dp), save, public, allocatable, dimension(:)     :: Oe
   real(dp), save, public, allocatable, dimension(:)     :: f
   real(dp), save, public, allocatable, dimension(:)     :: f_err
   real(dp), save, public, allocatable, dimension(:,:)   :: f_cov
   integer(int32),  save, public, allocatable, dimension(:) :: f_indx
   real(dp), save, public :: e_tot, e_tot_old
   real(dp), save, public :: e_err, e_var, e_err_old
   real(dp), save, public :: de_estimated
   real(dp), save, public :: force_cut, err_fct
   real(dp), save, public :: norm_da_vec, dev_max_f, covariance_prefactor
   real(dp), save, public :: norm_cut, norm_min, alpha_norm
   real(dp), save, public, allocatable, dimension(:) :: vec_wvfn_par_avg
   real(dp), save, public, allocatable, dimension(:) :: vec_wvfn_par_var
   real(dp), save, public, allocatable, dimension(:) :: vec_wvfn_par_var_old
   integer(int32), public, save   :: n_par_opt
   integer(int32), public, save   :: n_par_opt_s
   integer(int32), public, save   :: n_par_sgn
   integer(int32), public, save   :: i_opt_step, i_crs_step
   real(dp), allocatable, public, save, dimension(:,:) :: pip0
   real(dp), allocatable, public, save, dimension(:,:) :: pip1
   real(dp), allocatable, public, save, dimension(:)   :: pip2
   real(dp), allocatable, public, save, dimension(:,:) :: pip3
   public :: Init_WavefunctionOptimization_Variables
contains
   subroutine Init_WavefunctionOptimization_Variables()
      n_sampling_task= n_bin*bin_l
      allocate( pip0(n_par_opt_s,n_wlk), pip1(n_par_opt_s,n_wlk))
      allocate( pip2(n_wlk))
      allocate( pip3(n_bin,n_wlk))
      allocate( e_all(1:n_sampling_task,1:n_wlk) )                ; e_all = 0.0_dp
      allocate( O_all(1:n_par_opt_s,1:n_sampling_task,1:n_wlk) )  ; O_all = 0.0_dp
      allocate( f_b(1:n_par_opt_s,1:n_sampling_task,1:n_wlk) )  ; f_b = 0.0_dp
      allocate( O(1:n_par_opt_s) )  ; O = 0.0_dp
      allocate( Oe(1:n_par_opt_s) ) ; Oe = 0.0_dp
      allocate( f(1:n_par_opt_s) )  ; f = 0.0_dp
      allocate( f_err(1:n_par_opt_s) )   ; f_err = 0.0_dp
      allocate( f_indx(1:n_par_opt_s) )  ; f_indx = 0_int32
      if (opt_mthd.ne.'ams'.and.trim(lin_solv_mthd).ne.'cg') then
         allocate( f_cov(1:n_par_opt_s,1:n_par_opt_s) ) ; f_cov = 0.0_dp
      endif
      allocate( vec_wvfn_par_var(1:n_par_opt_s) ) ; vec_wvfn_par_var = 0.0_dp
      if (n_avg_step.ge.2_int32) then
         allocate( vec_wvfn_par_avg(1:n_par_opt_s) ) ; vec_wvfn_par_avg = 0.0_dp
      endif
      n_tot_sampling = dble(n_tot_wlk) * dble(n_bin) * dble(bin_l)
      n_red_sampling = n_tot_sampling - 1.0_dp
      n_tot_bins = dble(n_tot_wlk) * dble(n_bin)
      n_red_bins = n_tot_bins - 1.0_dp 
      i_opt_step = 0_int32
      i_crs_step = 0_int32
      de_estimated = 0.0_dp
   end subroutine Init_WavefunctionOptimization_Variables
end module wavefunction_optimization_v
