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
module wavefunction_optimization_m
   use fortran_kinds_v, only: int32, dp, stdout
   use openmp_mpi_m
   use write_lines_m
   use timings_m
   use molecular_system_v,         only: molsys_prs
   use qdo_system_v,               only: qdosys_prs
   use wavefunction_optimization_v
   use correlated_sampling_m,      only: Init_CorrelatedSampling, Fnlz_CorrelatedSampling, &
   & crrsmp, W_avg, Ovrlap_avg
   use io_chgoptpar_m,             only: read_chng_fle, write_chng_fle, delete_chng_fle
   use quantum_monte_carlo_v,      only: n_bin, n_wlk, bin_l, n_tot_wlk, n_bra, n_trm
   use fermionic_orbitals_m,       only: aoc_opt, aoe_opt, n_par_forbs
   use feporb_mod,                 only: poc_opt, poe_opt, n_par_eporbs
   use jastrow_orbitals_m,         only: jc_opt, je_opt, n_par_jorbs
   use drudonic_orbitals_m,        only: qc_opt, qe_opt, n_par_dorbs
   use jastrow_factors_v,                 only: jc1_opt, jc2_opt, jce_opt, jd1_opt, jd2_opt,&
   & n_par_jst, n_par_jstc, n_par_jstd
   use fermionic_wavefunction_v,   only: dl_opt, n_par_frm, n_par_det_e, n_par_det_p,&
   & n_par_det_ep, n_par_frm_wvfn_s, n_par_frm_wvfn
   use qdo_wavefunction_v,         only: ql_opt, n_par_drd, wf_type_d, n_par_drd_wvfn
   use jstqepar_var,               only: jqe_opt, n_par_jstqe
   use acceptance_rate_m,          only: upd_var_dt, cmp_acc_rates, write_acc_rates
   use optimization_averages_m,    only: Comp_Jackknife_Averages, Comp_Jackknife_Weighted_Averages, &
   & Comp_Jackknife_Covariances, Comp_MaxForce_Deviation,&
   & Simple_Sampling_Accumulation, Correlated_Sampling_Accumulation
   use fermionic_wavefunction_m,   only: save_fermionic_wavefunction_v, upd_fermionic_wavefunction_v
   use qdo_wavefunction_m,         only: save_qdowvf_par, upd_qdowvf_par
   implicit none
   public  :: Init_WavefunctionOptimization, Exec_WavefunctionOptimization, Fnlz_WavefunctionOptimization
   private :: comp_parameter_update, comp_CS_optimization_steps, cmp_norm_da_vec
contains
   subroutine Init_WavefunctionOptimization()
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Parameter step amplitude", da, var_name='da')
      call write_variable_line(stdout,0,mpi_rank,2,"Parameter for matrix regularization", eps, var_name='eps')
      call write_variable_line(stdout,0,mpi_rank,2,"Parameter for force regularization", eps_force_reg, var_name='eps_force_reg')
      call write_variable_line(stdout,0,mpi_rank,2,"Maximum force considered in optimization", force_cut, var_name='force_cut')
      call write_variable_line(stdout,0,mpi_rank,2,"Normalization limit", norm_cut, var_name='norm_cut')
      call write_variable_line(stdout,0,mpi_rank,2,"Minimum Normalization limit", norm_min, var_name='norm_min')
      call write_variable_line(stdout,0,mpi_rank,2,"Alpha normalization mix", alpha_norm, var_name='alpha_norm')
      if(molsys_prs) then
         call write_empty_line(stdout,0,mpi_rank)
         if(dl_opt)  call write_variable_line(stdout,0,mpi_rank,2,"Optimizing fermionic determinant", dl_opt, var_name='dl_opt')
         if(aoc_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing electronic orbital coefficients", aoc_opt, var_name='aoc_opt')
         if(aoe_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing electronic orbital exponents", aoe_opt, var_name='aoe_opt')
         if(poc_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing positronic orbital coefficients", poc_opt, var_name='poc_opt')
         if(poe_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing positronic orbital exponents", poe_opt, var_name='poe_opt')
         call write_empty_line(stdout,0,mpi_rank)
         if(jc1_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing One Body cusp", jc1_opt, var_name='jc1_opt')
         if(jc2_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing Two Body cusp", jc2_opt, var_name='jc2_opt')
         if(jce_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing cusps exponential parameters", jce_opt, var_name='jce_opt')
         if(jd1_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing dynamical one body Jastrow terms", jd1_opt, var_name='jd1_opt')
         if(jd2_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing dynamical two body Jastrow terms", jd2_opt, var_name='jd2_opt')
         if(jc_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing Jastrow orbitals coefficients", jc_opt, var_name='jc_opt')
         if(je_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing Jastrow orbitals exponentials", je_opt, var_name='je_opt')
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of elec. determinant/pfaffian parameters", n_par_det_e)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of posi. determinant/pfaffian parameters", n_par_det_p)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of elec.-posi. geminal parameters", n_par_det_ep)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of atomic orbitals parameters", n_par_forbs)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of positronic orbitals parameters", n_par_eporbs)
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of cusp parameters", n_par_jstc)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of dynamical Jastrow parameters", n_par_jstd)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of Jastrow orbitals parameters", n_par_jorbs)
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_line(stdout,0,mpi_rank,2,"Total number of ferimionic parameters", n_par_frm)
         call write_variable_line(stdout,0,mpi_rank,2,"Total number of Jastrow parameters", n_par_jst)
      endif
      if (qdosys_prs) then
         call write_empty_line(stdout,0,mpi_rank)
         select case ( abs(wf_type_d) )
          case(5)
            if(ql_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing drudonic product of orbitals", ql_opt, var_name='ql_opt')
            if(qc_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing drudonic orbital exponents", qc_opt, var_name='qc_opt')
            if(qe_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing drudonic orbital exponents", qe_opt, var_name='qe_opt')
          case(6)
            if(ql_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing drudonic dipole function", ql_opt, var_name='ql_opt')
          case default
         end select
         if(jqe_opt) call write_variable_line(stdout,0,mpi_rank,2,"Optimizing electron-drudons coupling", jqe_opt, var_name='jqe_opt')
         call write_empty_line(stdout,0,mpi_rank)
         select case ( abs(wf_type_d) )
          case(5)
            call write_variable_line(stdout,0,mpi_rank,2,"Number of linear drudonic parameters ", n_par_drd)
            call write_variable_line(stdout,0,mpi_rank,2,"Number of drudonic orbitals parameters", n_par_dorbs)
          case(6)
            call write_variable_line(stdout,0,mpi_rank,2,"Number of drudonic dipole wave function parameters", n_par_drd)
          case default
         end select
         call write_variable_line(stdout,0,mpi_rank,2,"Number of drudonic dipole wave function parameters", n_par_drd_wvfn)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of drudons-electrons coupling parameters", n_par_jstqe)
      endif
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of parameters to optimize", n_par_opt)
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of parameters after symmetrization", n_par_opt_s)
      if (n_bin*bin_l*n_wlk*n_mpi_tasks.lt.10*n_par_opt_s) then
         call write_variable_line(stdout,0,mpi_rank,2,"WARNING!!! n_bin < 10*n_par. Please increase n_bin to", 10*n_par_opt_s / (bin_l*n_wlk*n_mpi_tasks) + 1)
      endif
      call Init_WavefunctionOptimization_Variables()
      if ( corr_samp ) call Init_CorrelatedSampling( bin_l, n_bin, n_wlk)
      call write_chng_fle( )
   end subroutine Init_WavefunctionOptimization
   subroutine Exec_WavefunctionOptimization
      if (n_par_opt_s.eq.0_int32) then
         call write_empty_line(stdout,0,mpi_rank)
         call write_simple_line(stdout,0,mpi_rank,2,"l","WARNING!!! No parameters to optimize, stopping the run.")
         return
      endif
      do while( i_opt_step.lt.n_opt_step )
         call time_compute(t_opt_step)
         i_opt_step = i_opt_step + 1_int32
         call write_separator_line(stdout,0,mpi_rank,2,"_")
         call write_variable_line(stdout,0,mpi_rank,2,"Optimization step", i_opt_step)
         call Simple_Sampling_Accumulation( )
         call Comp_Jackknife_Averages( )
         call Comp_MaxForce_Deviation( )
         call Comp_Jackknife_Covariances(weighted_covariance=.false.)
         call comp_parameter_update( )
         call cmp_acc_rates()
         call write_acc_rates()
         call upd_var_dt(prnt=.true.)
         !call read_chng_fle( stop_opt )
         call time_compute_difference(t_opt_step,time_cpu_tot=t_opt_tot)
         call write_variable_line(stdout,0,mpi_rank,2,"Optimization step time",t_opt_step,units='sec.' )
         if( stop_opt ) exit
         if ( corr_samp ) call comp_CS_optimization_steps()
      enddo ! n_opt_steps (main optimization steps)
   end subroutine Exec_WavefunctionOptimization
   subroutine comp_CS_optimization_steps( )
      Ovrlap_avg = 1.0_dp
      i_crs_step = 0_int32
      do while ( (i_crs_step.lt.n_crs_step) .and. (i_opt_step.lt.n_opt_step) &
      &.and. (abs(Ovrlap_avg).ge.0.99) )
         call time_compute(t_crs_step)
         i_opt_step = i_opt_step + 1_int32
         i_crs_step = i_crs_step  + 1_int32
         call write_simple_line(stdout,0,mpi_rank,2,"r","___________________________________________________")
         call write_variable_line(stdout,0,mpi_rank,2,"Correlated sampling optimization step",i_opt_step)
         call Correlated_Sampling_Accumulation()
         call Comp_Jackknife_Weighted_Averages(  )
         call Comp_MaxForce_Deviation( )
         call Comp_Jackknife_Covariances( weighted_covariance=.true.)
         call comp_parameter_update( )
         call time_compute_difference(t_crs_step,time_cpu_tot=t_crs_tot)
         call write_variable_line(stdout,0,mpi_rank,2,"Optimization step time",t_crs_step,units='sec.' )
      enddo
      !call write_chng_fle()
   end subroutine comp_CS_optimization_steps
   subroutine Fnlz_WavefunctionOptimization( )
      call delete_chng_fle( )
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","OPTIMIZATION ENDED")
      call write_empty_line(stdout,0,mpi_rank)
      if (n_avg_step.ge.2_int32) then
         call write_variable_line(stdout,0,mpi_rank,2,"Averaging on the last parameters",n_avg_step,var_name='n_avg_step')
         if (molsys_prs) call upd_fermionic_wavefunction_v( vec_wvfn_par_avg(1:n_par_frm_wvfn_s) )
         if (qdosys_prs) call upd_qdowvf_par( vec_wvfn_par_avg(n_par_frm_wvfn_s+1:n_par_opt_s) )
      endif
      call write_simple_line(stdout,0,mpi_rank,2,"l","Saving wavefunction parameters...")
      if (mpi_rank.eq.0_int32) then
         if (molsys_prs) call save_fermionic_wavefunction_v()
         if (qdosys_prs) call save_qdowvf_par()
      endif
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"c","All well, farewell! See you soon!!")
      if ( corr_samp ) call Fnlz_CorrelatedSampling( n_wlk)
   end subroutine Fnlz_WavefunctionOptimization
   subroutine comp_parameter_update
      integer(int32) :: i1
      real(dp) :: r, da_tmp, f_da, de, de_err
      if (i_opt_step.gt.1_int32) then
         de = e_tot - e_tot_old
         de_err = sqrt(e_err**2+e_err_old**2)
      else
         de = e_tot
         de_err = e_err
      endif
      f_da = de / (err_fct*de_err)
      da_tmp = da
      if ( f_da.gt.1.0_dp .and.i_opt_step.gt.1_int32) then
         call write_simple_line(stdout,0,mpi_rank,2,"l","WARNING!!! Energy too high reverting optimization step")
         vec_wvfn_par_var = -0.99_dp * vec_wvfn_par_var_old
      else
         if (trim(lin_solv_mthd).eq.'cg') then
            call cg_lin_solv_mpi( cg_acc, n_par_opt_s, n_par_sgn, n_sampling_task, n_wlk, O_all(:,:,:), eps, f_err, f, vec_wvfn_par_var )
         else
            if( mpi_rank.eq.0_int32 ) then
               call lu_lin_solv( n_par_opt_s, n_par_sgn, f_cov(:,:), eps, f_err, f, vec_wvfn_par_var )
            endif
#if defined _MPI || defined _MPIh || defined _MPI08
            call mpi_barrier(MPI_COMM_WORLD, mpierr)
            call mpi_bcast(vec_wvfn_par_var(1:n_par_opt_s), n_par_opt_s, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
         endif
         call cmp_norm_da_vec( O_all(:,:,:) )
         if (norm_cut.gt.0.0_dp) then
            if ( norm_da_vec*da_tmp.ge.norm_cut) then
               call write_variable_line(stdout,0,mpi_rank,2,"WARNING!!! Slowing down to preserve norm. change", norm_cut / norm_da_vec)
               da_tmp = norm_cut / norm_da_vec
               if (da_tmp.lt.alpha_norm*da) then
                  da = (da + da_tmp) * 0.5_dp
                  if (da.lt.0.001_dp) da = 0.001_dp
                  call write_variable_line(stdout,0,mpi_rank,2,"Changing parameter step used", da, var_name='da' )
               endif
            else if (norm_da_vec*da_tmp.lt.alpha_norm*norm_cut.and. norm_cut.gt.norm_min) then
               norm_cut = (norm_cut+norm_da_vec*da_tmp) * 0.5_dp
               if (norm_cut.lt.norm_min) norm_cut = norm_min
               call write_variable_line(stdout,0,mpi_rank,2,"Changing norm. change cut-off to", norm_cut )
            endif
         endif
         de_estimated = -dot_product( f(1:n_par_sgn), vec_wvfn_par_var(1:n_par_sgn) )
         if (de_estimated.gt.0_int32) then
            call write_simple_line(stdout,0,mpi_rank,2,"l","WARNING!!! Rejecting move (DE > 0)")
            da_tmp = 0.0_dp
         endif
         call write_variable_line(stdout,0,mpi_rank,2,"Wavefunction norm change", norm_da_vec*da_tmp)
         call write_variable_line(stdout,0,mpi_rank,2,"Projected energy change", de_estimated*da_tmp  )
         call cmp_norm_da_vec( f_b(:,:,:) )
         call write_variable_line(stdout,0,mpi_rank,2,"Signal-to-noise ratio", abs(de_estimated) / norm_da_vec )
         call write_variable_line(stdout,0,mpi_rank,2,"Final parameter step used", da_tmp, var_name='da' )
         vec_wvfn_par_var(1:n_par_sgn) = vec_wvfn_par_var(1:n_par_sgn) * da_tmp
         de_estimated = de_estimated * da_tmp
         do i1 = 1_int32, n_par_sgn
            if (f_indx(i1).ne.i1) then
               vec_wvfn_par_var(f_indx(i1)) = vec_wvfn_par_var(i1)
               vec_wvfn_par_var(i1) = 0.0_dp
            endif
         enddo
      endif
      if (molsys_prs) call upd_fermionic_wavefunction_v( vec_wvfn_par_var(1:n_par_frm_wvfn_s) )
      if (qdosys_prs) call upd_qdowvf_par( vec_wvfn_par_var(n_par_frm_wvfn_s+1:n_par_opt_s) )
      if (n_avg_step.ge.2_int32) then
         if (i_opt_step.ge.n_opt_step-(n_avg_step-2)) then
            vec_wvfn_par_avg(1:n_par_opt_s) = vec_wvfn_par_avg(1:n_par_opt_s) &
            &- dble((n_avg_step-1)-n_opt_step+i_opt_step)/ dble(n_avg_step) * vec_wvfn_par_var(1:n_par_opt_s)
         endif
      endif
      e_tot_old = e_tot ; e_err_old = e_err
      vec_wvfn_par_var_old = vec_wvfn_par_var
   end subroutine comp_parameter_update
   subroutine cmp_norm_da_vec( this_O )
      real(dp), dimension(n_par_opt_s,n_sampling_task,n_wlk), intent(inout) :: this_O
      integer(int32) :: iw
      real(dp), dimension(n_sampling_task) :: par_tmp
      norm_da_vec = 0.0_dp
!$omp parallel default(shared) private(iw,par_tmp) reduction(+:norm_da_vec)
!$omp do
      do iw = 1_int32, n_wlk
         call dgemv( 'T', n_par_opt_s, n_sampling_task, 1.0_dp, this_O(1:n_par_opt_s,1:n_sampling_task,iw), n_par_opt_s, vec_wvfn_par_var(1:n_par_opt_s), 1_int32, 0.0_dp, par_tmp(:), 1_int32)
         norm_da_vec = norm_da_vec + dot_product ( par_tmp(:), par_tmp(:))
      enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE, norm_da_vec ,     1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      norm_da_vec = dsqrt(norm_da_vec)
   end subroutine cmp_norm_da_vec
end module wavefunction_optimization_m
