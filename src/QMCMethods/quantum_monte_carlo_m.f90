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
module quantum_monte_carlo_m
   use fortran_kinds_v,               only: dp, int32, stdout
   use openmp_mpi_m,                  only: mpi_rank, n_mpi_tasks, fnlz_ompmpi_env
   use write_lines_m
   use wavefunction_optimization_v,   only: opt_mthd
   use quantum_monte_carlo_v
   use monte_carlo_averages_m,        only: init_qmcsmp
   use molecular_system_v,            only: n_fe, n_at, molsys_prs, n_el, n_po
   use qdo_system_v,                  only: n_qdo, qdosys_prs
   use molecular_system_m,            only: frm_cnf
   use fermionic_positions_m,         only: save_frmdst
   use acceptance_rate_m,             only: ini_acc_rates, rst_acc_rates
   use variational_monte_carlo_m,     only: Init_VariationalMonteCarlo, &
   & Exec_VariationalMonteCarlo, Fnlz_VariationalMonteCarlo, Thrm_VariationalMonteCarlo
   use diffusion_monte_carlo_m,       only: Init_DiffusionMonteCarlo, &
   & Exec_DiffusionMonteCarlo, Fnlz_DiffusionMonteCarlo, Thrm_DiffusionMonteCarlo
   use wavefunction_optimization_m,   only: Init_WavefunctionOptimization, &
   & Exec_WavefunctionOptimization, Fnlz_WavefunctionOptimization
   use monte_carlo_checkpointfiles_m, only: load_checkpointfiles
   implicit none
   public :: init_qmc_run, exec_qmc_run
contains
   subroutine init_qmc_run()
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      select case ( opt_mthd )
       case('non')
         select case ( qmc_mthd )
          case('dmc')
            call write_simple_line(stdout,0,mpi_rank,2,"c","DIFFUSION MONTE CARLO")
          case('vmc')
            call write_simple_line(stdout,0,mpi_rank,2,"c","VARIATIONAL MONTE CARLO")
          case default
            call write_simple_line(stdout,0,mpi_rank,2,"c","WARNING!!! NOTHING TO BE RAN")
         end select
       case('src')
         call write_simple_line(stdout,0,mpi_rank,2,"c","STOCHASTIC RECONFIGURATION")
       case('snr')
         call write_simple_line(stdout,0,mpi_rank,2,"c","SIGNAL-TO-NOISE RATIO")
       case('ams')
         call write_simple_line(stdout,0,mpi_rank,2,"c","AMSGRAD")
       case('lma')
         call write_simple_line(stdout,0,mpi_rank,2,"c","LEVENBERG-MARQUARDT")
       case('sgd')
         call write_simple_line(stdout,0,mpi_rank,2,"c","STOCHASTIC GRADIENT DESCENT")
      end select
      if ( n_bra .eq. 0_int32 ) then
         n_bra_vmc = 2_int32
      else
         n_bra_vmc = n_bra
      endif
      if ( bin_l.eq.0_int32 ) bin_l = int((n_fe+n_qdo)/2)+1
      if ( n_scr.eq.0_int32 ) n_scr = n_fe + n_qdo
      n_tot_smpl =  bin_l * n_bin
      bin_l_old = 0_int32
      n_bin_old = 0_int32
      n_tot_smpl_old = 0_int32
      n_tot_wlk = n_wlk * n_mpi_tasks
      G_w_trgt  = dble( n_tot_wlk )
      call ini_acc_rates()
      call write_simple_line(stdout,0,mpi_rank,2,"c","________________________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","Variational Monte Carlo sampling")
      call write_empty_line(stdout,0,mpi_rank)
      select case( smp_type )
       case(0) ! Simple sampling 1d
         call write_simple_line(stdout,0,mpi_rank,2,"l","1D Uniform Transition prob.[-dt:dt] (random direction).")
       case(1) ! Gaussian move 1d
         call write_simple_line(stdout,0,mpi_rank,2,"l","1D Gaussian Transition prob. (dt variance) (random direction).")
       case(2) ! Simple sampling 3d
         call write_simple_line(stdout,0,mpi_rank,2,"l","3D Uniform Transition prob. in volume (2dt)^3.")
       case(3) ! Gaussian move 3d
         call write_simple_line(stdout,0,mpi_rank,2,"l","3D Gaussian Transition prob. (dt*I covariance).")
       case(10) ! Langevin dynamics.
         call write_simple_line(stdout,0,mpi_rank,2,"l","3D Transtion probability for Metropolis-adjusted Langevin.")
         if (nuc_corr) then
            call write_simple_line(stdout,0,mpi_rank,2,"l","Adding Umrigar Nightingal and Runge corrections to nuclei.")
         endif
      end select
      if (two_step_sampling) &
      & call write_simple_line(stdout,0,mpi_rank,2,"l","Two step sampling activated")
      if (spc_map) then
         call write_simple_line(stdout,0,mpi_rank,2,"l","Adding space mapping to dynamically adjust dt step.")
      endif
      call write_empty_line(stdout,0,mpi_rank)
      if ( n_el.gt.0_int32 ) &
      & call write_variable_line(stdout,0,mpi_rank,2,"Initial Metropolis step for electrons", dt_e,var_name="dt_e")
      if ( n_po.gt.0_int32 ) &
      & call write_variable_line(stdout,0,mpi_rank,2,"Initial Metropolis step for positrons", dt_e,var_name="dt_p")
      if ( n_qdo.gt.0_int32 ) &
      & call write_variable_line(stdout,0,mpi_rank,2,"Initial Metropolis step for drudons", dt_e,var_name="dt_q")
      call write_variable_line(stdout,0,mpi_rank,2,"Variable metropolis step", var_dt,var_name="var_dt")
      if ( opt_mthd.eq.'non' ) call init_qmcsmp( )
      select case ( opt_mthd )
       case('non')
         select case ( qmc_mthd )
          case('dmc')
            call Init_DiffusionMonteCarlo()
          case('vmc')
            call Init_VariationalMonteCarlo()
          case default
         end select
       case('src','snr','ams','lma','sgd','res')
         call Init_WavefunctionOptimization()
       case default
      end select
      if (opt_mthd.eq.'non') call load_checkpointfiles() 
      select case ( opt_mthd )
       case('non')
         if (save_status.eq.0_int32) then
            call Thrm_VariationalMonteCarlo()
         endif
         select case ( qmc_mthd )
          case('dmc')
            if (save_status.le.3_int32) then
               call Thrm_DiffusionMonteCarlo()
            endif
          case default
         end select
       case('src','snr','ams','lma','sgd')
         if (save_status.eq.0_int32) then
            call Thrm_VariationalMonteCarlo()
         endif
       case default
      end select
   end subroutine init_qmc_run
   subroutine exec_qmc_run()
      select case ( opt_mthd )
       case('non')
         select case ( qmc_mthd )
          case('vmc')
            if (save_status.eq.1_int32.or.save_status.eq.2_int32) then
               call Exec_VariationalMonteCarlo()
            endif
          case('dmc')
            if (save_status.eq.3_int32.or.save_status.eq.4_int32) then
               call Exec_DiffusionMonteCarlo()
            endif
          case default
         end select
       case('src','snr','ams','lma','sgd')
         call Exec_WavefunctionOptimization()
       case default
      end select
      select case ( opt_mthd )
       case('non')
         select case ( qmc_mthd )
          case('dmc')
            call Fnlz_DiffusionMonteCarlo()
          case('vmc')
            call Fnlz_VariationalMonteCarlo()
          case default
         end select
       case('src','snr','ams','lma','sgd','res')
         call Fnlz_WavefunctionOptimization()
       case default
      end select
      if (mpi_rank.eq.0_int32) then
        if( molsys_prs .and. n_at.gt.0_int32)  call save_frmdst(frm_cnf(1)%d_fn )
      endif
   end subroutine exec_qmc_run
end module quantum_monte_carlo_m
