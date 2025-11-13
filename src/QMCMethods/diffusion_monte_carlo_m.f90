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
module diffusion_monte_carlo_m
   use fortran_kinds_v,        only: int32, dp, stdout
   use openmp_mpi_m
   use write_lines_m
   use timings_m
   use molecular_system_v,     only: n_fe, molsys_prs, n_po, n_el
   use qdo_system_v,           only: qdosys_prs, n_qdo
   use molecular_system_m,     only: save_fermionic_confs
   use qdo_system_m,           only: save_drudonic_confs
   use quantum_monte_carlo_v,  only: n_bra, bin_l, n_bin, n_wlk, n_tot_wlk, dt, &
   & qmc_mthd, n_wlk_max, dt_e, dt_p, dt_d, var_dt, save_status, save_step, n_scr
   use branching_m,            only: Init_Branching, Fnlz_Branching, Exec_Branching, chk_brnchn
   use local_energy_m,         only: v_ext
   use monte_carlo_averages_m
   use diffusion_monte_carlo_v
   use dmc_samplings_m
   use acceptance_rate_m,      only: rst_acc_rates, upd_var_dt, cmp_acc_rates, write_acc_rates
   use vmc_sampling_m
   use mersenne_twister19937_m,       only: save_mt19937
   use monte_carlo_checkpointfiles_m, only: save_checkpointfiles, load_checkpointfiles, &
   & clean_checkpointfiles
   implicit none
   public  :: Init_DiffusionMonteCarlo, Thrm_DiffusionMonteCarlo, &
   & Exec_DiffusionMonteCarlo, Fnlz_DiffusionMonteCarlo
contains
   subroutine Init_DiffusionMonteCarlo()
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of walkers", n_tot_wlk,var_name="n_tot_wlk")
      call write_variable_line(stdout,0,mpi_rank,2,"Number of independent walkers per MPI task", n_wlk,var_name="n_wlk")
      call write_variable_line(stdout,0,mpi_rank,2,"Monte Carlo sampling per core", n_tot_smpl,var_name="n_tot_smpl")
      call write_variable_line(stdout,0,mpi_rank,2,"Particle sweeps before branching", n_bra,var_name="n_bra")
      call write_variable_line(stdout,0,mpi_rank,2,"Number of moves before matrix recomputation", n_scr,var_name="n_scr")
      dt_eff_avg = dt
      call write_simple_line(stdout,0,mpi_rank,2,"c","________________________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","Diffusion Monte Carlo parameters")
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of thermalization steps", n_trm_dmc*bin_l, var_name='n_trm_dmc*bin_l')
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of sampling steps", n_tot_smpl, var_name='bin_l*n_blocks')
      call write_variable_line(stdout,0,mpi_rank,2,"Time discretization step", dt, var_name='dt', units='au')
      call write_variable_line(stdout,0,mpi_rank,2,"Thermalization time", dt*dble(n_trm_dmc*bin_l), units='au')
      call write_variable_line(stdout,0,mpi_rank,2,"Sampling time", dt*dble(n_tot_smpl), units='au')
      call write_variable_line(stdout,0,mpi_rank,2,"Total target weight of the configurations", G_w_trgt)
      call write_variable_line(stdout,0,mpi_rank,2,"Reference energy rescaling factor", g_r,units="Eh")
      select case(brnch_type)
       case(1)
         call write_variable_line(stdout,0,mpi_rank,2,"Branching with variable population", brnch_type,var_name='brnch_type')
       case(2)
         call write_variable_line(stdout,0,mpi_rank,2,"Stochastic Reconfiguration with fixed population", brnch_type,var_name='brnch_type')
       case default
         call write_variable_line(stdout,0,mpi_rank,2,"No branching ", brnch_type,var_name='brnch_type')
      end select
      if (brnch_type.gt.0_int32) then
         call write_variable_line(stdout,0,mpi_rank,2,"Number of configurational updates between branching", n_bra,var_name='n_bra')
      else
         call write_simple_line(stdout,0,mpi_rank,2,"l","Automatic branching between configurational updates")
      endif
      select case(encrr_type)
       case(1)
         call write_variable_line(stdout,0,mpi_rank,2,"ZSGMA 2016 energy cut-off",encrr_type,var_name='encrr_type')
         call write_variable_line(stdout,0,mpi_rank,2,"Cut-off rescaling factor",g_e,var_name='g_e')
       case(2)
         call write_variable_line(stdout,0,mpi_rank,2,"ZSGMA 2016 energy cut-off, based on variance",encrr_type,var_name='encrr_type')
       case(3)
         call write_variable_line(stdout,0,mpi_rank,2,"UNR 1993 energy cut-off",encrr_type,var_name='encrr_type')
       case(4)
         call write_variable_line(stdout,0,mpi_rank,2,"AU 2021 energy cut-off",encrr_type,var_name='encrr_type')
       case(11)
         call write_variable_line(stdout,0,mpi_rank,2,"ZSGMA 2016 energy cut-off (p.p.)",encrr_type,var_name='encrr_type')
         call write_variable_line(stdout,0,mpi_rank,2,"Cut-off rescaling factor",g_e,var_name='g_e')
       case(12)
         call write_variable_line(stdout,0,mpi_rank,2,"ZSGMA 2016 energy cut-off, based on variance (p.p.)",encrr_type,var_name='encrr_type')
       case(13)
         call write_variable_line(stdout,0,mpi_rank,2,"UNR 1993 energy cut-off (p.p.)",encrr_type,var_name='encrr_type')
       case(14)
         call write_variable_line(stdout,0,mpi_rank,2,"AU 2021 energy cut-off (p.p.)",encrr_type,var_name='encrr_type')
       case(0)
         call write_variable_line(stdout,0,mpi_rank,2,"No energy cut-off in weight",encrr_type,var_name='encrr_type')
      end select
      call allocate_dmcmthd_var( n_fe+n_qdo, n_wlk, n_wlk_max, n_pp_obs, bin_l )
      call Init_Branching()
   end subroutine Init_DiffusionMonteCarlo
   subroutine Thrm_DiffusionMonteCarlo()
      integer(int32) :: i1, i2, i3
      logical  :: do_branching
      real(dp) :: G_w_tot
      obs_b(:,:) = 0.0_dp
      if (allocated( obs_pp_b ) ) obs_pp_b(:,:,:) = 0.0_dp
      do i1 = 1_int32, bin_l
         call vmc_move_all_particles_and_walkers()
         if (encrr_type.ge.10_int32) then
            call bld_obsvec( obs_b(:,i1), obs_pp_b(:,:,i1))
            call accu_qmcavg(obs_b(:,i1), obs_pp_b(:,:,i1))
         else
            call bld_obsvec( obs_b(:,i1))
            call accu_qmcavg( obs_b(:,i1) )
         endif
      enddo
      var_dt = .false.
      dt_e = dt; dt_p = dt; dt_d = dt
      e_best       => obs_avg(2)
      e_best_var   => obs_var(2)
      if (encrr_type.ge.20_int32) then
         spe_best     => obs_pp_avg(:,2)
         spe_best_var => obs_pp_var(:)
      else if (encrr_type.ge.10_int32) then
         spe_best     => obs_pp_avg(:,1)
         spe_best_var => obs_pp_var(:)
      endif
      call rst_acc_rates()
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_error_dble_line(stdout,0,mpi_rank,2,"E_VMC",e_avg,sqrt(e_var / dble(n_tot_wlk)/ dble(bin_l)),units="Eh")
      call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Var[E_VMC]",e_var,e_var*sqrt(2.0_dp /(dble(n_tot_wlk)*dble(bin_l)-1)),units="Eh^2")
      if (encrr_type.ge.20_int32) then
         call write_variable_error_dble_line(stdout,0,mpi_rank,2,"PP E_VMC",sum(spe_best)/dble(n_pp_obs) , sqrt(sum(spe_best_var * obs_pp_avg(:,3) / obs_pp_avg(:,1)**2)/dble(n_pp_obs) ),units="Eh")
         call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Var[PP E_VMC]",sum(spe_best_var)/dble(n_pp_obs) , sum(spe_best_var*sqrt( 2.0_dp / ( obs_pp_avg(:,1)**2 / obs_pp_avg(:,3) - 1.0_dp) ))/dble(n_pp_obs),units="Eh^2")
      else if (encrr_type.ge.10_int32) then
         call write_variable_error_dble_line(stdout,0,mpi_rank,2,"PP E_VMC",sum(spe_best)/dble(n_pp_obs) , sqrt(sum(spe_best_var)/dble(n_pp_obs) * w2_avg / w_avg**2  ),units="Eh")
         call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Var[PP E_VMC]",sum(spe_best_var)/dble(n_pp_obs) , sum(spe_best_var)/dble(n_pp_obs)*sqrt( 2.0_dp / ( w_avg**2 / w2_avg - 1.0_dp) ),units="Eh^2")
      endif
      e_r        = e_best
      e_r_avg    = e_r
      call dmc_initialize_all_particles_and_walkers()
      w(1:n_wlk) = 1.0_dp
      if(n_wlk_max.gt.n_wlk) w(n_wlk+1:n_wlk_max) = 0.0_dp
      obs_b(:,:) = 0.0_dp
      if (allocated( obs_pp_b) ) obs_pp_b(:,:,:) = 0.0_dp
      call init_dmc_rene()
      do i1 = 1_int32, n_trm_dmc
         call time_compute(t_dmc_step)
         call write_separator_line(stdout,0,mpi_rank,2,"_")
         call write_variable_line(stdout,0,mpi_rank,2,"Thermalization block", i1)
         t_bra = 0_int32
         do i2 = 1_int32, bin_l
            call dmc_move_all_particles_and_walkers()
            if (encrr_type.ge.10_int32) then
               call bld_obsvec( obs_b(:,i2), obs_pp_b(:,:,i2) )
            else
               call bld_obsvec( obs_b(:,i2) )
            endif
            if (n_bra.eq.0_int32) then
               call chk_brnchn( do_branching )
            else
               do_branching=.false.
               if ( mod(i2,n_bra).eq.0_int32) do_branching=.true.
            endif
            call cmp_acc_rates()
            if ( do_branching ) then
               if (brnch_type.eq.1_int32) then
                  call Exec_Branching( 2_int32 )
               else
                  call Exec_Branching( brnch_type )
               endif
               call dmc_initialize_all_particles_and_walkers()
               t_bra = t_bra + 1_int32
            endif
            w_blck     = sum(obs_b(1,1:bin_l))
            w2_blck    = sum(obs_b(n_obs-1,1:bin_l))
            e_blck     = sum(obs_b(2,1:bin_l)) / w_blck
            e_var_blck = sum(obs_b(n_obs,1:bin_l)) / w_blck
            e_var_blck = (e_var_blck - e_blck**2) * w_blck**2 / (w_blck**2 - w2_blck )
            if ( i1.gt.1_int32 ) dt_eff_avg = sum(obs_b(n_obs-2,1:bin_l)) / w_blck
            if (encrr_type.ge.10_int32) then
               spe_blck(:,:) = 0.0_dp
               spe_var_blck(:) = 0.0_dp
               do i3 = 1_int32, bin_l
                  spe_blck(:,:) = spe_blck(:,:) + obs_pp_b(:,:,i3)
               enddo
               if (encrr_type.ge.20_int32) then
                  spe_blck(:,2) = spe_blck(:,2) / spe_blck(:,1)
                  spe_blck(:,4) = spe_blck(:,4) / spe_blck(:,1)
                  spe_var_blck(:) = (spe_blck(:,4)  - spe_blck(:,2)**2) * spe_blck(:,1)**2 / (spe_blck(:,1)**2 - spe_blck(:,3) )
               else
                  spe_blck(:,:) = spe_blck(:,:) / w_blck
                  spe_var_blck(:) = (spe_blck(:,2)  - spe_blck(:,1)**2) * w_blck**2/ (w_blck**2 - w2_blck )
               endif
            endif
            call cmp_dmc_rene( e_r_avg, obs_b(1,i2), 1.0_dp )
            e_r_vec(i2) = e_r * obs_b(1,i2)
            e_r_avg = sum( e_r_vec(1:bin_l) ) / w_blck
         enddo ! i2 (bin_l)
         if ( i1.eq.1_int32 ) then
            e_best       => e_blck
            e_best_var   => e_var_blck
            if (encrr_type.ge.20_int32) then
               spe_best     => spe_blck(:,2)
               spe_best_var => spe_var_blck(:)
            else if (encrr_type.ge.10_int32) then
               spe_best     => spe_blck(:,1)
               spe_best_var => spe_var_blck(:)
            endif
         endif
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_error_dble_line(stdout,0,mpi_rank,2,"E_b",e_blck,sqrt(e_var_blck*w2_blck/w_blck**2),units="Eh")
         call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Var[E_b]",e_var_blck ,e_var_blck*sqrt(2.0_dp/(w_blck**2/w2_blck-1.0_dp)),units="Eh^2")
         if (encrr_type.ge.20_int32) then
            call write_empty_line(stdout,0,mpi_rank)
            call write_variable_error_dble_line(stdout,0,mpi_rank,2,"PP E_b",sum(spe_best)/dble(n_pp_obs) , sqrt(sum(spe_best_var * spe_blck(:,3) / spe_blck(:,1)**2)/dble(n_pp_obs) ),units="Eh")
            call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Var[PP E_b]",sum(spe_best_var)/dble(n_pp_obs) , sum(spe_best_var*sqrt( 2.0_dp / ( spe_blck(:,1)**2 / spe_blck(:,3) - 1.0_dp) ))/dble(n_pp_obs),units="Eh^2")
         else if (encrr_type.ge.10_int32) then
            call write_empty_line(stdout,0,mpi_rank)
            call write_variable_error_dble_line(stdout,0,mpi_rank,2,"PP E_b",sum(spe_best)/dble(n_pp_obs) , sqrt(sum(spe_best_var)/dble(n_pp_obs) * w2_blck / w_blck**2  ),units="Eh")
            call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Var[PP E_b]",sum(spe_best_var)/dble(n_pp_obs) , sum(spe_best_var)/dble(n_pp_obs)*sqrt( 2.0_dp / ( w_blck**2 / w2_blck - 1.0_dp) ),units="Eh^2")
         endif
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of branchings",t_bra)
         call write_variable_line(stdout,0,mpi_rank,2,"Reference Energy",e_r_avg,units='Eh')
         call write_variable_line(stdout,0,mpi_rank,2,"Average effective time step",dt_eff_avg,units='au')
         call write_variable_line(stdout,0,mpi_rank,2,"Total number of configurations",n_tot_wlk)
         call write_variable_line(stdout,0,mpi_rank,2,"Rescaled total weight",obs_b(1,bin_l))
         call write_acc_rates( )
         call rst_acc_rates()
         call time_compute_difference(t_dmc_step,time_cpu_tot=t_therm_dmc)
         call write_variable_line(stdout,0,mpi_rank,2,"Block time",t_dmc_step,units='sec.')
      enddo
      call save_checkpointfiles( save_status_in = 3_int32, step_in=0_int32, create_backup=.true. )
      call empty_qmcsmp()
   end subroutine Thrm_DiffusionMonteCarlo
   subroutine Exec_DiffusionMonteCarlo()
      integer(int32) :: i1, i2, i3
      logical :: do_branching
      restart_walker(:) = .true.
      call dmc_initialize_all_particles_and_walkers()
      e_best       => obs_avg(2)
      e_best_var   => obs_var(2)
      if (encrr_type.ge.20_int32) then
         spe_best     => obs_pp_avg(1:n_pp_obs,2)
         spe_best_var => obs_pp_var(1:n_pp_obs)
      elseif (encrr_type.ge.10_int32) then
         spe_best     => obs_pp_avg(1:n_pp_obs,1)
         spe_best_var => obs_pp_var(1:n_pp_obs)
      endif
      do i1 = 1_int32, n_bin
         call time_compute(t_dmc_step)
         call write_separator_line(stdout,0,mpi_rank,2,"_")
         call write_variable_line(stdout,0,mpi_rank,2,"Accumulation block", i1)
         t_bra = 0_int32
         do i2 = 1_int32, bin_l
            call dmc_move_all_particles_and_walkers()
            if (encrr_type.gt.10_int32) then
               call bld_obsvec( obs_b(:,i2), obs_pp_b(:,:,i2) )
            else
               call bld_obsvec( obs_b(:,i2) )
            endif
            if (n_bra.eq.0_int32) then
               call chk_brnchn( do_branching )
            else
               do_branching=.false.
               if ( mod(i2,n_bra).eq.0_int32) do_branching=.true.
            endif
            call cmp_acc_rates()
            if ( do_branching ) then
               call Exec_Branching( brnch_type )
               call dmc_initialize_all_particles_and_walkers()
               t_bra = t_bra + 1_int32
            endif
            call cmp_dmc_rene( e_r_avg, obs_b(1,i2), 1.0_dp )
            if (encrr_type.gt.10_int32) then
               call accu_qmcavg(obs_b(:,i2), obs_pp_b(:,:,i2))
            else
               call accu_qmcavg(obs_b(:,i2) )
            endif
            w_blck     = sum(obs_b(1,1:bin_l))
            e_r_vec(i2) = e_r * obs_b(1,i2)
            e_r_avg = sum( e_r_vec(1:bin_l) ) / w_blck
            dt_eff_avg = sum(obs_b(n_obs-2,1:bin_l)) / w_blck
         enddo ! i2 (bin_l)
         w_blck     = sum(obs_b(1,1:bin_l))
         w2_blck    = sum(obs_b(n_obs-1,1:bin_l))
         e_blck     = sum(obs_b(2,1:bin_l)) / w_blck
         e_var_blck = sum(obs_b(n_obs,1:bin_l)) / w_blck
         e_var_blck = (e_var_blck - e_blck**2) * w_blck**2 / (w_blck**2 - w2_blck )
         if (encrr_type.ge.10_int32) then
            spe_blck(:,:) = 0.0_dp
            spe_var_blck(:) = 0.0_dp
            do i3 = 1_int32, bin_l
               spe_blck(:,:) = spe_blck(:,:) + obs_pp_b(:,:,i3)
            enddo
            if (encrr_type.ge.20_int32) then
               spe_var_blck(:) = (spe_blck(:,3) * spe_blck(:,1) - spe_blck(:,1)**2) / (spe_blck(:,1)**2 - spe_blck(:,4) )
            else
               spe_var_blck(:) = (spe_blck(:,2) * w_blck - spe_blck(:,1)**2) / (w_blck**2 - w2_blck )
            endif
         endif
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_error_dble_line(stdout,0,mpi_rank,2,"E_b",e_blck,sqrt(e_var_blck*w2_blck/w_blck**2),units="Eh")
         call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Var[E_b]",e_var_blck ,e_var_blck*sqrt(2.0_dp/(w_blck**2/w2_blck-1.0_dp)),units="Eh^2")
         if (encrr_type.ge.20_int32) then
            call write_empty_line(stdout,0,mpi_rank)
            call write_variable_error_dble_line(stdout,0,mpi_rank,2,"PP E_MA",sum(spe_best)/dble(n_pp_obs) , sqrt(sum(spe_best_var(:) * obs_pp_avg(:,3) / obs_pp_avg(:,1)**2)/dble(n_pp_obs) ),units="Eh")
            call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Var[PP E_MA]",sum(spe_best_var)/dble(n_pp_obs) , sum(spe_best_var(:)*sqrt( 2.0_dp / ( obs_pp_avg(:,1)**2 / obs_pp_avg(:,3) - 1.0_dp) ))/dble(n_pp_obs),units="Eh^2")
         else if (encrr_type.ge.10_int32) then
            call write_empty_line(stdout,0,mpi_rank)
            call write_variable_error_dble_line(stdout,0,mpi_rank,2,"E[PP E_MA]",sum(spe_best)/dble(n_pp_obs) , sqrt(sum(spe_best_var)/dble(n_pp_obs) * w2_avg / w_avg**2  ),units="Eh")
            call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Var[PP E_MA]",sum(spe_best_var)/dble(n_pp_obs), sum(spe_best_var)/dble(n_pp_obs)*sqrt( 2.0_dp / ( w_avg**2 / w2_avg - 1.0_dp) ),units="Eh^2")
         endif
         call write_variable_error_dble_line(stdout,0,mpi_rank,2,"E_MA",e_avg,sqrt(e_var*w2_avg/w_avg**2),units="Eh")
         call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Var[E_MA]",e_var,e_var*sqrt(2.0_dp/(w_avg**2/w2_avg-1.0_dp)),units="Eh^2")
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_line(stdout,0,mpi_rank,2,"Number of branchings",t_bra)
         call write_variable_line(stdout,0,mpi_rank,2,"Reference Energy",e_r_avg,units='Eh')
         call write_variable_line(stdout,0,mpi_rank,2,"Average effective time step",dt_eff_avg,units='au')
         call write_variable_line(stdout,0,mpi_rank,2,"Total number of configurations",n_tot_wlk)
         call write_variable_line(stdout,0,mpi_rank,2,"Rescaled total weight",obs_b(1,bin_l))
         call write_acc_rates( )
         call time_compute_difference(t_dmc_step,time_cpu_tot=t_dmc_tot)
         call write_variable_line(stdout,0,mpi_rank,2,"Block time",t_dmc_step,units='sec.')
         call save_checkpointfiles( save_status_in = 4_int32, step_in=i1, create_backup=.false.  )
      enddo
   end subroutine Exec_DiffusionMonteCarlo
   subroutine Fnlz_DiffusionMonteCarlo()
      if (mpi_rank.eq.0_int32) call write_qmcavg( )
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","DIFFUSION MONTE CARLO CALCULATION ENDED")
      call deallocate_dmcmthd_var()
      call Fnlz_Branching()
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"c","All well, farewell! See you soon!!")
      call clean_checkpointfiles()
   end subroutine Fnlz_DiffusionMonteCarlo
end module diffusion_monte_carlo_m
