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
module variational_monte_carlo_m
   use fortran_kinds_v, only: int32, dp, stdout
   use openmp_mpi_m
   use timings_m
   use vmc_sampling_m
   use write_lines_m
   use quantum_monte_carlo_v,  only: bin_l, n_bin, n_trm, n_wlk, n_wlk_max, n_tot_wlk, &
   & qmc_mthd, var_dt, dt_e, dt_p, dt_d, arate, save_status, G_w_trgt, n_tot_smpl, n_tot_smpl_old, n_scr
   use molecular_system_v,     only: n_fe, molsys_prs, n_po, n_el
   use qdo_system_v,           only: qdosys_prs, n_qdo
   use acceptance_rate_m
   use monte_carlo_averages_m, only: bld_obsvec, obs_b, &
   & accu_qmcavg, write_qmcavg, n_obs, n_obs2, dealloc_qmcsmp
   use diffusion_monte_carlo_v, only: w
   use monte_carlo_checkpointfiles_m, only: save_checkpointfiles, clean_checkpointfiles
   implicit none
   public :: Init_VariationalMonteCarlo, Exec_VariationalMonteCarlo, &
   & Fnlz_VariationalMonteCarlo, Thrm_VariationalMonteCarlo
contains
   subroutine Init_VariationalMonteCarlo()
      allocate(w(1:n_wlk_max)) ; w(1:n_wlk) = 1.0_dp / G_w_trgt
      if ( n_wlk .lt. n_wlk_max ) w(n_wlk+1:n_wlk_max) = 0.0_dp
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of walkers", n_tot_wlk,var_name="n_tot_wlk")
      call write_variable_line(stdout,0,mpi_rank,2,"Number of independent walkers per MPI task", n_wlk,var_name="n_wlk")
      call write_variable_line(stdout,0,mpi_rank,2,"Monte Carlo sampling per core", n_tot_smpl,var_name="n_tot_smpl")
      call write_variable_line(stdout,0,mpi_rank,2,"Particle sweeps before derivatives computation", n_bra_vmc,var_name="n_bra_vmc")
      call write_variable_line(stdout,0,mpi_rank,2,"Number of moves before matrix recomputation", n_scr,var_name="n_scr")
   end subroutine Init_VariationalMonteCarlo
   subroutine Thrm_VariationalMonteCarlo( )
      integer(int32) :: i1, i2
      call time_compute(t_therm_vmc)
      call vmc_initialize_all_particles_and_walkers()
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","Variational Monte Carlo thermalization")
      call write_empty_line(stdout,0,mpi_rank)
      if ( mpi_rank.eq.0_int32 ) then
         write(*,'(3X,"          Block N°      Acc. rate. [%]      Time step [au]      ")')
         if (molsys_prs.and.qdosys_prs) then
            if (n_po.gt.0_int32.and.n_el.gt.0_int32) then
               write(*,'(3X,"           ele.     pos.       drd.      ele.      pos.      drd.    ")')
            elseif( n_el.gt.0_int32 ) then
               write(*,'(3X,"                        ele.     drd.       ele.      drd.      ")')
            else
               write(*,'(3X,"                        pos.     drd.       pos.      drd.      ")')
            endif
         else
            if ( molsys_prs ) then
               if (n_po.gt.0_int32.and.n_el.gt.0_int32) then
                  write(*,'(3X,"                        ele.     pos.       ele.      pos.      ")')
               endif
            endif
         endif
      endif 
      acc_rate_e = arate
      acc_rate_p = arate
      acc_rate_d = arate
      if (n_el.gt.0_int32)  acc_rate_e = 2.0
      if (n_po.gt.0_int32)  acc_rate_p = 2.0
      if (n_qdo.gt.0_int32) acc_rate_d = 2.0
      i2 = 1_int32
      do while ( i2.lt.n_trm )
         do i1 = 1_int32, bin_l
            call vmc_move_all_particles_and_walkers()
         enddo 
         call cmp_acc_rates()
         call upd_var_dt( prnt=.false. )
         if ( mpi_rank.eq.0_int32 ) then
            if (molsys_prs.and.qdosys_prs) then
               if (n_po.gt.0_int32.and.n_el.gt.0_int32) then
                  write(*,'(7X,I10,8X,F7.3,2X,F7.3,2X,F7.3,4X,E7.3,2X,E7.3,2X,E7.3)') &
                  & i2, acc_rate_e*100.0_dp, acc_rate_p*100.0_dp, acc_rate_d*100.0_dp, dt_e, dt_p, dt_d
               elseif( n_el.gt.0_int32 ) then
                  write(*,'(7X,I10,8X,F7.3,2X,F7.3,4X,F7.4,2X,F7.4)') i2, acc_rate_e*100.0_dp, acc_rate_d*100.0_dp, dt_e, dt_d
               else
                  write(*,'(7X,I10,8X,F7.3,2X,F7.3,4X,F7.4,2X,F7.4)') i2, acc_rate_p*100.0_dp, acc_rate_d*100.0_dp, dt_p, dt_d
               endif
            else
               if ( molsys_prs ) then
                  if (n_po.gt.0_int32.and.n_el.gt.0_int32) then
                     write(*,'(7X,I10,8X,F7.3,2X,F7.3,4X,F7.4,2X,F7.4)') i2, acc_rate_e*100.0_dp, acc_rate_p*100.0_dp, dt_e, dt_p
                  elseif (n_el.gt.0_int32) then
                     write(*,'(6X,I10,5X,F15.3,10X,E15.8)') i2, acc_rate_e * 100.0_dp, dt_e
                  else
                     write(*,'(6X,I10,5X,F15.3,10X,E15.8)') i2, acc_rate_p * 100.0_dp, dt_p
                  endif
               endif
               if (qdosys_prs) then
                  write(*,'(6X,I10,5X,F15.3,10X,E15.8)') i2, acc_rate_d * 100.0_dp, dt_d
               endif
            endif
         endif
         i2 = i2 + 1_int32
      enddo
      call write_empty_line(stdout,0,mpi_rank)
      if(n_el.gt.0_int32) &
      & call write_variable_line(stdout,0,mpi_rank,2,"Final Time step for electrons", dt_e, var_name="dt_e")
      if(n_po.gt.0_int32) &
      & call write_variable_line(stdout,0,mpi_rank,2,"Final Time step for electrons", dt_p, var_name="dt_p")
      if(n_qdo.gt.0_int32) &
      & call write_variable_line(stdout,0,mpi_rank,2,"Final Time step for electrons", dt_d, var_name="dt_d")
      call time_compute_difference(t_therm_vmc)
      call save_checkpointfiles( save_status_in = 1_int32, step_in=0_int32, create_backup=.true. )
   end subroutine Thrm_VariationalMonteCarlo
   subroutine Exec_VariationalMonteCarlo()
      integer(int32) :: i1, i2
      character(32) :: block_strng
      call vmc_initialize_all_particles_and_walkers()
      var_dt = .false.
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","Variational Monte Carlo sampling")
      call write_empty_line(stdout,0,mpi_rank)
      do i1 = 1_int32, n_bin
         obs_b(:,:) = 0.0_dp
         call time_compute(t_sampl_step)
         do i2 = 1_int32, bin_l
            call vmc_move_all_particles_and_walkers()
            call bld_obsvec( obs_b(:,i2) )
            call accu_qmcavg(obs_b(:,i2) )
         enddo 
         call time_compute_difference(t_sampl_step,time_cpu_tot=t_sampl_tot)
         write(block_strng,'("Block Num.",I11," sampled in")') i1
         call write_variable_line(stdout,0,mpi_rank,2,block_strng, t_sampl_step,units='sec.')
         call save_checkpointfiles( save_status_in = 2_int32, step_in=i1, create_backup=.false.  )
      enddo 
      call write_empty_line(stdout,0,mpi_rank)
      call cmp_acc_rates()
      call write_acc_rates()
   end subroutine Exec_VariationalMonteCarlo
   subroutine Fnlz_VariationalMonteCarlo()
      if (mpi_rank.eq.0_int32) call write_qmcavg( )
      call dealloc_qmcsmp( )
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","VARIATIONAL MONTE CARLO CALCULATION ENDED")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"c","All well, farewell! See you soon!!")
      call clean_checkpointfiles()
   end subroutine Fnlz_VariationalMonteCarlo
end module variational_monte_carlo_m
