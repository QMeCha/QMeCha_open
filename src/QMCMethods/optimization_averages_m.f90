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
module optimization_averages_m
   use fortran_kinds_v, only: dp, int32
   use merge_sort_m, only: mrg_srt_dbl
   use openmp_mpi_m
   use write_lines_m
   use timings_m
   use quantum_monte_carlo_v,      only: n_bin, n_wlk, bin_l, n_trm, n_bra, n_tot_wlk
   use correlated_sampling_m,      only: Save_CorrelatedSampling, crrsmp, W_avg, W_err, W2_avg, &
   & Ovrlap_avg, logshift
   use wavefunction_optimization_v
   use molecular_system_v,         only: molsys_prs, n_fe
   use qdo_system_v,               only: qdosys_prs, n_qdo
   use fermionic_config_c,         only: frm_cnf
   use drudonic_config_c,          only: drd_cnf
   use fermionic_sampling_m,       only: smpl_mtrpls_frm
   use drudonic_sampling_m,        only: smpl_mtrpls_drd
   use local_energy_m,             only: cmp_locene, e_l
   use fermionic_wavefunction_m,   only: wvfn, sym_wvfn_par
   use qdo_wavefunction_m,         only: wvfnq, sym_qdowvf_par
   use fermionic_wavefunction_v,   only: n_par_frm_wvfn_s, n_par_frm_wvfn
   use qdo_wavefunction_v,         only: n_par_drd_wvfn
   implicit none
   real(dp), external :: regularization_f
   public :: Simple_Sampling_Accumulation, Correlated_Sampling_Accumulation, &
   &Comp_Jackknife_Averages, Comp_Jackknife_Weighted_Averages, Comp_Jackknife_Covariances,&
   &Comp_MaxForce_Deviation
contains
   subroutine Simple_Sampling_Accumulation( )
      logical        :: move_acc
      integer(int32) :: iw
      integer(int32) :: i1, i2, i3, i4
      real(dp)       :: force_reg
      call time_compute(t_sampl_step)
!$omp parallel default(shared) private(iw,i1,i2,i3,i4,move_acc,force_reg)
!$omp do
      do iw = 1_int32, n_wlk
         e_all(:,iw)     = 0.0_dp
         O_all(:,:,iw)   = 0.0_dp
         if (molsys_prs) call wvfn%cmp(iw)
         if (qdosys_prs) call wvfnq(iw)%cmp(iw)
         do i1 = 1_int32, n_trm * n_bra
            if (molsys_prs) then
               do i2 = 1_int32, n_fe
                  call smpl_mtrpls_frm( iw, move_acc, count=.false. )
                  if (move_acc) then
                     call frm_cnf(iw)%upd( )
                     if (qdosys_prs) call drd_cnf(iw)%upd(iw)
                     call wvfn%upd( iw )
                     if (qdosys_prs) call wvfnq(iw)%upd(iw)
                  endif
               enddo ! i2 ( n_fe  )
            endif
            if (qdosys_prs) then
               do i2 = 1_int32, n_qdo
                  call smpl_mtrpls_drd( iw, move_acc, count=.false. )
                  if (move_acc) then
                     call drd_cnf(iw)%upd(iw)
                     call wvfnq(iw)%upd(iw)
                  endif
               enddo ! i2 ( n_qdo )
            endif
         enddo ! i3 n_trm * n_bra
         call cmp_locene(iw, K_cmp=.true.)
         do i1 = 1_int32, n_bin * bin_l
            do i2 = 1_int32, n_bra
               if (molsys_prs) then
                  do i3 = 1_int32, n_fe
                     call smpl_mtrpls_frm( iw, move_acc, count=.true. )
                     if (move_acc) then
                        call frm_cnf(iw)%upd( )
                        if (qdosys_prs) call drd_cnf(iw)%upd(iw)
                        call wvfn%upd( iw )
                        if (qdosys_prs) call wvfnq(iw)%upd(iw)
                     endif
                  enddo ! i3 ( n_fe  )
               endif
               if (qdosys_prs) then
                  do i3 = 1_int32, n_qdo
                     call smpl_mtrpls_drd( iw, move_acc, count=.true. )
                     if (move_acc) then
                        call drd_cnf(iw)%upd(iw)
                        call wvfnq(iw)%upd(iw)
                     endif
                  enddo ! i3 ( n_qdo )
               endif
            enddo ! i2 n_bra
            if ( corr_samp ) call Save_CorrelatedSampling( iw, i1 )
            call cmp_locene(iw, K_cmp=.true.)
            e_all(i1,iw) = e_l(iw)
            call cmp_V2_avg( iw, force_reg )
            force_reg = 1.0_dp/sqrt(force_reg)
            force_reg = regularization_f( force_reg, eps_force_reg )
            if (n_par_frm_wvfn.gt.0_int32) then
               call wvfn%cmp_dp(iw)
               wvfn%dp_ln_wvfn(:,iw)= wvfn%dp_ln_wvfn(:,iw) * force_reg
               call sym_wvfn_par( wvfn%dp_ln_wvfn(:,iw),  O_all(1:n_par_frm_wvfn_s,i1,iw) )
            endif
            if (n_par_drd_wvfn.gt.0_int32) then
               call wvfnq(iw)%cmp_dp(iw)
               wvfnq(iw)%dp_ln_wvfnq(:) = wvfnq(iw)%dp_ln_wvfnq(:) * force_reg
               call sym_qdowvf_par(wvfnq(iw)%dp_ln_wvfnq(:), O_all(n_par_frm_wvfn_s+1:n_par_opt_s,i1,iw))
            end if
         enddo ! i1 ( n_bin * bin_l )
      enddo ! iw ( n_wlk )
!$omp end do
!$omp end parallel
      call time_compute_difference(t_sampl_step,time_cpu_tot=t_sampl_tot)
   end subroutine Simple_Sampling_Accumulation
   subroutine Correlated_Sampling_Accumulation( )
      real(dp)       :: r, spsi, logpsi, force_reg
      integer(int32) :: iw
      integer(int32) :: i1
      call time_compute(t_sampl_step)
      logshift = 0.0_dp
!$omp parallel default(shared) private(iw,logpsi) reduction(+:logshift)
!$omp do
      do iw = 1_int32, n_wlk
         if (molsys_prs) then
            frm_cnf(iw)%r_fe(:,:) = crrsmp(iw)%r_fe_sav(:,:,1)
            call frm_cnf(iw)%cmp()
         endif
         if (qdosys_prs) then
            drd_cnf(iw)%r_drd(:,:) = crrsmp(iw)%r_drd_sav(:,:,1)
            call drd_cnf(iw)%cmp(iw)
         endif
         logpsi = 0.0_dp
         if (molsys_prs) then
            call wvfn%cmp( iw )
            logpsi = logpsi +  wvfn%logwvfn(iw)
         endif
         if (qdosys_prs) then
            call wvfnq(iw)%cmp( iw )
            logpsi = logpsi + wvfnq(iw)%logwvfnq
         endif
         logshift  = logshift + crrsmp(iw)%logpsi(1)-logpsi
      enddo ! iw ( n_wlk )
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE, logshift, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      logshift = logshift / dble(n_wlk*n_mpi_tasks)
      Ovrlap_avg   = 0.0_dp
!$omp parallel default(shared) private(iw,i1,r,logpsi,spsi,force_reg) reduction(+:Ovrlap_avg)
!$omp do
      do iw = 1_int32, n_wlk
         O_all(:,:,iw)   = 0.0_dp
         do i1 = 1_int32, n_bin * bin_l
            if (molsys_prs) then
               frm_cnf(iw)%r_fe(:,:) = crrsmp(iw)%r_fe_sav(:,:,i1)
               call frm_cnf(iw)%cmp()
            endif
            if (qdosys_prs) then
               drd_cnf(iw)%r_drd(:,:) = crrsmp(iw)%r_drd_sav(:,:,i1)
               call drd_cnf(iw)%cmp(iw)
            endif
            spsi   = 1.0_dp
            logpsi = 0.0_dp
            if (molsys_prs) then
               call wvfn%cmp( iw )
               spsi   = spsi*wvfn%swvfn(iw)
               logpsi = logpsi +  wvfn%logwvfn(iw)
            endif
            if (qdosys_prs) then
               call wvfnq(iw)%cmp( iw )
               spsi   = spsi*wvfnq(iw)%swvfnq
               logpsi = logpsi + wvfnq(iw)%logwvfnq
            endif
            r = dexp(logpsi-crrsmp(iw)%logpsi(i1)+logshift)
            Ovrlap_avg = Ovrlap_avg + crrsmp(iw)%spsi(i1) * spsi * r
            crrsmp(iw)%ratio2(i1) = r**2
            call cmp_locene(iw, K_cmp=.true.)
            call cmp_V2_avg( iw, force_reg )
            force_reg = 1.0_dp/sqrt(force_reg)
            force_reg = regularization_f( force_reg, eps_force_reg )
            e_all(i1,iw) = e_l(iw)
            if ( n_par_frm_wvfn.gt.0_int32 ) then
               call wvfn%cmp_dp(iw)
               wvfn%dp_ln_wvfn(:,iw) = wvfn%dp_ln_wvfn(:,iw) * force_reg
               call sym_wvfn_par( wvfn%dp_ln_wvfn(:,iw),  O_all(1:n_par_frm_wvfn_s,i1,iw) )
            endif
            if ( n_par_drd_wvfn.gt.0_int32 ) then
               call wvfnq(iw)%cmp_dp(iw)
               wvfnq(iw)%dp_ln_wvfnq(:) = wvfnq(iw)%dp_ln_wvfnq(:) * force_reg
               call sym_qdowvf_par(wvfnq(iw)%dp_ln_wvfnq(:), O_all(n_par_frm_wvfn_s+1:n_par_opt_s,i1,iw))
            endif
         enddo ! i1 ( n_bin )
      enddo ! iw ( n_wlk )
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE, Ovrlap_avg,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      call time_compute_difference(t_sampl_step,time_cpu_tot=t_resampl_tot)
   end subroutine Correlated_Sampling_Accumulation
   subroutine Comp_Jackknife_Averages( )
      integer(int32) :: iw, ib, il
      real(dp) :: factor
      e_tot = 0.0_dp ; O(:) = 0.0_dp ; Oe(:) = 0.0_dp
!$omp parallel default(shared) private(iw, ib) reduction(+:e_tot,O,Oe)
!$omp do
      do iw = 1_int32, n_wlk
         e_tot = e_tot + sum(e_all(:,iw))
         do ib = 1_int32, n_sampling_task
            O(:) = O(:)   + O_all(:,ib,iw)
            Oe(:) = Oe(:) + O_all(:,ib,iw) * e_all(ib,iw)
         enddo ! ib ( n_bin * bin_l )
      enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE, e_tot,     1_int32, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE, O(:),  n_par_opt_s, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE, Oe(:), n_par_opt_s, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      factor = n_tot_sampling - dble(bin_l)
      f(:) = 0.0_dp; e_err = 0.0_dp
!$omp parallel default(shared) private(iw, ib, il) reduction(+:f,e_err)
!$omp do
      do iw = 1_int32, n_wlk
         do ib = 1_int32, n_bin
            pip0(:,iw) = 0.0_dp
            pip1(:,iw) = 0.0_dp
            pip2(iw) = 0.0_dp
            do il = 1_int32,  bin_l
               pip0(:,iw) = pip0(:,iw) + O_all(:,bin_l*(ib-1_int32)+il,iw)
               pip1(:,iw) = pip1(:,iw) + O_all(:,bin_l*(ib-1_int32)+il,iw) * e_all(bin_l*(ib-1_int32)+il,iw)
               pip2(iw)   = pip2(iw) + e_all(bin_l*(ib-1_int32)+il,iw)
            enddo ! il ( bin_l )
            pip0(:,iw) = (O(:) - pip0(:,iw)) / factor
            pip1(:,iw) = (Oe(:) - pip1(:,iw)) / factor
            pip2(iw)   = (e_tot - pip2(iw)) / factor
            e_err      =  e_err + (pip2(iw) - e_tot/n_tot_sampling )**2
            f_b(:,ib,iw) = -2.0_dp *( pip1(:,iw) - pip0(:,iw) * pip2(iw) )
            f(:) = f(:) + f_b(:,ib,iw)
         enddo ! ib ( n_bin )
      enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE,  f(:), n_par_opt_s, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE, e_err,     1_int32, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      e_tot = e_tot / n_tot_sampling
      O(:)  = O(:)  / n_tot_sampling
      Oe(:) = Oe(:) / n_tot_sampling
      f(:)  = f(:)  / n_tot_bins
      e_err = e_err * n_red_bins / n_tot_bins
      e_err = sqrt(e_err)
      f_err(:) = 0.0_dp
!$omp parallel default(shared) private(iw,ib) reduction(+:f_err)
!$omp do
      do iw = 1_int32, n_wlk
         do ib = 1_int32, n_bin
            f_err(:) = f_err(:) + (f_b(:,ib,iw) - f(:))**2
         enddo
      enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE,  f_err(:),  n_par_opt_s, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      factor = n_red_bins / n_tot_bins
      f_err(:) = sqrt(f_err(:)*factor)    ! Standard Error of forces
      f(:) = 0.0_dp ; e_var = 0.0_dp
      factor = n_tot_sampling/ n_red_sampling
!$omp parallel default(shared) private(iw, ib) reduction(+:f,e_var)
!$omp do
      do iw = 1_int32, n_wlk
         do ib = 1_int32, n_sampling_task
            pip1(:,iw) = (Oe(:) - O_all(:,ib,iw) * e_all(ib,iw)/n_tot_sampling) * factor
            pip2(iw)   = (e_tot - e_all(ib,iw)/n_tot_sampling) * factor
            O_all(:,ib,iw) = (O(:) - O_all(:,ib,iw)/n_tot_sampling ) * factor
            e_var     =  e_var + (pip2(iw)- e_tot)**2
            f_b(:,ib,iw) = -2.0_dp *( pip1(:,iw) - O_all(:,ib,iw) * pip2(iw) )
            f(:) = f(:) + f_b(:,ib,iw)
         enddo ! ib ( n_sampling_task )
      enddo ! iw ( n_wlk )
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE,  f(:), n_par_opt_s, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE, e_var,     1_int32, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      e_var = e_var*n_red_sampling
      f(:) = f(:) / n_tot_sampling
      if (trim(opt_mthd).ne.'lma') then
!$omp parallel default(shared) private(iw, ib)
!$omp do
         do iw = 1_int32, n_wlk
            do ib = 1_int32, n_sampling_task
               O_all(:,ib,iw) = ( O_all(:,ib,iw) - O(:) )
               f_b(:,ib,iw)   = ( f_b(:,ib,iw)   - f(:) )
            enddo ! ib ( n_sampling_task )
         enddo ! iw ( n_wlk )
!$omp end do
!$omp end parallel
      else
         factor = n_tot_sampling/n_red_sampling
!$omp parallel default(shared) private(iw, ib)
!$omp do
         do iw = 1_int32, n_wlk
            do ib = 1_int32, n_sampling_task
               O_all(:,ib,iw) = ( O_all(:,ib,iw) - O(:) )
               f_b(:,ib,iw)   = factor*f(:)-f_b(:,ib,iw) ! Divided by N -1
            enddo ! ib ( n_sampling_task )
         enddo ! iw ( n_wlk )
!$omp end do
!$omp end parallel
      endif
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Energy", e_tot, e_err,units="Eh")
      call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Variance",e_var, e_err/sqrt( 2.0_dp * n_red_sampling),units="Eh^2")
      f(:) = -2.0_dp* n_tot_sampling * ( Oe(:) - O(:)*e_tot) - n_red_sampling * f(:)
   end subroutine Comp_Jackknife_Averages
   subroutine Comp_Jackknife_Weighted_Averages( )
      integer(int32) :: iw, ib, il, id
      real(dp) :: factor
      W_avg = 0.0_dp
!$omp parallel default(shared) private(iw, ib) reduction(+:W_avg)
!$omp do
      do iw = 1_int32, n_wlk
         W_avg = W_avg + sum(crrsmp(iw)%ratio2(:))
      enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE, W_avg ,      1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      Ovrlap_avg = Ovrlap_avg / sqrt(W_avg*n_tot_sampling) ! <R> / sqrt(<R^2>)
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Wave functions overlap",Ovrlap_avg )
      e_tot = 0.0_dp ; ; O(:) = 0.0_dp ; Oe(:) = 0.0_dp
!$omp parallel default(shared) private(iw, ib) reduction(+:e_tot,O,Oe)
!$omp do
      do iw = 1_int32, n_wlk
         do ib = 1_int32, n_sampling_task
            crrsmp(iw)%ratio2(ib) = crrsmp(iw)%ratio2(ib) / W_avg
            e_tot = e_tot + e_all(ib,iw) * crrsmp(iw)%ratio2(ib)
            O(:) = O(:)   + O_all(:,ib,iw) * crrsmp(iw)%ratio2(ib)
            Oe(:) = Oe(:) + O_all(:,ib,iw) * e_all(ib,iw)* crrsmp(iw)%ratio2(ib)
         enddo ! ib ( n_bin * bin_l )
      enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE, e_tot,            1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE, O(:) ,  n_par_opt_s, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE, Oe(:) , n_par_opt_s, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      W_avg = 1.0_dp
      f(:) = 0.0_dp; e_err = 0.0_dp; W2_avg = 0.0_dp
!$omp parallel default(shared) private(iw, ib, il, id) reduction(+:f,e_err,W2_avg)
!$omp do
      do iw = 1_int32, n_wlk
         do ib = 1_int32, n_bin
            pip0(:,iw)  = 0.0_dp
            pip1(:,iw)  = 0.0_dp
            pip2(iw)    = 0.0_dp
            pip3(ib,iw) = 0.0_dp
            do il = 1_int32,  bin_l
               id = bin_l*(ib-1_int32)+il
               pip0(:,iw)  = pip0(:,iw) + O_all(:,id,iw) * crrsmp(iw)%ratio2(id)
               pip1(:,iw)  = pip1(:,iw) + O_all(:,id,iw) * e_all(id,iw) * crrsmp(iw)%ratio2(id)
               pip2(iw)    = pip2(iw) + e_all(id,iw) * crrsmp(iw)%ratio2(id)
               pip3(ib,iw) = pip3(ib,iw) + crrsmp(iw)%ratio2(id)
            enddo ! il ( bin_l )
            W2_avg = W2_avg + (pip3(ib,iw))**2
            pip3(ib,iw) = W_avg - pip3(ib,iw)
            pip0(:,iw) = (W_avg*O(:) - pip0(:,iw)) / pip3(ib,iw)
            pip1(:,iw) = (W_avg*Oe(:) - pip1(:,iw)) / pip3(ib,iw)
            pip2(iw)   = (W_avg*e_tot - pip2(iw)) / pip3(ib,iw)
            e_err      =  e_err + (pip2(iw) - e_tot )**2 *pip3(ib,iw)**2/ (W_avg-pip3(ib,iw))
            f_b(:,ib,iw) = -2.0_dp *( pip1(:,iw) - pip0(:,iw) * pip2(iw) )
            f(:) = f(:) + f_b(:,ib,iw) * pip3(ib,iw)/ n_red_bins
         enddo ! ib ( n_bin )
      enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE,   f(:), n_par_opt_s, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE,  e_err,     1_int32, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE, W2_avg,     1_int32, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      e_tot = e_tot / W_avg
      O(:)  = O(:)  / W_avg
      Oe(:) = Oe(:) / W_avg
      f(:) = f(:) / W_avg
      e_err = sqrt(e_err / (W_avg**2 -W2_avg) * W2_avg / W_avg  )
      f_err(:) = 0.0_dp
!$omp parallel default(shared) private(iw,ib) reduction(+:f_err)
!$omp do
      do iw = 1_int32, n_wlk
         do ib = 1_int32, n_bin
            f_err(:) = f_err(:) + (f_b(:,ib,iw) - f(:))**2 *pip3(ib,iw)**2/ (W_avg-pip3(ib,iw))
         enddo
      enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE,  f_err(:),  n_par_opt_s, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      f_err(:) = f_err(:) / (W_avg**2 -W2_avg) * W2_avg / W_avg ! Leslie Kish factor N_eff = W_avg**2/W2_avg
      f_err(:) = sqrt(f_err(:))
      f(:) = 0.0_dp ; e_var = 0.0_dp; W2_avg = 0.0_dp
!$omp parallel default(shared) private(iw, ib, factor) reduction(+:f,e_var,W2_avg)
!$omp do
      do iw = 1_int32, n_wlk
         do ib = 1_int32, n_sampling_task
            W2_avg = W2_avg + crrsmp(iw)%ratio2(ib)**2
            factor = W_avg - crrsmp(iw)%ratio2(ib)
            pip1(:,iw) = (W_avg*Oe(:) - crrsmp(iw)%ratio2(ib)*O_all(:,ib,iw) * e_all(ib,iw) )/factor
            O_all(:,ib,iw) = (W_avg*O(:) - crrsmp(iw)%ratio2(ib)*O_all(:,ib,iw))/factor
            pip2(iw)   = (W_avg*e_tot - crrsmp(iw)%ratio2(ib)* e_all(ib,iw))/factor
            f_b(:,ib,iw) = -2.0_dp *( pip1(:,iw) - O_all(:,ib,iw) * pip2(iw) )
            e_var     =  e_var + (pip2(iw)- e_tot)**2 *factor**2/crrsmp(iw)%ratio2(ib)
            f(:) = f(:) + f_b(:,ib,iw) * factor / n_red_sampling
         enddo ! ib ( n_bin )
      enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE,  f(:), n_par_opt_s, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE, e_var,     1_int32, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE, W2_avg,    1_int32, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      f(:)  = f(:) / W_avg
      e_var = e_var *W_avg / (W_avg**2 -W2_avg)
      if (trim(opt_mthd).ne.'lma') then
!$omp parallel default(shared) private(iw, ib)
!$omp do
         do iw = 1_int32, n_wlk
            do ib = 1_int32, n_sampling_task
               O_all(:,ib,iw) = ( O_all(:,ib,iw) - O(:) )
               f_b(:,ib,iw)   = ( f_b(:,ib,iw)   - f(:) )
            enddo ! ib ( n_sampling_task )
         enddo ! iw ( n_wlk )
!$omp end do
!$omp end parallel
      else
!$omp parallel default(shared) private(iw, ib)
!$omp do
         do iw = 1_int32, n_wlk
            do ib = 1_int32, n_sampling_task
               pip2(iw) = W_avg/ (W_avg - crrsmp(iw)%ratio2(ib))
               O_all(:,ib,iw) = ( O_all(:,ib,iw) - O(:) )
               f_b(:,ib,iw)   = pip2(iw)*f(:)-f_b(:,ib,iw) ! multiplied by w_i/(W -w_i)
            enddo ! ib ( n_sampling_task )
         enddo ! iw ( n_wlk )
!$omp end do
!$omp end parallel
      endif
      call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Energy", e_tot, e_err,units="Eh")
      call write_variable_error_dble_line(stdout,0,mpi_rank,2,"Variance",e_var, e_err/sqrt( 2.0_dp *(W_avg**2 -W2_avg)),units="Eh^2")
      f(:) = - 2.0_dp *( W_avg**2 / W2_avg ) * ( Oe(:) - O(:)* e_tot) - ( W_avg**2 / W2_avg -1.0_dp ) * f(:)
   end subroutine Comp_Jackknife_Weighted_Averages
   subroutine Comp_Jackknife_Covariances( weighted_covariance )
      logical, intent(in) :: weighted_covariance
      integer(int32) :: ip, ib, iw
      real(dp) :: covariance_prefactor2
      if (.not.weighted_covariance) then
         covariance_prefactor = sqrt(n_red_sampling)
         if (trim(opt_mthd).eq.'snr'.or.trim(opt_mthd).eq.'lma') then 
            covariance_prefactor2 = covariance_prefactor / (2.0_dp *sqrt(e_var) )
         else 
            covariance_prefactor2 = covariance_prefactor
         endif
!$omp parallel default(shared) private(iw,ip,ib)
!$omp do
         do iw = 1_int32, n_wlk
            do ib = 1_int32, n_sampling_task
               do ip = 1_int32, n_par_sgn
                  f_b(ip,ib,iw) = f_b(f_indx(ip),ib,iw) * covariance_prefactor2
                  O_all(ip,ib,iw) = O_all(f_indx(ip),ib,iw) * covariance_prefactor
               enddo
            enddo
            O_all(n_par_sgn+1:n_par_opt_s,:,iw) = 0.0_dp
            f_b(n_par_sgn+1:n_par_opt_s,:,iw) = 0.0_dp
         enddo
!$omp end do
!$omp end parallel
      else
         if (trim(opt_mthd).eq.'lma') then ! For LM Algorithm
            covariance_prefactor = sqrt(W_avg / (W_avg**2 -W2_avg))
            !$omp parallel default(shared) private(iw,ip,ib)
            !$omp do
            do iw = 1_int32, n_wlk
               do ib = 1_int32, n_sampling_task
                  pip2(iw) = (W_avg - crrsmp(iw)%ratio2(ib)) / sqrt(crrsmp(iw)%ratio2(ib)) * covariance_prefactor
                  do ip = 1_int32, n_par_sgn
                     f_b(ip,ib,iw) = f_b(f_indx(ip),ib,iw) * pip2(iw) * sqrt(W_avg) / (2.0_dp *sqrt(e_var) )
                     O_all(ip,ib,iw) = O_all(f_indx(ip),ib,iw) * pip2(iw)
                  enddo
               enddo
               O_all(n_par_sgn+1:n_par_opt_s,:,iw) = 0.0_dp
               f_b(n_par_sgn+1:n_par_opt_s,:,iw) = 0.0_dp
            enddo
            !$omp end do
            !$omp end parallel
         else 
            covariance_prefactor = sqrt(W_avg / (W_avg**2 -W2_avg))
            !$omp parallel default(shared) private(iw,ip,ib)
            !$omp do
            do iw = 1_int32, n_wlk
               do ib = 1_int32, n_sampling_task
                  pip2(iw) = (W_avg - crrsmp(iw)%ratio2(ib)) / sqrt(crrsmp(iw)%ratio2(ib)) * covariance_prefactor
                  do ip = 1_int32, n_par_sgn
                     if ( trim(opt_mthd).eq.'snr') then
                        f_b(ip,ib,iw) = f_b(f_indx(ip),ib,iw) *pip2(iw)/ (2.0_dp *sqrt(e_var) )
                     else 
                        f_b(ip,ib,iw) = f_b(f_indx(ip),ib,iw) *pip2(iw)
                     endif
                     O_all(ip,ib,iw) = O_all(f_indx(ip),ib,iw) *pip2(iw)
                  enddo
               enddo
               O_all(n_par_sgn+1:n_par_opt_s,:,iw) = 0.0_dp
               f_b(n_par_sgn+1:n_par_opt_s,:,iw) = 0.0_dp
            enddo
            !$omp end do
            !$omp end parallel
         endif
      endif
      if (trim(opt_mthd).ne.'ams'.and.trim(lin_solv_mthd).eq.'lu') then
         f_cov = 0.0_dp
         if (trim(opt_mthd).eq.'snr'.or.trim(opt_mthd).eq.'lma') then
            call dsyrk( 'L', 'N', n_par_opt_s, n_sampling_task, 1.0_dp, f_b(1:n_par_opt_s,:,1), n_par_opt_s, &
            & 0.0_dp, f_cov(1:n_par_opt_s,1:n_par_opt_s), n_par_opt_s )
            if (n_wlk.gt.1_int32) then
               do iw = 2_int32, n_wlk
                  call dsyrk( 'L', 'N', n_par_opt_s, n_sampling_task, 1.0_dp, f_b(1:n_par_opt_s,:,iw), n_par_opt_s, &
                  & 1.0_dp, f_cov(1:n_par_opt_s,1:n_par_opt_s), n_par_opt_s )
               enddo
            endif
         else
            call dsyrk( 'L', 'N', n_par_opt_s, n_sampling_task, 1.0_dp, O_all(1:n_par_opt_s,:,1), n_par_opt_s, &
            & 0.0_dp, f_cov(1:n_par_opt_s,1:n_par_opt_s), n_par_opt_s )
            if (n_wlk.gt.1_int32) then
               do iw = 2_int32, n_wlk
                  call dsyrk( 'L', 'N', n_par_opt_s, n_sampling_task, 1.0_dp, O_all(1:n_par_opt_s,:,iw), n_par_opt_s, &
                  & 1.0_dp, f_cov(1:n_par_opt_s,1:n_par_opt_s), n_par_opt_s )
               enddo
            endif
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         if(mpi_rank.eq.0_int32) then
            call mpi_reduce(MPI_IN_PLACE, f_cov(:,:) , n_par_opt_s**2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpierr )
         else
            call mpi_reduce( f_cov(:,:), f_cov(:,:), n_par_opt_s**2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpierr )
         endif
#endif
      else
         f_err = 0.0_dp
         if (trim(opt_mthd).eq.'snr'.or.trim(opt_mthd).eq.'lma') then
!$omp parallel default(shared) private(iw,ip) reduction(+:f_err)
!$omp do
            do iw = 1_int32, n_wlk
               do ip = 1_int32, n_par_sgn
                  f_err(ip) = f_err(ip) + dot_product(f_b(ip,:,iw),f_b(ip,:,iw) )
               enddo
            enddo
!$omp end do
!$omp end parallel
         else
!$omp parallel default(shared) private(iw,ip) reduction(+:f_err)
!$omp do
            do iw = 1_int32, n_wlk
               do ip = 1_int32, n_par_sgn
                  f_err(ip) = f_err(ip) + dot_product(O_all(ip,:,iw),O_all(ip,:,iw) )
               enddo
            enddo
!$omp end do
!$omp end parallel
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_allreduce(MPI_IN_PLACE,  f_err(1:n_par_sgn), n_par_sgn, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      endif
   end subroutine Comp_Jackknife_Covariances
   subroutine Comp_MaxForce_Deviation( )
      integer(int32) :: ip, ib
      real(dp) :: alpha
      n_par_sgn = 0_int32
      alpha = 1.0_dp / 0.75_dp
      f_indx = 0_int32
      if (force_cut.eq.0.0_dp) then
         n_par_sgn = n_par_opt_s
         do ib = 1_int32, n_par_opt_s
            f_indx(ib) = ib
         enddo
      else
         do while (n_par_sgn.lt.1_int32)
            alpha = alpha * 0.75_dp
            n_par_sgn = 0_int32
            f_indx = 0_int32
            ib = 1_int32; ip = n_par_opt_s
            do while(ib.le.ip)
               if ( abs(f(ib)).gt.force_cut .and. abs(f(ib)).gt.alpha*f_err(ib) ) then
                  f_indx(ib) = ib
                  ib = ib + 1_int32
                  n_par_sgn = n_par_sgn + 1_int32
               else
                  if ( abs(f(ip)).gt.force_cut .and. abs(f(ip)).gt.alpha*f_err(ip)) then
                     f_indx(ib) = ip
                     ib = ib + 1_int32
                     n_par_sgn = n_par_sgn + 1_int32
                  endif
                  ip = ip - 1_int32
               endif
            enddo
         enddo
         do ib = 1_int32, n_par_sgn
            if (f_indx(ib).ne.ib) then
               f(ib) = f(f_indx(ib))
               f_err(ib) = f_err(f_indx(ib))
               f(f_indx(ib)) = 0.0_dp
               f_err(f_indx(ib)) = 0.0_dp
            endif
         enddo
      endif
      if (n_par_sgn.gt.0.0_dp) then
         dev_max_f = maxval(abs(f(1:n_par_sgn)/f_err(1:n_par_sgn)))
      else 
         dev_max_f = 0.0_dp
      endif
      call write_variable_line(stdout,0,mpi_rank,2,"Number of signal active parameters",n_par_sgn )
      call write_variable_line(stdout,0,mpi_rank,2,"Maximum force deviation",dev_max_f )
   end subroutine Comp_MaxForce_Deviation
end module optimization_averages_m
