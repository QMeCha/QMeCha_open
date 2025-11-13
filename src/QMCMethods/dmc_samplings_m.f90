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
module dmc_samplings_m
   use fortran_kinds_v,          only: int32, dp
   use openmp_mpi_m
   use qdo_system_v,             only: qdosys_prs, n_qdo
   use molecular_system_v,       only: molsys_prs, n_fe
   use quantum_monte_carlo_v,    only: n_wlk, dt
   use diffusion_monte_carlo_v,  only: dt_eff, dr2, pdr2, e_l_old, w, wlk_life, &
   & max_wlk_life, e_l_pp_old, restart_walker, e_best, e_best_var, w_max, log_w_max,&
   & encrr_type, V2, V2_pp, e_r, g_e, fe_i_at_old
   use local_energy_m,           only: e_l, e_l_pp, cmp_locene
   use fermionic_wavefunction_m, only: wvfn
   use qdo_wavefunction_m,       only: wvfnq
   use fermionic_config_c,       only: frm_cnf
   use drudonic_config_c,        only: drd_cnf
   use kinetic_energy_m,         only: upd_kinene, upd_kinene_emb
   use electronic_properties_m,  only: cmp_mlt_mmn, elec_prp
   use drudonic_sampling_m,      only: smpl_lngvin_drd
   use fermionic_sampling_m,     only: smpl_lngvin_frm
   implicit none
   public  :: dmc_move_all_particles_and_walkers, dmc_initialize_all_particles_and_walkers
   private :: dmc_move_all_frm, dmc_move_all_drd
contains
   subroutine dmc_move_all_frm( iw )
      integer(int32), intent(in) :: iw
      logical        :: move_acc, wlk_moved
      integer(int32) :: i1
      wlk_moved = .false.
      do i1 = 1_int32, n_fe
         call smpl_lngvin_frm( iw, move_acc, .true., dr2=dr2(iw), pdr2=pdr2(iw) )
         if( move_acc ) then
            wlk_moved = .true.
            call frm_cnf(iw)%upd( )
            if (qdosys_prs) call drd_cnf(iw)%upd( iw )
            call wvfn%upd( iw )
            if (qdosys_prs) call wvfnq(iw)%upd( iw )
            call upd_kinene(iw)
            if (elec_prp) call cmp_mlt_mmn(iw)
            e_l_old(iw) = e_l(iw)
            e_l_pp_old(:,iw) = e_l_pp(:,iw)
         endif
      enddo 
      if ( .not.wlk_moved ) then
         wlk_life(iw) = wlk_life(iw) + 1_int32
      else
         wlk_life(iw) = 1_int32
      endif
   end subroutine dmc_move_all_frm
   subroutine dmc_move_all_drd( iw )
      integer(int32), intent(in) :: iw
      logical        :: move_acc, wlk_moved
      integer(int32) :: i1
      wlk_moved = .false.
      do i1 = 1_int32, n_qdo
         call smpl_lngvin_drd( iw, move_acc, .true., dr2=dr2(iw), pdr2=pdr2(iw) )
         if( move_acc ) then
            wlk_moved = .true.
            call drd_cnf(iw)%upd( iw )
            call wvfnq(iw)%upd( iw )
            call upd_kinene_emb(iw)
            if (elec_prp) call cmp_mlt_mmn(iw)
         endif
      enddo
      if ( .not.wlk_moved ) then
         wlk_life(iw) = wlk_life(iw) + 1_int32
      else
         wlk_life(iw) = 1_int32
      endif
   end subroutine dmc_move_all_drd
   subroutine dmc_move_all_particles_and_walkers()
      integer(int32) :: iw
!$omp parallel default(shared) private(iw)
!$omp do
      do iw = 1_int32, n_wlk
         pdr2(iw) = 0.0_dp ; dr2(iw) = 0.0_dp
         if (molsys_prs) call dmc_move_all_frm( iw )
         if (qdosys_prs) call dmc_move_all_drd( iw )
         if (encrr_type.lt.10_int32) then
            dt_eff(iw) = pdr2(iw) / dr2(iw)
         else 
            dt_eff(iw) = dt 
         endif
         call cmp_locene(iw, K_cmp=.false.)
         call cmp_dmc_wght( iw, e_l(iw), e_l_old(iw), e_l_pp(:,iw), e_l_pp_old(:,iw), w(iw) )
         if ( max_wlk_life.gt.1_int32.and.wlk_life(iw).ge.max_wlk_life )  then
            w(iw) = min( 0.5_dp, w(iw) )
         endif
         if (molsys_prs .and. (encrr_type.ge.20_int32)) &
         & fe_i_at_old(1:n_fe,iw) = frm_cnf(iw)%fe_i_at(1:n_fe)
         e_l_old(iw) = e_l(iw)
         e_l_pp_old(:,iw) = e_l_pp(:,iw)
      enddo
!$omp end do
!$omp end parallel
   end subroutine dmc_move_all_particles_and_walkers
   subroutine dmc_initialize_all_particles_and_walkers(  )
      integer(int32) :: iw
!$omp parallel default(shared) private(iw)
!$omp do
      do iw = 1_int32, n_wlk
         if ( restart_walker(iw) ) then
            if (molsys_prs) call frm_cnf(iw)%cmp( )
            if (qdosys_prs) call drd_cnf(iw)%cmp( iw )
            if (molsys_prs) call wvfn%cmp( iw )
            if (qdosys_prs) call wvfnq(iw)%cmp( iw )
            call cmp_locene(iw, K_cmp=.true.)
            select case(encrr_type)
             case(3)
               call cmp_V2_V2_bar_ratio(iw, V2(iw))
             case(4,5)
               call cmp_V2_avg(iw, V2(iw))
             case(13)
               call cmp_f_UNR93_vec(iw, V2_pp(:,iw))
             case(14,24)
               call cmp_V2_vec(iw, V2_pp(:,iw))
             case(15,25)
               call cmp_F_vec(iw, dt, V2_pp(:,iw))
             case default
            end select
            e_l_old(iw) = e_l(iw)
            if (molsys_prs .and. (encrr_type.ge.20_int32)) &
            & fe_i_at_old(1:n_fe,iw) = frm_cnf(iw)%fe_i_at(1:n_fe)
            e_l_pp_old(:,iw) = e_l_pp(:,iw)
            if (elec_prp) call cmp_mlt_mmn(iw)
            restart_walker(iw) = .false.
            wlk_life(iw) = 1_int32
         endif
      enddo
!$omp end do
!$omp end parallel
   end subroutine dmc_initialize_all_particles_and_walkers
end module dmc_samplings_m
