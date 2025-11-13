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
module vmc_sampling_m
   use fortran_kinds_v,          only: int32, dp
   use openmp_mpi_m
   use quantum_monte_carlo_v,    only: n_wlk, smp_type, n_bra_vmc
   use fermionic_wavefunction_m, only: wvfn
   use qdo_wavefunction_m, only: wvfnq
   use fermionic_config_c, only: frm_cnf
   use drudonic_config_c, only: drd_cnf
   use local_energy_m, only: cmp_locene, e_l
   use kinetic_energy_m, only: upd_kinene, upd_kinene_emb
   use qdo_system_v, only: qdosys_prs, n_qdo
   use molecular_system_v, only: molsys_prs, n_fe
   use electronic_properties_m, only: cmp_mlt_mmn, elec_prp
   use fermionic_sampling_m, only: smpl_lngvin_frm, smpl_mtrpls_frm
   use drudonic_sampling_m, only: smpl_lngvin_drd, smpl_mtrpls_drd
   implicit none
   public  :: vmc_initialize_all_particles_and_walkers, vmc_move_all_particles_and_walkers
   private :: vmc_move_all_frm, vmc_move_all_drd
contains
   subroutine vmc_move_all_frm( iw )
      integer(int32), intent(in) :: iw
      logical        :: move_acc
      integer(int32) :: i1
      do i1 = 1_int32, n_fe
         if ( smp_type.eq.10 ) then
            call smpl_lngvin_frm( iw, move_acc, count=.true. )
         else
            call smpl_mtrpls_frm( iw, move_acc, count=.true. )
         endif
         if( move_acc ) then
            call frm_cnf(iw)%upd( )
            if (qdosys_prs) call drd_cnf(iw)%upd( iw )
            call wvfn%upd( iw )
            if (qdosys_prs) call wvfnq(iw)%upd( iw )
            if ( smp_type.eq.10 ) call upd_kinene(iw)
         endif
      enddo
   end subroutine vmc_move_all_frm
   subroutine vmc_move_all_drd( iw )
      integer(int32), intent(in) :: iw
      logical        :: move_acc
      integer(int32) :: i1
      do i1 = 1_int32, n_qdo
         if (smp_type.eq.10 ) then
            call smpl_lngvin_drd( iw, move_acc, count=.true. )
         else
            call smpl_mtrpls_drd( iw, move_acc, count=.true. )
         endif
         if( move_acc ) then
            call drd_cnf(iw)%upd( iw )
            call wvfnq(iw)%upd( iw )
            if ( smp_type.eq.10  ) call upd_kinene_emb(iw)
         endif
      enddo
   end subroutine vmc_move_all_drd
   subroutine vmc_initialize_all_particles_and_walkers()
      integer(int32) :: iw
!$omp parallel do default(shared) private(iw)
      do iw = 1_int32, n_wlk
         if (molsys_prs) call frm_cnf(iw)%cmp( )
         if (qdosys_prs) call drd_cnf(iw)%cmp( iw )
         if (molsys_prs) call wvfn%cmp( iw )
         if (qdosys_prs) call wvfnq(iw)%cmp( iw )
         call cmp_locene(iw, K_cmp=.true.)
         if (elec_prp) call cmp_mlt_mmn(iw)
      enddo
!$omp end parallel do
   end subroutine vmc_initialize_all_particles_and_walkers
   subroutine vmc_move_all_particles_and_walkers()
      integer(int32) :: iw, i1
!$omp parallel do default(shared) private(iw,i1)
      do iw = 1_int32, n_wlk
         do i1 = 1_int32, n_bra_vmc
            if (molsys_prs) call vmc_move_all_frm( iw )
            if (qdosys_prs) call vmc_move_all_drd( iw )
         enddo
         if ( smp_type.eq.10 ) then
            call cmp_locene(iw, K_cmp=.false.)
         else
            call cmp_locene(iw, K_cmp=.true.)
         endif
         if (elec_prp) call cmp_mlt_mmn(iw)
      enddo
!$omp end parallel do
   end subroutine vmc_move_all_particles_and_walkers
end module vmc_sampling_m
