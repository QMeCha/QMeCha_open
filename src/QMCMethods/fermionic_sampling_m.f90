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
module fermionic_sampling_m
   use fortran_kinds_v,          only: int32, dp
   use physical_constants_v,     only: pi, twopi
   use openmp_mpi_m
   use quantum_monte_carlo_v,    only: dt_e, dt_p, var_dt, spc_map, smp_type, &
   & nuc_corr, a_drft, two_step_sampling
   use acceptance_rate_m
   use molecular_system_v,       only: n_po, n_el, m_e, m_p, r_at
   use fermionic_config_c,       only: frm_cnf
   use fermionic_wavefunction_m, only: wvfn
   use qdo_wavefunction_m,       only: wvfnq
   use drudonic_config_c,        only: drd_cnf
   use qdo_system_v,             only: qdosys_prs
   use mersenne_twister19937_m,  only: rndn
   use jstqepar_var,             only: jstqe_prs
   use jstqe_mod,                only: jstqe_fct
   implicit none
   public :: smpl_mtrpls_frm, smpl_lngvin_frm
contains
   subroutine smpl_mtrpls_frm( iw, move_acc, count )
      integer(int32), intent(in)  :: iw
      logical,        intent(out) :: move_acc
      logical,        intent(in)  :: count
      real(dp)               :: g, g_fq
      integer(int32)         :: i_dir
      real(dp), dimension(3) :: dr
      real(dp)               :: dt_i, dt_f, dexp_fun, T_ratio
      call frm_cnf(iw)%chs( iw )
      if (frm_cnf(iw)%i_fe.le.n_el) then
         dt_i = dt_e / m_e
         if (spc_map ) then
            dt_i = dt_i * frm_cnf(iw)%sm_el(frm_cnf(iw)%i_fe)
         endif
      else
         dt_i = dt_p / m_p
      endif
      select case( smp_type )
       case(0) ! Simple sampling 1d
         dr(:) = 0.0_dp
         i_dir = int( 3.0_dp * rndn(iw)%rndo() + 1.0_dp )
         dr(i_dir) = dt_i * (2.0_dp * rndn(iw)%rndo() - 1.0_dp)
       case(1) ! Gaussian move 1d
         dr(:) = 0.0_dp
         i_dir = int( 3.0_dp * rndn(iw)%rndo() + 1.0_dp )
         dr(i_dir) = sqrt(dt_i) * rndn(iw)%rndbm()
       case(2) ! Simple sampling 3d
         dr(1) = dt_i * (2.0_dp * rndn(iw)%rndo() - 1.0_dp)
         dr(2) = dt_i * (2.0_dp * rndn(iw)%rndo() - 1.0_dp)
         dr(3) = dt_i * (2.0_dp * rndn(iw)%rndo() - 1.0_dp)
       case(3) ! Gaussian move 3d
         dr(1) = sqrt(dt_i) * rndn(iw)%rndbm()
         dr(2) = sqrt(dt_i) * rndn(iw)%rndbm()
         dr(3) = sqrt(dt_i) * rndn(iw)%rndbm()
      end select
      frm_cnf(iw)%r_fe_new(:) = frm_cnf(iw)%r_fe(:,frm_cnf(iw)%i_fe) + dr(1:3)
      call frm_cnf(iw)%new( )
      if (frm_cnf(iw)%i_fe.le.n_el) then
         dt_f = dt_e / m_e
         if ( spc_map ) then
            dt_f = dt_f * frm_cnf(iw)%sm_el_new
            select case( smp_type ) ! T_f / T_i
             case(0) ! Simple sampling 1d
               T_ratio = dt_i / dt_f
             case(1) ! Gaussian move 1d
               dexp_fun = -0.5_dp * dr(i_dir)**2 * ( 1.0/ dt_f - 1.0/ dt_i )
               if (dexp_fun.lt.-30.0_dp) then
                  dexp_fun = 0.0_dp
               else
                  dexp_fun = exp(dexp_fun)
               endif
               T_ratio = sqrt( dt_i / dt_f ) * dexp_fun
             case(2) ! Simple sampling 3d
               T_ratio = (dt_i / dt_f )**3
             case(3) ! Gaussian move 3d
               dexp_fun = -0.5_dp * sum(dr(1:3)**2) * ( 1.0/ dt_f - 1.0/ dt_i )
               if (dexp_fun.lt.-30.0_dp) then
                  dexp_fun = 0.0_dp
               else
                  dexp_fun = exp(dexp_fun)
               endif
               T_ratio = ( dt_i / dt_f )**1.5 * dexp_fun
            end select
         else
            T_ratio = 1.0_dp
         endif
      else
         T_ratio = 1.0_dp
         dt_f = dt_p / m_p
      endif
      if ( two_step_sampling ) then
         call wvfn%rat1b( iw, g )
      else
         call wvfn%ratio( iw, g )
      endif
      if( qdosys_prs ) then
         call drd_cnf(iw)%new( iw )
         call wvfnq(iw)%ratio( iw, g_fq )
      else
         g_fq = 1.0_dp
      endif
      move_acc = .false.
      if ( two_step_sampling ) then
         g = g**2  * T_ratio
      else
         g = ( g * g_fq )**2 * T_ratio
      endif
      if( g.ge.1.0_dp ) then
         move_acc = .true.
      else
         if( g.ge.rndn(iw)%rndo() ) move_acc = .true.
      endif
      if ( two_step_sampling ) then 
         if (  move_acc ) then
            call wvfn%rat2b( iw, g )
            g = (g*g_fq)**2
            move_acc = .false.
            if(g.ge.1.0_dp) then
               move_acc = .true.
            else
               if( g.ge.rndn(iw)%rndo() ) move_acc = .true.
            endif
         endif
      endif
      if( count ) then
         if( move_acc ) then
            if ( frm_cnf(iw)%i_fe.le.n_el ) then
               n_acc_mv_pw_e(iw) = n_acc_mv_pw_e(iw) + 1_int32
            else
               n_acc_mv_pw_p(iw) = n_acc_mv_pw_p(iw) + 1_int32
            endif
         endif
         if ( frm_cnf(iw)%i_fe.le.n_el ) then
            n_tot_mv_pw_e(iw) = n_tot_mv_pw_e(iw) + 1_int32
         else
            n_tot_mv_pw_p(iw) = n_tot_mv_pw_p(iw) + 1_int32
         endif
      endif
   end subroutine smpl_mtrpls_frm
   subroutine smpl_lngvin_frm( iw, move_acc, count, dr2, pdr2 )
      integer(int32),           intent(in)    :: iw
      logical,                  intent(in)    :: count
      logical,                  intent(out)   :: move_acc
      real(dp), optional,       intent(inout) :: dr2, pdr2
      real(dp), dimension(6)    :: vec_v
      real(dp), dimension(3)    :: d_r, dr_diff, r_fin
      real(dp)                    :: g, g_eq, dt, q_tilde, a, zeta
      integer(int32)              :: i_at_min
      external                    :: nearest_atom
      external                    :: a_vel_rscl
      external                    :: rsc_drift_vel
      external                    :: smpl_drift_diff
      external                    :: UNRr_drift_diff
      real(dp), external          :: gaus_func_3d, exp_func_3d, trns_prob
      call frm_cnf(iw)%chs( iw )
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         dt = dt_e / m_e
      else
         dt = dt_p / m_p
      endif
      vec_v(1:3) = wvfn%DD2_ln_wvfn(4*frm_cnf(iw)%i_fe-3:4*frm_cnf(iw)%i_fe-1,iw)
      if ( jstqe_prs ) then
         vec_v(1:3) = vec_v(1:3) + jstqe_fct(iw)%DD2e_jstqe(4*frm_cnf(iw)%i_fe-3:4*frm_cnf(iw)%i_fe-1)
      endif
      if (nuc_corr.and.frm_cnf(iw)%i_fe.le.n_el) then
         call nearest_atom( frm_cnf(iw)%d_fn(:,frm_cnf(iw)%i_fe), i_at_min )
      else
         i_at_min = 0_int32
      endif
      if (i_at_min.ne.0_int32) then
         call a_vel_rscl( i_at_min, frm_cnf(iw)%d_fn(:,frm_cnf(iw)%i_fe), vec_v(1:3), a )
         call rsc_drift_vel( dt, a, vec_v(1:3) )
         call UNRc_drift( i_at_min, dt, frm_cnf(iw)%d_fn(:,frm_cnf(iw)%i_fe), vec_v(1:3), &
         & zeta, q_tilde, d_r )
         call UNRc_diff( iw, i_at_min, dt, zeta, q_tilde, d_r, dr_diff, r_fin )
      else
         a = a_drft
         q_tilde = 0.0_dp
         call rsc_drift_vel( dt, a, vec_v(1:3) )
         call smpl_drift( dt, vec_v(1:3), frm_cnf(iw)%r_fe(:,frm_cnf(iw)%i_fe), d_r )
         call gaus_dist_3d( iw, sqrt(dt), dr_diff )
         r_fin = d_r + dr_diff
      endif
      frm_cnf(iw)%r_fe_new = r_fin
      call frm_cnf(iw)%new( )
      if (qdosys_prs) call drd_cnf(iw)%new( iw )
      call wvfn%ratio( iw, g )
      if( qdosys_prs ) then
         call wvfnq(iw)%ratio( iw, g_eq )
         g = g * g_eq
      endif
      if ( g.ge.0.0_dp ) then
         g = g**2
         if (i_at_min.eq.0_int32) then
            g = g / gaus_func_3d( dt, dr_diff )
         else
            g = g / trns_prob( dt, zeta, q_tilde, r_fin - d_r, r_fin - r_at(:,i_at_min) )
         endif
         call wvfn%drift( iw )
         vec_v(4:6) = wvfn%DD2_ln_wvfn_new(1:3,iw)
         if ( jstqe_prs ) then
            call wvfnq(iw)%drift( iw )
            vec_v(4:6) = vec_v(4:6) + wvfnq(iw)%DD2_ln_wvfnq_new(1:3)
         endif
         if ( nuc_corr .and.frm_cnf(iw)%i_fe.le.n_el) then
            call nearest_atom( frm_cnf(iw)%d_fn_new(:), i_at_min )
         else
            i_at_min = 0_int32
         endif
         if (i_at_min.ne.0_int32) then
            call a_vel_rscl( i_at_min, frm_cnf(iw)%d_fn_new, vec_v(4:6), a )
            call rsc_drift_vel( dt, a, vec_v(4:6) )
            call UNRc_drift( i_at_min, dt, frm_cnf(iw)%d_fn_new, vec_v(4:6), zeta, q_tilde, d_r )
            g = g * trns_prob( dt, zeta, q_tilde, frm_cnf(iw)%r_fe(:,frm_cnf(iw)%i_fe) - d_r, &
            & frm_cnf(iw)%r_fe(:,frm_cnf(iw)%i_fe) -r_at(:,i_at_min) )
         else
            a = a_drft
            q_tilde = 0.0_dp
            call rsc_drift_vel( dt, a, vec_v(4:6) )
            call smpl_drift( dt, vec_v(4:6), frm_cnf(iw)%r_fe_new, d_r )
            g = g * gaus_func_3d( dt, frm_cnf(iw)%r_fe(:,frm_cnf(iw)%i_fe) - d_r )
         endif
         move_acc = .false.
         if(g.ge.1.0_dp) then
            move_acc = .true.
            g = 1.0_dp
         else
            if( g.ge.rndn(iw)%rndo() ) move_acc = .true.
         endif
      else ! g < 0
         move_acc = .false.
         g = 0.0_dp
      endif
      if (present(dr2)) then
         if ( frm_cnf(iw)%i_fe.le.n_el ) then
            a = sum(dr_diff(1:3)**2) * m_e
            pdr2 = pdr2 + dt_e * g * a
         else
            a = sum(dr_diff(1:3)**2) * m_p
            pdr2 = pdr2 + dt_p * g * a
         endif
         dr2  = dr2 + a
      endif
      if( count ) then
         if( move_acc ) then
            if ( frm_cnf(iw)%i_fe.le.n_el ) then
               n_acc_mv_pw_e(iw) = n_acc_mv_pw_e(iw) + 1_int32
            else
               n_acc_mv_pw_p(iw) = n_acc_mv_pw_p(iw) + 1_int32
            endif
         endif
         if ( frm_cnf(iw)%i_fe.le.n_el ) then
            n_tot_mv_pw_e(iw) = n_tot_mv_pw_e(iw) + 1_int32
         else
            n_tot_mv_pw_p(iw) = n_tot_mv_pw_p(iw) + 1_int32
         endif
      endif
   end subroutine smpl_lngvin_frm
end module fermionic_sampling_m
