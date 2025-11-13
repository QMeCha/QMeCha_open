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
subroutine cmp_dmc_wght( iw, e_l_new, e_l_old, e_l_pp_new, e_l_pp_old, w )
   use fortran_kinds_v,         only: dp, int32
   use physical_constants_v,    only: pi
   use quantum_monte_carlo_v,   only: dt
   use diffusion_monte_carlo_v, only: encrr_type, e_r, e_best, e_best_var, g_e,&
   & V2, log_w_max, w_max, V2_pp, spe_best, spe_best_var, trotter_order, dt_eff,&
   & fe_i_at_old
   use molecular_system_v,      only: n_fe, n_el, n_po, n_at
   use qdo_system_v,            only: n_qdo
   use drudonic_config_c,       only: drd_cnf
   use fermionic_config_c,      only: frm_cnf
   implicit none
   integer(int32),                    intent(in)    :: iw
   real(dp),                          intent(in)    :: e_l_new
   real(dp),                          intent(in)    :: e_l_old
   real(dp), dimension(n_fe+n_qdo), intent(in)    :: e_l_pp_new
   real(dp), dimension(n_fe+n_qdo), intent(in)    :: e_l_pp_old
   real(dp),                          intent(inout) :: w
   real(dp), external :: erf_1
   real(dp)       :: e_tmp, dt_rescaled, e_cut, s_r, e_abs, e_sgn
   integer(int32) :: i1, i2, i_min
   if (encrr_type.eq.0) then
      if (trotter_order.eq.1) then
         s_r = e_r - e_l_new
      else if (trotter_order.eq.2) then
         s_r = 2.0_dp *e_r - (e_l_old + e_l_new)
      endif
   else 
      dt_rescaled = 1.0_dp / sqrt(dt_eff(iw))
      if (trotter_order.eq.1) then
         s_r = (e_r - e_best)
      else if (trotter_order.eq.2) then
         s_r = 2.0_dp * (e_r - e_best)
      endif
      select case(encrr_type)
       case(1) ! Zen et. al, Ye Luo et al., De Pasquale et al.
         e_cut = g_e * sqrt( dble(n_fe+n_qdo)/dt_eff(iw) )
         if (trotter_order.eq.2) then
            e_tmp = e_best - e_l_old
            if ( e_tmp.gt.e_cut) e_tmp = e_cut
            s_r = s_r + e_tmp
         endif
         e_tmp = e_best - e_l_new
         if ( e_tmp.gt.e_cut) e_tmp = e_cut
         s_r = s_r + e_tmp
       case(3) ! Umrigar et al. 1993
         if (trotter_order.eq.2) then
            e_tmp = e_best - e_l_old
            s_r = s_r + e_tmp * sqrt(V2(iw))
         endif
         call cmp_V2_V2_bar_ratio(iw, V2(iw))
         e_tmp = e_best - e_l_new
         s_r = s_r + e_tmp * sqrt(V2(iw))
       case(4) ! Anderson Umrigar 2021
         if (trotter_order.eq.2) then
            e_tmp = e_best - e_l_old
            s_r = s_r + e_tmp /(1.0 + (V2(iw) * dt_eff(iw)**2 ))
         endif
         call cmp_V2_avg(iw, V2(iw))
         e_tmp = e_best - e_l_new
         s_r = s_r + e_tmp /(1.0 + (V2(iw) * dt_eff(iw)**2 ))
       case(5) ! Anderson Umrigar 2021
         if (trotter_order.eq.2) then
            e_tmp = e_best - e_l_old
            s_r = s_r + e_tmp * erf_1(1.0/(sqrt(V2(iw)) * dt_eff(iw) ) )
         endif
         call cmp_V2_avg(iw, V2(iw))
         e_tmp = e_best - e_l_new
         s_r = s_r + e_tmp * erf_1(1.0/(sqrt(V2(iw)) * dt_eff(iw)**2 ) )
       case(11)
         e_cut = g_e / sqrt(dt_eff(iw))
         do i1 = 1_int32, n_fe + n_qdo
            if (trotter_order.eq.2) then
               e_tmp = spe_best(i1) - e_l_pp_old(i1)
               if ( e_tmp.gt.e_cut) e_tmp = e_cut
               s_r = s_r + e_tmp
            endif
            e_tmp = spe_best(i1) - e_l_pp_new(i1)
            if ( e_tmp.gt.e_cut) e_tmp = e_cut
            s_r = s_r + e_tmp
         enddo
       case(13) ! Umrigar et al. Particle by Particle 1993
         if (trotter_order.eq.2) then
            do i1 = 1_int32, n_fe + n_qdo
               e_tmp = spe_best(i1) - e_l_pp_old(i1)
               s_r = s_r + e_tmp * V2_pp(i1,iw)
            enddo
         endif
         call cmp_f_UNR93_vec(iw, V2_pp(:,iw))
         do i1 = 1_int32, n_fe + n_qdo
            e_tmp = spe_best(i1) - e_l_pp_new(i1)
            s_r = s_r + e_tmp * V2_pp(i1,iw)
         enddo
       case(14) ! Umrigar et al. Particle by Particle 2021
         if (trotter_order.eq.2) then
            do i1 = 1_int32, n_fe + n_qdo
               e_tmp = spe_best(i1) - e_l_pp_old(i1)
               s_r = s_r + e_tmp / (1.0 + V2_pp(i1,iw)*dt_eff(iw)**2)
            enddo
         endif
         call cmp_V2_vec(iw, V2_pp(:,iw))
         do i1 = 1_int32, n_fe + n_qdo
            e_tmp = spe_best(i1) - e_l_pp_new(i1)
            s_r = s_r + e_tmp / (1.0 + V2_pp(i1,iw)*dt_eff(iw)**2)
         enddo
       case(15) ! Umrigar et al. Particle by Particle 2021
         if (trotter_order.eq.2) then
            do i1 = 1_int32, n_fe + n_qdo
               e_tmp = spe_best(i1) - e_l_pp_old(i1)
               s_r = s_r + e_tmp  *V2_pp(i1,iw)
            enddo
         endif
         call cmp_F_vec(iw, dt_eff(iw), V2_pp(:,iw))
         do i1 = 1_int32, n_fe + n_qdo
            e_tmp = spe_best(i1) - e_l_pp_new(i1)
            s_r = s_r + e_tmp  *V2_pp(i1,iw)
         enddo
       case(21)
         e_cut = g_e / sqrt(dt_eff(iw))
         i2 = 0_int32
         if (n_el.gt.0_int32) then
            do i1 = 1_int32, n_el
               if (trotter_order.eq.2) then
                  i_min = fe_i_at_old(i1,iw)
                  e_tmp = spe_best(i_min) - e_l_pp_old(i1)
                  if ( e_tmp.gt.e_cut) e_tmp = e_cut
                  s_r = s_r + e_tmp
               endif
               i_min = frm_cnf(iw)%fe_i_at(i1)
               e_tmp = spe_best(i_min) - e_l_pp_new(i1)
               if ( e_tmp.gt.e_cut) e_tmp = e_cut
               s_r = s_r + e_tmp
            enddo
            i2 = i2 + n_at
         endif
         if (n_po.gt.0_int32) then
            do i1 = n_el+1_int32, n_fe
               if (trotter_order.eq.2) then
                  i_min = i2 + fe_i_at_old(i1,iw)
                  e_tmp = spe_best(i_min) - e_l_pp_old(i1)
                  if ( e_tmp.gt.e_cut) e_tmp = e_cut
                  s_r = s_r + e_tmp
               endif
               i_min = i2 + frm_cnf(iw)%fe_i_at(i1)
               e_tmp = spe_best(i_min) - e_l_pp_new(i1)
               if ( e_tmp.gt.e_cut) e_tmp = e_cut
               s_r = s_r + e_tmp
            enddo
            i2 = i2 + n_at
         endif
       case(24) ! Umrigar et al. Particle by Particle 2021
         if (trotter_order.eq.2) then
            do i1 = 1_int32, n_el
               i_min = fe_i_at_old(i1,iw)
               e_tmp = spe_best(i_min) - e_l_pp_old(i1)
               s_r = s_r + e_tmp  / (1.0 + V2_pp(i1,iw)*dt_eff(iw)**2)
               s_r = s_r + e_tmp  *V2_pp(i1,iw)
            enddo
         endif
         call cmp_V2_vec(iw, V2_pp(:,iw))
         do i1 = 1_int32, n_el
            i_min = frm_cnf(iw)%fe_i_at(i1)
            e_tmp = spe_best(i_min) - e_l_pp_new(i1)
            s_r = s_r + e_tmp  / (1.0 + V2_pp(i1,iw)*dt_eff(iw)**2)
         enddo
       case(25) ! Umrigar et al. Particle by Particle 2021
         if (trotter_order.eq.2) then
            do i1 = 1_int32, n_el
               i_min = fe_i_at_old(i1,iw)
               e_tmp = spe_best(i_min) - e_l_pp_old(i1)
               s_r = s_r + e_tmp *V2_pp(i1,iw)
            enddo
         endif
         call cmp_F_vec(iw, dt, V2_pp(:,iw))
         do i1 = 1_int32, n_el
            i_min = frm_cnf(iw)%fe_i_at(i1)
            e_tmp = spe_best(i_min) - e_l_pp_new(i1)
            s_r = s_r + e_tmp *V2_pp(i1,iw)
         enddo
       case default
      end select
   endif
   if (trotter_order.eq.1) then
      s_r = s_r *dt_eff(iw) + log(w)
   else if (trotter_order.eq.2) then
      s_r = 0.5 * s_r *dt_eff(iw) + log(w)
   endif
   if ( s_r.gt. log_w_max ) then
      w = w_max
   else
      if (s_r.gt.-34.0) then
         w = dexp( s_r )
      else
         w = 0.0_dp
      endif
   endif
end subroutine cmp_dmc_wght
subroutine cmp_V2_avg( iw, V2 )
   use fortran_kinds_v,          only: dp, int32
   use molecular_system_v,       only: n_fe, molsys_prs
   use qdo_system_v,             only: n_qdo, qdosys_prs
   use jstqepar_var,             only: jstqe_prs
   use fermionic_wavefunction_m, only: wvfn
   use qdo_wavefunction_m,               only: wvfnq
   use jstqe_mod,                only: jstqe_fct
   implicit none
   integer(int32), intent(in)  :: iw
   real(dp),       intent(out) :: V2
   real(dp)       :: vec_v(1:3)
   integer(int32) :: i1
   V2 = 0.0_dp
   if (molsys_prs) then
      do i1 = 1_int32, n_fe
         vec_v(1:3) = wvfn%DD2_ln_wvfn(4*i1-3:4*i1-1,iw)
         if ( jstqe_prs ) then
            vec_v(1:3) = vec_v(1:3) + jstqe_fct(iw)%DD2e_jstqe(4*i1-3:4*i1-1)
         endif
         V2 = V2 + sum(vec_v(1:3)**2)
      enddo
   endif
   if (qdosys_prs) then
      do i1 = 1_int32, n_qdo
         V2 = V2 + sum( wvfnq(iw)%DD2_ln_wvfnq( (4*i1-3):(4*i1-1) )**2)
      enddo
   endif
end subroutine cmp_V2_avg
subroutine cmp_V2_vec( iw, V2_pp )
   use fortran_kinds_v,          only: dp, int32
   use molecular_system_v,       only: n_fe, molsys_prs
   use qdo_system_v,             only: n_qdo, qdosys_prs
   use jstqepar_var,             only: jstqe_prs
   use fermionic_wavefunction_m, only: wvfn
   use qdo_wavefunction_m,               only: wvfnq
   use jstqe_mod,                only: jstqe_fct
   implicit none
   integer(int32),                    intent(in)  :: iw
   real(dp), dimension(n_fe+n_qdo), intent(out) :: V2_pp
   real(dp)       :: vec_v(1:3)
   integer(int32) :: i1
   V2_pp = 0.0_dp
   if (molsys_prs) then
      do i1 = 1_int32, n_fe
         vec_v(1:3) = wvfn%DD2_ln_wvfn(4*i1-3:4*i1-1,iw)
         if ( jstqe_prs ) then
            vec_v(1:3) = vec_v(1:3) + jstqe_fct(iw)%DD2e_jstqe(4*i1-3:4*i1-1)
         endif
         V2_pp(i1) = V2_pp(i1) + sum(vec_v(1:3)**2)
      enddo
   endif
   if (qdosys_prs) then
      do i1 = 1_int32, n_qdo
         V2_pp(n_fe+i1) = V2_pp(n_fe+i1) + sum( wvfnq(iw)%DD2_ln_wvfnq( (4*i1-3):(4*i1-1) )**2)
      enddo
   endif
end subroutine cmp_V2_vec
subroutine cmp_F_vec( iw, eps, F )
   use fortran_kinds_v,          only: dp, int32
   use molecular_system_v,       only: n_fe, molsys_prs
   use qdo_system_v,             only: n_qdo, qdosys_prs
   use jstqepar_var,             only: jstqe_prs
   use fermionic_wavefunction_m, only: wvfn
   use qdo_wavefunction_m,       only: wvfnq
   use jstqe_mod,                only: jstqe_fct
   implicit none
   integer(int32),                    intent(in)  :: iw
   real(dp),                          intent(in)  :: eps
   real(dp), dimension(n_fe+n_qdo), intent(out) :: F
   real(dp), external :: regularization_f, erf_1
   real(dp)       :: vec_v(1:3)
   integer(int32) :: i1
   if (molsys_prs) then
      do i1 = 1_int32, n_fe
         vec_v(1:3) = wvfn%DD2_ln_wvfn(4*i1-3:4*i1-1,iw)
         if ( jstqe_prs ) then
            vec_v(1:3) = vec_v(1:3) + jstqe_fct(iw)%DD2e_jstqe(4*i1-3:4*i1-1)
         endif
         F(i1) = 1.0_dp / sum(vec_v(1:3)**2)
      enddo
   endif
   if (qdosys_prs) then
      do i1 = 1_int32, n_qdo
         F(n_fe+i1) = 1.0_dp / sum( wvfnq(iw)%DD2_ln_wvfnq( (4*i1-3):(4*i1-1) )**2)
      enddo
   endif
   F(:) = sqrt(F(:)) / eps
   do i1 = 1_int32, n_fe + n_qdo
      F(i1) = erf_1(  F(i1)  )
   enddo
end subroutine cmp_F_vec
subroutine cmp_f_UNR93_vec( iw, V2_pp )
   use fortran_kinds_v,          only: dp, int32
   use molecular_system_v,       only: n_el, n_fe, molsys_prs, n_at
   use fermionic_config_c,       only: frm_cnf
   use qdo_system_v,             only: n_qdo, qdosys_prs
   use jstqepar_var,             only: jstqe_prs
   use fermionic_wavefunction_m, only: wvfn
   use qdo_wavefunction_m,               only: wvfnq
   use jstqe_mod,                only: jstqe_fct
   use quantum_monte_carlo_v,    only: dt, nuc_corr, a_drft
   implicit none
   integer(int32),                    intent(in)  :: iw
   real(dp), dimension(n_fe+n_qdo), intent(out) :: V2_pp
   real(dp)        :: vec_v(1:3), a
   integer(int32)  :: i1, i_at_min
   V2_pp = 0.0_dp
   if (molsys_prs) then
      do i1 = 1_int32, n_fe
         vec_v(1:3) = wvfn%DD2_ln_wvfn(4*i1-3:4*i1-1,iw)
         if ( jstqe_prs ) then
            vec_v(1:3) = vec_v(1:3) + jstqe_fct(iw)%DD2e_jstqe(4*i1-3:4*i1-1)
         endif
         a = a_drft
         if (i1.le.n_el.and.nuc_corr) then
            i_at_min = 0_int32
            call nearest_atom( frm_cnf(iw)%d_fn(1:n_at,i1), i_at_min )
            call a_vel_rscl( i_at_min, frm_cnf(iw)%d_fn(1:n_at,i1), vec_v(1:3), a )
         endif
         a = a*sum(vec_v(1:3)**2)*dt
         V2_pp(i1) = (sqrt(1.0_dp+2.0_dp*a) -1.0_dp ) / a
      enddo
   endif
   if (qdosys_prs) then
      do i1 = 1_int32, n_qdo
         vec_v(1:3) =  wvfnq(iw)%DD2_ln_wvfnq( (4*i1-3):(4*i1-1) )
         a = a_drft
         a = a*sum(vec_v(1:3)**2)*dt
         V2_pp(n_fe+i1) = (sqrt(1.0_dp+2.0_dp*a) -1.0_dp ) / a
      enddo
   endif
end subroutine cmp_f_UNR93_vec
subroutine cmp_V2_V2_bar_ratio( iw, V2_bar )
   use fortran_kinds_v, only: dp, int32
   use quantum_monte_carlo_v, only: dt, nuc_corr, a_drft
   use molecular_system_v, only: n_fe, molsys_prs, n_at
   use fermionic_config_c, only: frm_cnf
   use qdo_system_v, only: n_qdo, qdosys_prs
   use jstqepar_var, only: jstqe_prs
   use fermionic_wavefunction_m, only: wvfn
   use qdo_wavefunction_m, only: wvfnq
   use jstqe_mod,  only: jstqe_fct
   implicit none
   integer(int32), intent(in)  :: iw
   real(dp),       intent(out) :: V2_bar
   real(dp)       :: vec_v(1:3), a, V2
   integer(int32) :: i1, i_at_min
   V2 = 0.0_dp ; V2_bar = 0.0_dp
   if (molsys_prs) then
      do i1 = 1_int32, n_fe
         vec_v(1:3) = wvfn%DD2_ln_wvfn(4*i1-3:4*i1-1,iw)
         if ( jstqe_prs ) then
            vec_v(1:3) = vec_v(1:3) + jstqe_fct(iw)%DD2e_jstqe(4*i1-3:4*i1-1)
         endif
         V2 = V2 + sum(vec_v(1:3)**2)
         a = a_drft
         if (nuc_corr) then
            i_at_min = 0_int32
            call nearest_atom( frm_cnf(iw)%d_fn(1:n_at,i1), i_at_min )
            call a_vel_rscl( i_at_min, frm_cnf(iw)%d_fn(1:n_at,i1), vec_v(1:3), a )
         endif
         call rsc_drift_vel( dt, a, vec_v(1:3) )
         V2_bar = V2_bar + sum(vec_v(1:3)**2)
      enddo
   endif
   if (qdosys_prs) then
      do i1 = 1_int32, n_qdo
         vec_v(1:3) = wvfnq(iw)%DD2_ln_wvfnq( (4*i1-3):(4*i1-1) )
         V2 = V2 + sum(vec_v(1:3)**2)
         a= a_drft
         call rsc_drift_vel( dt, a, vec_v(1:3) )
         V2_bar = V2_bar + sum(vec_v(1:3)**2)
      enddo
   endif
   V2_bar = V2_bar / V2
end subroutine cmp_V2_V2_bar_ratio
subroutine init_dmc_rene( )
   use fortran_kinds_v,         only: dp, int32
   use openmp_mpi_m
   use physical_constants_v,    only: pi
   use write_lines_m,           only: write_variable_line
   use quantum_monte_carlo_v,   only: dt, n_wlk, n_tot_wlk
   use diffusion_monte_carlo_v, only: encrr_type, e_r, e_r_avg, e_best, e_best_var, &
   & g_e, V2, V2_pp, spe_best, spe_best_var, w, dt_eff_avg
   use molecular_system_v,      only: n_fe, n_el, n_po, n_at
   use local_energy_m,          only: e_l, e_l_pp
   use qdo_system_v,            only: n_qdo
   use drudonic_config_c,       only: drd_cnf
   use fermionic_config_c,      only: frm_cnf
   implicit none
   real(dp)       :: e_tmp, e_cut, s_r, e_abs, e_sgn, rescaled_weight
   integer(int32) :: i1, i2, i_min
   integer(int32) :: iw
   rescaled_weight = 0.0_dp
!$omp parallel default(shared) private(iw,e_tmp,e_cut,e_abs, e_sgn, s_r) reduction(+:rescaled_weight)
!$omp do
   do iw = 1_int32, n_wlk
      if (encrr_type.eq.0) then 
         s_r = (e_best - e_l(iw))
      else
         select case(encrr_type)
          case(1) ! Zen et. al, Ye Luo et al., De Pasquale et al.
            e_cut = g_e * sqrt( dble(n_fe+n_qdo)/dt )
            e_tmp = e_best - e_l(iw)
            if ( e_tmp.gt.e_cut) e_tmp = e_cut
            s_r = e_tmp
          case(3) ! Umrigar et al. 1993
            call cmp_V2_V2_bar_ratio(iw, V2(iw))
            e_tmp = e_best - e_l(iw)
            s_r = e_tmp * sqrt(V2(iw))
          case(4) ! Anderson Umrigar 2021
            call cmp_V2_avg(iw, V2(iw))
            e_tmp = e_best - e_l(iw)
            s_r = e_tmp /(1.0 + V2(iw) * dt**2)
          case(5) ! Anderson Umrigar 2021
            call cmp_V2_avg(iw, V2(iw))
            e_tmp = e_best - e_l(iw)
            s_r = e_tmp * erf(1.0/(sqrt(V2(iw)) * dt ) )
          case(11)
            s_r = 0.0_dp
            e_cut = g_e / sqrt(dt)
            do i1 = 1_int32, n_fe + n_qdo
               e_tmp = spe_best(i1) - e_l_pp(i1,iw)
               if ( e_tmp.gt.e_cut) e_tmp = e_cut
               s_r = s_r + e_tmp
            enddo
          case(13) ! Umrigar et al. Particle by Particle 1993
            s_r = 0.0_dp
            call cmp_f_UNR93_vec(iw, V2_pp(:,iw))
            do i1 = 1_int32, n_fe + n_qdo
               e_tmp = spe_best(i1) - e_l_pp(i1,iw)
               s_r = s_r + e_tmp * V2_pp(i1,iw)
            enddo
          case(14) ! Umrigar et al. Particle by Particle 2021
            s_r = 0.0_dp
            call cmp_V2_vec(iw, V2_pp(:,iw))
            do i1 = 1_int32, n_fe + n_qdo
               e_tmp = spe_best(i1) - e_l_pp(i1,iw)
               s_r = s_r + e_tmp / (1.0 + V2_pp(i1,iw)*dt**2)
            enddo
          case(15) ! Umrigar et al. Particle by Particle 2021
            s_r = 0.0_dp
            call cmp_F_vec(iw, dt ,V2_pp(:,iw))
            do i1 = 1_int32, n_fe + n_qdo
               e_tmp = spe_best(i1) - e_l_pp(i1,iw)
               s_r = s_r + e_tmp * V2_pp(i1,iw)
            enddo
          case(21)
            s_r = 0.0_dp
            e_cut = g_e / sqrt(dt)
            i2 = 0_int32
            if (n_el.gt.0_int32) then
               do i1 = 1_int32, n_el
                  i_min = frm_cnf(iw)%fe_i_at(i1)
                  e_tmp = spe_best(i_min) - e_l_pp(i1,iw)
                  if ( e_tmp.gt.e_cut) e_tmp = e_cut
                  s_r = s_r + e_tmp
               enddo
               i2 = i2 + n_at
            endif
            if (n_po.gt.0_int32) then
               do i1 = n_el+1_int32, n_fe
                  i_min = i2 + frm_cnf(iw)%fe_i_at(i1)
                  e_tmp = spe_best(i_min) - e_l_pp(i1,iw)
                  if ( e_tmp.gt.e_cut) e_tmp = e_cut
                  s_r = s_r + e_tmp
               enddo
               i2 = i2 + n_at
            endif
          case(24) ! Umrigar et al. Particle by Particle 2021
            s_r = 0.0_dp
            call cmp_V2_vec(iw, V2_pp(:,iw))
            do i1 = 1_int32, n_fe
               i_min = frm_cnf(iw)%fe_i_at(i1)
               e_tmp = spe_best(i_min) - e_l_pp(i1,iw)
               s_r = s_r + e_tmp / (1.0 + V2_pp(i1,iw)*dt**2)
            enddo
          case(25) ! Umrigar et al. Particle by Particle 2021
            s_r = 0.0_dp
            call cmp_F_vec(iw, dt, V2_pp(:,iw))
            do i1 = 1_int32, n_fe + n_qdo
               i_min = frm_cnf(iw)%fe_i_at(i1)
               e_tmp = spe_best(i_min) - e_l_pp(i1,iw)
               s_r = s_r + e_tmp * V2_pp(i1,iw)
            enddo
          case default
         end select
      endif
      rescaled_weight = rescaled_weight + exp(dt*s_r)
   enddo
!$omp end do
!$omp end parallel
#if defined _MPI || defined _MPIh || defined _MPI08
   if (n_mpi_tasks.gt.1_int32) then
      call mpi_allreduce(MPI_IN_PLACE, rescaled_weight, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
   endif
#endif
   e_r = e_best - 1.0_dp /dt * dlog ( rescaled_weight / dble(n_tot_wlk) )
   e_r_avg = e_r
   call write_empty_line(stdout,0,mpi_rank)
   call write_variable_line(stdout,0,mpi_rank,2,"Initialized reference Energy",e_r_avg,units='Eh')
end subroutine init_dmc_rene
subroutine cmp_dmc_rene( e_ini, G_w_tot, G_w_trgt )
   use fortran_kinds_v,         only: dp
   use diffusion_monte_carlo_v, only: g_r, e_r, dt_eff_avg
   use quantum_monte_carlo_v,   only: dt
   implicit none
   real(dp), intent(in) :: e_ini
   real(dp), intent(in) :: G_w_tot
   real(dp), intent(in) :: G_w_trgt
   e_r = e_ini - g_r * log ( G_w_tot / G_w_trgt )
end subroutine cmp_dmc_rene
