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
subroutine a_vel_rscl( i_at_min, d_fn, vec_v, a )
   use fortran_kinds_v, only: dp, int32
   use molecular_system_v, only: n_at, atoms, dist_t
   implicit none
   integer(int32),                  intent(in)    :: i_at_min
   type(dist_t), dimension(n_at), intent(in)    :: d_fn
   real(dp), dimension(3), intent(in) :: vec_v
   real(dp)                , intent(out) :: a
   real(dp)                 :: mod_v, Zz2
   real(dp), dimension(3) :: hat_v
   real(dp), dimension(3) :: hat_d_fn
   mod_v = sqrt(sum(vec_v(1:3)**2))
   if ( mod_v.ne.0.0_dp ) then
      hat_v = vec_v / mod_v
   else
      hat_v = 0.0_dp
   endif
   if ( d_fn(i_at_min)%m.ne.0.0_dp ) then
      hat_d_fn = d_fn(i_at_min)%v / d_fn(i_at_min)%m
   else
      hat_d_fn = 0.0_dp
   endif
   Zz2 = (dble(atoms(i_at_min)%atm_z) * d_fn(i_at_min)%m)**2
   a = 0.1_dp * Zz2 / ( 4.0_dp + Zz2 ) + 0.5_dp * (1.0_dp + dot_product( hat_v, hat_d_fn ) )
end subroutine a_vel_rscl
subroutine rsc_drift_vel( dt, a, vec_v )
   use fortran_kinds_v, only: dp
   implicit none
   real(dp),                 intent(in)    :: dt
   real(dp),                 intent(inout) :: a
   real(dp), dimension(3), intent(inout) :: vec_v
   a = a*sum(vec_v(1:3)**2)*dt
   if (a.gt.0.0_dp) then
      vec_v = vec_v * (sqrt(1.0_dp+2.0_dp*a) -1.0_dp ) / a
   endif
end subroutine rsc_drift_vel
subroutine smpl_drift( dt, vec_v, r_ini, r_fin )
   use fortran_kinds_v, only: dp
   implicit none
   real(dp),                       intent(in)    :: dt
   real(dp),       dimension(3), intent(in)    :: vec_v
   real(dp),       dimension(3), intent(in)    :: r_ini
   real(dp),       dimension(3), intent(inout) :: r_fin
   r_fin(:) = r_ini(:) + dt * vec_v(:)
end subroutine smpl_drift
subroutine UNRc_drift( i_at_min, dt, d_fn, vec_v, zeta, q_tilde, d_r )
   use fortran_kinds_v, only: dp, int32
   use molecular_system_v, only: n_at, atoms, dist_t, r_at
   implicit none
   integer(int32),                     intent(in)    :: i_at_min
   real(dp),                           intent(in)    :: dt
   type(dist_t),       dimension(n_at), intent(in) :: d_fn
   real(dp),           dimension(3), intent(inout) :: vec_v
   real(dp),                           intent(inout) :: zeta
   real(dp),                           intent(inout) :: q_tilde
   real(dp),           dimension(3), intent(inout) :: d_r
   real(dp), dimension(3)                :: hat_d_fn
   real(dp) :: mod_v_p
   real(dp) :: mod_v_o
   real(dp), dimension(3) :: hat_v_o
   real(dp) :: mod_dr_p, mod_dr_o
   if ( d_fn(i_at_min)%m.ne.0.0_dp ) then
      hat_d_fn = d_fn(i_at_min)%v / d_fn(i_at_min)%m
   else
      hat_d_fn = 0.0_dp
   endif
   zeta = dsqrt( atoms(i_at_min)%atm_z**2 + 1.0_dp / dt )
   mod_v_p = dot_product( vec_v, hat_d_fn )
   hat_v_o = vec_v - mod_v_p * hat_d_fn
   mod_v_o = sqrt(sum(hat_v_o(1:3)**2))
   if ( mod_v_o.ne.0.0_dp ) then
      hat_v_o = hat_v_o / mod_v_o
   else
      hat_v_o = 0.0_dp
   endif
   mod_dr_p = d_fn(i_at_min)%m + mod_v_p * dt
   if ( mod_dr_p / sqrt(2.0_dp*dt) .gt. 25.0_dp ) then
      q_tilde = 0.0_dp
   else
      q_tilde = 0.5_dp*erfc(mod_dr_p / sqrt(2.0_dp*dt))
   endif
   mod_dr_p = max(mod_dr_p,0.0_dp)
   if( mod_dr_p.ne.0.0_dp ) then
      mod_dr_o = 2.0_dp*mod_dr_p/(d_fn(i_at_min)%m + mod_dr_p) * mod_v_o * dt
   else
      mod_dr_o = 0.0_dp
   endif
   d_r = r_at(:,i_at_min) + mod_dr_o * hat_v_o + mod_dr_p * hat_d_fn
end subroutine UNRc_drift
subroutine UNRc_diff( iw, i_at_min, dt, zeta, q_tilde, d_r, xi, r_fin )
   use fortran_kinds_v, only: dp, int32
   use molecular_system_v, only: dist_t, r_at
   use mersenne_twister19937_m, only: rndn
   implicit none
   integer(int32),                     intent(in)    :: iw
   integer(int32),                     intent(in)    :: i_at_min
   real(dp),                           intent(in)    :: dt
   real(dp),                           intent(in)    :: zeta
   real(dp),                           intent(in)    :: q_tilde
   real(dp),           dimension(3), intent(in)    :: d_r
   real(dp),           dimension(3), intent(inout) :: r_fin
   real(dp),           dimension(3), intent(out) :: xi
   if ( (1.0_dp - rndn(iw)%rndo()).lt.q_tilde ) then
      call exp_dist_3d( iw, zeta, xi )
      r_fin = r_at(:,i_at_min) + xi
   else
      call gaus_dist_3d( iw, sqrt(dt), xi )
      r_fin = d_r + xi
   endif
end subroutine UNRc_diff
function trns_prob( dt, zeta, q_tilde, dr_gauss, dr_exp ) result( t_prob )
   use fortran_kinds_v, only: dp
   implicit none
   real(dp),                           intent(in)    :: dt
   real(dp),                           intent(in)    :: zeta
   real(dp),                           intent(in)    :: q_tilde
   real(dp), dimension(3), intent(in)    :: dr_gauss
   real(dp), dimension(3), intent(in)    :: dr_exp
   real(dp)                                :: t_prob
   real(dp), external :: gaus_func_3d, exp_func_3d
   if (q_tilde.lt.10d-16) then
      t_prob = gaus_func_3d(dt,dr_gauss)
   else
      t_prob = (1.0_dp - q_tilde) * gaus_func_3d(dt,dr_gauss) &
      & + q_tilde * exp_func_3d(zeta, dr_exp)
   endif
end function trns_prob
subroutine exp_dist_3d( iw, zeta, r )
   use fortran_kinds_v, only: dp, int32
   use mersenne_twister19937_m, only: rndn
   use physical_constants_v, only: twopi
   implicit none
   integer(int32),           intent(in)    :: iw
   real(dp),                 intent(in)    :: zeta
   real(dp), dimension(3), intent(inout) :: r
   real(dp) :: rand1,rand2,rand3
   real(dp) :: r_mod, cos_theta, sin_theta,phi
   rand1 = 1.0_dp - rndn(iw)%rndo()
   rand2 = 1.0_dp - rndn(iw)%rndo()
   rand3 = 1.0_dp - rndn(iw)%rndo()
   r_mod = -log (rand1*rand2*rand3) / (2.0_dp * zeta)
   cos_theta = 1.0_dp - 2.0_dp * rndn(iw)%rndo()
   sin_theta = sqrt(1.0_dp-cos_theta**2)
   phi       = rndn(iw)%rndo() * twopi
   r(1) = r_mod * sin_theta * cos(phi)
   r(2) = r_mod * sin_theta * sin(phi)
   r(3) = r_mod * cos_theta
end subroutine exp_dist_3d
subroutine gaus_dist_3d( iw, sqrt_dt, r )
   use fortran_kinds_v, only: dp, int32
   use mersenne_twister19937_m, only: rndn
   implicit none
   integer(int32),           intent(in)    :: iw
   real(dp),                 intent(in)    :: sqrt_dt
   real(dp), dimension(3), intent(inout) :: r
   r(1) = rndn(iw)%rndbm()
   r(2) = rndn(iw)%rndbm()
   r(3) = rndn(iw)%rndbm()
   r(:) = sqrt_dt * r(:)
end subroutine gaus_dist_3d
function gaus_func_3d( dt, r ) result(g)
   use fortran_kinds_v, only: dp
   use physical_constants_v, only: twopi
   implicit none
   real(dp),                 intent(in) :: dt
   real(dp), dimension(3), intent(in) :: r
   real(dp)                             :: g
   g = dexp (-0.5_dp * sum(r(1:3)**2) / dt) / dsqrt(twopi * dt) ** 3
end function gaus_func_3d
function exp_func_3d( zeta, r ) result(g)
   use fortran_kinds_v, only: dp
   use physical_constants_v, only: pi
   implicit none
   real(dp),                 intent(in) :: zeta
   real(dp), dimension(3), intent(in) :: r
   real(dp)                             :: g
   g = dexp ( -2.0_dp*zeta*sqrt(r(1)**2+r(2)**2+r(3)**2) )  * zeta ** 3 / pi
end function exp_func_3d
subroutine nearest_atom( d_fn, i_at_min )
   use fortran_kinds_v, only: dp, int32
   use molecular_system_v, only: n_at, atoms, dist_t
   implicit none
   type(dist_t), dimension(n_at), intent(in)  :: d_fn
   integer(int32),                intent(out) :: i_at_min
   real(dp)                                   :: d_fn_min
   integer(int32) :: i1
   i_at_min = 0_int32 ; d_fn_min = 0.0_dp
   if ( atoms(1)%atm_z.gt.0_int32 ) then
      d_fn_min = d_fn(1)%m
      i_at_min = 1_int32
   endif
   if (n_at.gt.1_int32) then
      do i1 = 2_int32, n_at
         if( atoms(i1)%atm_z.gt.0_int32 .and. d_fn(i1)%m.le.d_fn_min ) then
            i_at_min = i1
            d_fn_min = d_fn(i1)%m
         endif
      enddo
   endif
end subroutine nearest_atom
subroutine nearest_qdo( d_dq, i_qdo_min )
   use fortran_kinds_v, only: dp, int32
   use qdo_system_v, only: n_qdo, qdos, dist_t
   implicit none
   type(dist_t), dimension(n_qdo), intent(in)  :: d_dq
   integer(int32),                 intent(out) :: i_qdo_min
   real(dp)                                       :: d_dq_min
   integer(int32) :: i1
   d_dq_min = d_dq(1)%m
   i_qdo_min = 1_int32
   if (n_qdo.gt.1_int32) then
      do i1 = 2_int32, n_qdo
         if( d_dq(i1)%m.le.d_dq_min ) then
            i_qdo_min = i1
            d_dq_min = d_dq(i1)%m
         endif
      enddo
   endif
end subroutine nearest_qdo
