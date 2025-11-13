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
module drudonic_sampling_m
   use fortran_kinds_v,         only: int32, dp, stdout
   use openmp_mpi_m,            only: mpi_rank
   use write_lines_m
   use quantum_monte_carlo_v,   only: dt_d, var_dt, smp_type, a_drft
   use mersenne_twister19937_m, only: rndn
   use acceptance_rate_m
   implicit none
   public :: smpl_mtrpls_drd, smpl_lngvin_drd
contains
   subroutine smpl_mtrpls_drd( iw, move_acc, count )
      use qdo_wavefunction_m, only: wvfnq
      use drudonic_config_c, only: drd_cnf
      use qdo_system_v, only: qdos
      integer(int32), intent(in)  :: iw
      logical,        intent(in)  :: count
      logical,        intent(out) :: move_acc
      real(dp)                    :: g
      integer(int32)              :: i_dir
      real(dp), dimension(3)    :: dr
      real(dp)                    :: dt
      call drd_cnf(iw)%chs( iw )
      dt = dt_d / qdos(drd_cnf(iw)%i_drd)%qdo_m
      select case( smp_type )
       case(0,1,5)
         i_dir = int( 3.0_dp * rndn(iw)%rndo() + 1.0_dp )
      end select
      dr(:) = 0.0_dp
      select case( smp_type )
       case(0) ! Simple sampling 1d
         dr(i_dir) = dt * (2.0_dp * rndn(iw)%rndo() - 1.0_dp)
       case(1,5) ! Gaussian move 1d
         dr(i_dir) = sqrt(dt) * rndn(iw)%rndbm()
       case(2) ! Simple sampling 3d
         dr(1) = dt * (2.0_dp * rndn(iw)%rndo() - 1.0_dp)
         dr(2) = dt * (2.0_dp * rndn(iw)%rndo() - 1.0_dp)
         dr(3) = dt * (2.0_dp * rndn(iw)%rndo() - 1.0_dp)
       case(3,6) ! Gaussian move 3d
         dr(1) = sqrt(dt) * rndn(iw)%rndbm()
         dr(2) = sqrt(dt) * rndn(iw)%rndbm()
         dr(3) = sqrt(dt) * rndn(iw)%rndbm()
      end select
      call drd_cnf(iw)%new( iw, r_new=( drd_cnf(iw)%r_drd(:,drd_cnf(iw)%i_drd) + dr(1:3) ) )
      call wvfnq(iw)%ratio( iw, g ) ; g = g**2
      move_acc = .false.
      if(g.ge.1.0_dp) then
         move_acc = .true.
      else
         if( g.ge.rndn(iw)%rndo() ) move_acc = .true.
      endif
      if( count ) then
         if( move_acc ) then
            n_acc_mv_pw_d(iw) = n_acc_mv_pw_d(iw) + 1_int32
         endif
         n_tot_mv_pw_d(iw) = n_tot_mv_pw_d(iw) + 1_int32
      endif
   end subroutine smpl_mtrpls_drd
   subroutine smpl_lngvin_drd( iw, move_acc, count, dr2, pdr2 )
      use drudonic_config_c,  only: drd_cnf
      use qdo_wavefunction_m, only: wvfnq
      use qdo_system_v,       only: qdos, n_qdo
      integer(int32),           intent(in)    :: iw
      logical,                  intent(in)    :: count
      logical,                  intent(out)   :: move_acc
      real(dp),       optional, intent(inout) :: dr2, pdr2
      real(dp), dimension(6)    :: vec_v
      real(dp), dimension(3)    :: d_r, dr_diff
      real(dp)                    :: g, dt, a
      external                    :: rsc_drift_vel
      external                    :: smpl_drift_diff
      real(dp), external          :: gaus_func_3d, trns_prob
      call drd_cnf(iw)%chs( iw )
      dt = dt_d / qdos(drd_cnf(iw)%i_drd)%qdo_m
      vec_v(1:3) = wvfnq(iw)%DD2_ln_wvfnq( (4*drd_cnf(iw)%i_drd-3):(4*drd_cnf(iw)%i_drd-1) )
      a = a_drft
      call rsc_drift_vel( dt, a, vec_v(1:3) )
      call smpl_drift( dt, vec_v(1:3), drd_cnf(iw)%r_drd(:,drd_cnf(iw)%i_drd), d_r )
      call gaus_dist_3d( iw, sqrt(dt), dr_diff )
      call drd_cnf(iw)%new( iw, r_new=(d_r + dr_diff) )
      call wvfnq(iw)%ratio( iw, g )
      if ( g.ge.0.0_dp ) then
         g = g**2
         g = g / gaus_func_3d( dt, dr_diff )
         call wvfnq(iw)%drift( iw )
         vec_v(4:6) = wvfnq(iw)%DD2_ln_wvfnq_new(1:3)
         a = a_drft
         call rsc_drift_vel( dt, a, vec_v(4:6) )
         call smpl_drift( dt, vec_v(4:6), drd_cnf(iw)%r_drd_new, d_r )
         g = g * gaus_func_3d( dt, drd_cnf(iw)%r_drd(:,drd_cnf(iw)%i_drd) - d_r )
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
         a=sum(dr_diff(1:3)**2) * qdos(drd_cnf(iw)%i_drd)%qdo_m
         pdr2  = pdr2 + dt_d*g * a
         dr2   = dr2 + a
      endif
      if( count ) then
         if( move_acc ) n_acc_mv_pw_d(iw) = n_acc_mv_pw_d(iw) + 1_int32
         n_tot_mv_pw_d(iw) = n_tot_mv_pw_d(iw) + 1_int32
      endif
   end subroutine smpl_lngvin_drd
end module drudonic_sampling_m
