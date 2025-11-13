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
module jastrow_factors_m
   use fortran_kinds_v, only: dp, int32, stdout
   use openmp_mpi_m, only: mpi_rank
   use write_lines_m
   use molecular_system_v, only: n_fe
   use fermionic_config_c, only: frm_cnf
   use quantum_monte_carlo_v, only: restart, grad_updat, n_wlk_max
   use jastrow_factors_v
   use cusps_functions_params_m, only: b_en, b_pn, b_ee, b_pp, b_ep, jcen_prs, jcpn_prs, jcee_prs, jcpp_prs, &
   & jcep_prs, n_par_jcen, n_par_jcpn, n_par_jcee, n_par_jcpp, n_par_jcep
   use jastrow_orbitals_m, only: n_jorbs, n_par_jorbs, jorbs_mat, n_par_jorbs_s, jc_opt, je_opt
   use jeporb_mod, only: n_jepos, n_par_jepos, jepos_mat, n_par_jepos_s, jpoc_opt, jpoe_opt
   use cusps_factors_m, only: init_cusps, comp_cusps, variation_1body_cusps, variation_2body_cusps, updt_cusps, &
   & DD2_cusps, new_DD2_cusps, updt_DD2_cusps, dp_cusps
   use orbitalsbased_jastrow_factors_m, only: init_orbitalsbased_jastrow_factors, &
   & comp_orbitalsbased_jastrow_factors, variation_orbitalsbased_1body_factors, &
   & variation_orbitalsbased_2body_factors, updt_orbitalsbased_jastrow_factors, &
   & DD2_orbitalsbased_jastrow_factors, new_DD2_orbitalsbased_jastrow_factors, &
   & updt_DD2_orbitalsbased_jastrow_factors, dp_orbitalsbased_jastrow_factors
   use jastrow_factor_params_m, only: read_jastrow_factor_params
   implicit none
   real(dp), allocatable, dimension(:)     :: jst
   real(dp), allocatable, dimension(:)     :: jst_var_1b
   real(dp), allocatable, dimension(:)     :: jst_var_2b
   real(dp), allocatable, dimension(:,:)   :: DD2_jst
   real(dp), allocatable, dimension(:,:)   :: new_DD2_jst
   real(dp), allocatable, dimension(:,:)   :: vec_1b
   real(dp), allocatable, dimension(:)     :: new_1b
   public :: init_jastrow_factors, comp_jastrow_factors, variation_jastrow_factors, &
   & variation_1body_jastrow_factors, variation_2body_jastrow_factors, updt_jastrow_factors, &
   & DD2_jastrow_factors, new_DD2_jastrow_factors, updt_DD2_jastrow_factors, &
   & dp_jastrow_factors
contains
   subroutine init_jastrow_factors( )
      jst_prs = .false.
      call write_simple_line(stdout,0,mpi_rank,2,"l","Initializing Jastrow factor.")
      call write_empty_line(stdout,0,mpi_rank)
      call init_cusps()
      call init_orbitalsbased_jastrow_factors()
      if ( restart ) call read_jastrow_factor_params()
      n_par_jst   = n_par_jstc + n_par_jstd   + n_par_jorbs
      n_par_jst_s = n_par_jstc + n_par_jstd_s + n_par_jorbs_s
      n_par_jst   = n_par_jst   + n_par_jepd + n_par_jepa + n_par_jepos
      n_par_jst_s = n_par_jst_s + n_par_jepd_s + n_par_jepa_s + n_par_jepos_s
      if ( jcen_prs .or. jcpn_prs .or. jc2_prs .or. jd1_prs .or. jd2_prs .or. jdep_prs.or.jaep_prs) jst_prs = .true.
      if( jst_prs ) then
#ifdef _PRNTALL
         call write_variable_line(stdout,0,mpi_rank,2,"Total number of Jastrow symmetrized parameters", n_par_jst_s,var_name="n_par_jst_s")
#endif
         allocate(jst(1:n_wlk_max) ) ; jst = 0.0_dp
         allocate( new_1b(1:n_wlk_max) ) ; new_1b = 0.0_dp
         allocate(jst_var_1b(1:n_wlk_max)) ; jst_var_1b = 0.0_dp
         allocate(jst_var_2b(1:n_wlk_max)) ; jst_var_2b = 0.0_dp
         allocate( vec_1b(1:n_fe,1:n_wlk_max) ) ; vec_1b = 0.0_dp
         allocate( DD2_jst(1:4*n_fe,1:n_wlk_max) ) ; DD2_jst = 0.0_dp
         allocate( new_DD2_jst(1:4,1:n_wlk_max) ) ; new_DD2_jst = 0.0_dp
      else
         call write_simple_line(stdout,0,mpi_rank,2,"l","WARNING!!! No Jastrow factor is included.")
      endif
   end subroutine init_jastrow_factors
   subroutine comp_jastrow_factors( iw )
      integer(int32),  intent(in)    :: iw
      jst(iw) = 0.0_dp
      vec_1b(:,iw) = 0.0_dp
      call comp_cusps( iw, vec_1b(:,iw), jst(iw) )
      call comp_orbitalsbased_jastrow_factors( iw, jst(iw) )
   end subroutine comp_jastrow_factors
   subroutine variation_jastrow_factors(  iw, g )
      integer(int32),  intent(in)  :: iw
      real(dp),        intent(out) :: g
      real(dp)                     :: g_j1b, g_j2b
      call variation_1body_jastrow_factors( iw, g_j1b )
      call variation_2body_jastrow_factors( iw, g_j2b )
      g = g_j1b * g_j2b
   end subroutine variation_jastrow_factors
   subroutine variation_1body_jastrow_factors( iw, g_j1b )
      integer(int32), intent(in)  :: iw
      real(dp),       intent(out) :: g_j1b
      jst_var_1b(iw) = 0.0_dp
      new_1b(iw) = 0.0_dp
      call variation_1body_cusps( iw, vec_1b(:,iw), new_1b(iw), jst_var_1b(iw) )
      call variation_orbitalsbased_1body_factors( iw, jst_var_1b(iw) )
      g_j1b = exp(jst_var_1b(iw))
   end subroutine variation_1body_jastrow_factors
   subroutine variation_2body_jastrow_factors( iw, g_j2b )
      integer(int32), intent(in)  :: iw
      real(dp),       intent(out) :: g_j2b
      jst_var_2b(iw) = 0.0_dp
      call variation_2body_cusps( iw, vec_1b(:,iw), new_1b(iw), jst_var_2b(iw) )
      call variation_orbitalsbased_2body_factors( iw, jst_var_2b(iw) )
      g_j2b = exp(jst_var_2b(iw))
   end subroutine variation_2body_jastrow_factors
   subroutine updt_jastrow_factors( iw )
      integer(int32), intent(in) :: iw
      jst(iw) = jst(iw) + jst_var_1b(iw) + jst_var_2b(iw)
      if( jcen_prs.or.jcpn_prs .or. jd1_prs ) vec_1b(frm_cnf(iw)%i_fe,iw) = new_1b(iw)
      call updt_cusps( iw )
      call updt_orbitalsbased_jastrow_factors( iw )
   end subroutine updt_jastrow_factors
   subroutine DD2_jastrow_factors( iw )
      integer(int32), intent(in) :: iw
      DD2_jst(:,iw) = 0.0_dp
      call DD2_cusps( iw, DD2_jst(:,iw) )
      call DD2_orbitalsbased_jastrow_factors( iw, DD2_jst(:,iw) )
   end subroutine DD2_jastrow_factors
   subroutine new_DD2_jastrow_factors( iw )
      integer(int32), intent(in) :: iw
      new_DD2_jst(:,iw) = 0.0_dp
      call new_DD2_cusps( iw, new_DD2_jst(:,iw) )
      call new_DD2_orbitalsbased_jastrow_factors( iw, new_DD2_jst(:,iw) )
   end subroutine new_DD2_jastrow_factors
   subroutine updt_DD2_jastrow_factors( iw )
      integer(int32), intent(in) :: iw
      DD2_jst(:,iw) = 0.0_dp
      if ( .not.grad_updat ) call new_DD2_jastrow_factors( iw )
      DD2_jst(4*frm_cnf(iw)%i_fe-3:4*frm_cnf(iw)%i_fe,iw) = new_DD2_jst(1:4,iw)
      call updt_DD2_cusps(iw, DD2_jst(:,iw) )
      call updt_DD2_orbitalsbased_jastrow_factors(iw, DD2_jst(:,iw) )
   end subroutine updt_DD2_jastrow_factors
   subroutine dp_jastrow_factors( iw, dp_ln_jst )
      integer(int32),                   intent(in) :: iw
      real(dp), dimension(n_par_jst), intent(inout) :: dp_ln_jst
      call dp_cusps( iw, dp_ln_jst(1:n_par_jstc) )
      if (n_par_jst.gt.n_par_jstc) call dp_orbitalsbased_jastrow_factors( iw, dp_ln_jst(n_par_jstc+1:n_par_jst) )
   end subroutine dp_jastrow_factors
end module jastrow_factors_m
