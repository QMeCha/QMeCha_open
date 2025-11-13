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
module qdo_wavefunction_c
   use fortran_kinds_v,   only: dp, int32
   use qdo_wavefunction_v,   only: ql_opt, wf_type_d, n_par_drd, n_par_drd_s, &
   & n_par_drd_wvfn, n_par_drd_wvfn_s
   use quantum_monte_carlo_v, only: grad_updat
   use qdo_system_v, only: n_qdo
   use prdqcfs_mod,  only: prdq_c_d
   use prdqfun_cls,  only: prdqfun_t
   use dipqcfs_mod,  only: dipq_c_d
   use dipqfun_cls,  only: dipqfun_t
   use drudonic_config_c,   only: drd_cnf
   use drudonic_orbitals_m,   only: dorbs_mat, n_dorbs, qc_opt, n_par_dorbs
   use jstqe_mod,    only: jstqe_fct
   use fermionic_wavefunction_m,   only: wvfn
   use fermionic_wavefunction_v,   only: n_par_frm_wvfn
   use qdo_jastrow_m, only: qdo_jastrow, qdo_cusp_prs, qdo_cusp_opt, n_qdo_jst_par
   use jstqepar_var,  only: jstqe_prs, jqe_opt, n_par_jstqe
   implicit none
   type, public :: wavfunq_t
      real(dp)                  :: wvfnq, logwvfnq, swvfnq
      type(prdqfun_t)              :: prdq_d
      type(dipqfun_t)              :: dipq_d
      real(dp)                     :: g_d
      real(dp)                     :: g_qe
      real(dp), allocatable, dimension(:) ::  DD2_ln_wvfnq
      real(dp), allocatable, dimension(:) ::  DD2_ln_wvfnq_new
      real(dp), allocatable, dimension(:) ::  dp_ln_wvfnq
   contains
      procedure :: ini    => ini_wvfnq
      procedure :: cmp    => cmp_wvfnq
      procedure :: ratio  => ratio_wvfnq
      procedure :: drift  => drift_wvfnq
      procedure :: upd    => upd_wvfnq
      procedure :: cmp_D  => cmp_DD2_wvfnq
      procedure :: upd_D  => upd_DD2_wvfnq
      procedure :: cmp_dp => cmp_dp_wvfnq
   end type wavfunq_t
   private :: ini_wvfnq, cmp_wvfnq, ratio_wvfnq, drift_wvfnq, upd_wvfnq, &
   & cmp_DD2_wvfnq, upd_DD2_wvfnq, cmp_dp_wvfnq
contains
   subroutine ini_wvfnq( obj )
      class(wavfunq_t), intent(inout) :: obj
      obj%wvfnq  = 1.0_dp
      obj%g_d    = 1.0_dp
      obj%g_qe   = 1.0_dp
      select case ( abs(wf_type_d) )
       case(5)
         call obj%prdq_d%ini( n_qdo )
       case(6)
         call obj%dipq_d%ini( n_qdo )
       case default
      end select
      allocate( obj%DD2_ln_wvfnq(1:4*n_qdo) )       ; obj%DD2_ln_wvfnq     = 0.0_dp
      if ( grad_updat ) then
         allocate( obj%DD2_ln_wvfnq_new(1:4) ) ; obj%DD2_ln_wvfnq_new = 0.0_dp
      endif
      if ( n_par_drd_wvfn.gt.0_int32 ) then
         allocate( obj%dp_ln_wvfnq(1:n_par_drd_wvfn) )   ; obj%dp_ln_wvfnq = 0.0_dp
      endif
   end subroutine ini_wvfnq
   subroutine cmp_wvfnq( obj, iw )
      class(wavfunq_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      obj%wvfnq  = 1.0_dp
      select case ( abs(wf_type_d) )
       case(5)
         if ( n_dorbs.ne.0_int32 ) call dorbs_mat(iw)%cmp()
         call obj%prdq_d%cmp( n_qdo, prdq_c_d, dorbs_mat(iw)%o )
         obj%wvfnq  = obj%wvfnq * obj%prdq_d%psi_prd
       case(6)
         call obj%dipq_d%cmp( n_qdo, dipq_c_d, drd_cnf(iw) )
         obj%wvfnq  = obj%wvfnq * obj%dipq_d%psi_dip
       case default
      end select
      obj%swvfnq   = sign(1.0_dp, obj%wvfnq)
      obj%logwvfnq = log(abs(obj%wvfnq))
      if ( jstqe_prs) then
         call jstqe_fct(iw)%cmp( iw )
         obj%logwvfnq = obj%logwvfnq + log(abs( jstqe_fct(iw)%jstqe ))
      endif
      if (qdo_cusp_prs) then
         call qdo_jastrow(iw)%cmp( )
         obj%logwvfnq = obj%logwvfnq + qdo_jastrow(iw)%jst
      endif
   end subroutine cmp_wvfnq
   subroutine ratio_wvfnq( obj, iw, g )
      class(wavfunq_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(out)   :: g
      real(dp) :: g_1b, g_2b
      if ( drd_cnf(iw)%i_type.eq.1_int32 ) then
         select case ( abs(wf_type_d) )
          case(5)
            if ( grad_updat ) then
               call dorbs_mat(iw)%new_D( )
            else
               call dorbs_mat(iw)%new( )
            endif
            call obj%prdq_d%ratio( drd_cnf(iw)%i_drd, prdq_c_d, dorbs_mat(iw)%o_new, obj%g_d )
          case(6)
            call obj%dipq_d%ratio( n_qdo, dipq_c_d, drd_cnf(iw), obj%g_d )
          case default
            obj%g_d = 1.0_dp
         end select
      else
         obj%g_d = 1.0_dp
      endif
      if( jstqe_prs ) then
         call jstqe_fct(iw)%ratio( iw, obj%g_qe )
      else
         obj%g_qe = 1.0_dp
      endif
      if (qdo_cusp_prs) then
         call qdo_jastrow(iw)%r1b( g_1b )
         call qdo_jastrow(iw)%r2b( g_2b )
      else
         g_1b = 1.0_dp
         g_2b = 1.0_dp
      endif
      g = obj%g_d * obj%g_qe * g_1b * g_2b
   end subroutine ratio_wvfnq
   subroutine drift_wvfnq( obj, iw )
      class(wavfunq_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      if ( drd_cnf(iw)%i_type.eq.1_int32 ) then
         select case ( abs(wf_type_d) )
          case(5)
            call obj%prdq_d%new_D( drd_cnf(iw)%i_drd, prdq_c_d, dorbs_mat(iw)%DD2_o_new )
            obj%DD2_ln_wvfnq_new = obj%prdq_d%DD2_ln_prdq_new
          case(6)
            call obj%dipq_d%new_D( drd_cnf(iw)%i_drd, dipq_c_d )
            obj%DD2_ln_wvfnq_new = obj%dipq_d%DD2_ln_dipq_new
          case default
         end select
      else
         obj%DD2_ln_wvfnq_new = 0.0_dp
      endif
      if ( jstqe_prs ) then
         call jstqe_fct(iw)%new_D( iw )
         obj%DD2_ln_wvfnq_new = obj%DD2_ln_wvfnq_new + jstqe_fct(iw)%new_DD2_jstqe
      endif
      if (qdo_cusp_prs) then
         call qdo_jastrow(iw)%nDD2( obj%DD2_ln_wvfnq_new )
      endif
   end subroutine drift_wvfnq
   subroutine upd_wvfnq( obj, iw )
      class(wavfunq_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      if ( drd_cnf(iw)%i_type.eq.1_int32 ) then
         if ( abs(wf_type_d).eq.5 ) then
            call dorbs_mat(iw)%upd()
         endif
         select case ( abs(wf_type_d) )
          case(5)
            call obj%prdq_d%upd( drd_cnf(iw)%i_drd )
          case(6)
            call obj%dipq_d%upd(  )
          case default
         end select ! wf_type_d
      endif
      obj%wvfnq = obj%wvfnq * obj%g_d
      obj%swvfnq   = sign(1.0_dp, obj%wvfnq)
      obj%logwvfnq = log(abs(obj%wvfnq))
      if ( jstqe_prs ) then
         call jstqe_fct(iw)%upd( iw )
         obj%logwvfnq = obj%logwvfnq + log(abs( jstqe_fct(iw)%jstqe ))
      endif
      if (qdo_cusp_prs) then
         call qdo_jastrow(iw)%upd( )
         obj%logwvfnq = obj%logwvfnq + qdo_jastrow(iw)%jst
      endif
   end subroutine upd_wvfnq
   subroutine cmp_DD2_wvfnq( obj, iw )
      class(wavfunq_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      select case ( abs(wf_type_d) )
       case(5)
         call obj%prdq_d%cmp_D( n_qdo, prdq_c_d, dorbs_mat(iw)%DD2_o )
         obj%DD2_ln_wvfnq(:) = obj%prdq_d%DD2_ln_prdq(:)
       case(6)
         call obj%dipq_d%cmp_D( n_qdo, dipq_c_d )
         obj%DD2_ln_wvfnq(:) = obj%dipq_d%DD2_ln_dipq(:)
       case default
      end select
      if ( jstqe_prs )  then
         call jstqe_fct(iw)%cmp_D()
         obj%DD2_ln_wvfnq = obj%DD2_ln_wvfnq + jstqe_fct(iw)%DD2d_jstqe
      endif
      if (qdo_cusp_prs) then
         call qdo_jastrow(iw)%cDD2( obj%DD2_ln_wvfnq )
      endif
   end subroutine cmp_DD2_wvfnq
   subroutine upd_DD2_wvfnq( obj, iw )
      class(wavfunq_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      if ( drd_cnf(iw)%i_type.eq.1_int32 ) then
         select case ( abs(wf_type_d) )
          case(5)
            call obj%prdq_d%upd_D( n_qdo, drd_cnf(iw)%i_drd, prdq_c_d, dorbs_mat(iw)%DD2_o )
            obj%DD2_ln_wvfnq(:) = obj%prdq_d%DD2_ln_prdq(:)
          case(6)
            call obj%dipq_d%upd_D(n_qdo )
            obj%DD2_ln_wvfnq(:) = obj%dipq_d%DD2_ln_dipq(:)
          case default
         end select
      else
         select case ( abs(wf_type_d) )
          case(5)
            obj%DD2_ln_wvfnq(:) = obj%prdq_d%DD2_ln_prdq(:)
          case(6)
            obj%DD2_ln_wvfnq(:) = obj%dipq_d%DD2_ln_dipq(:)
          case default
         end select
      endif
      if ( jstqe_prs )  then
         call jstqe_fct(iw)%upd_D( iw )
         obj%DD2_ln_wvfnq = obj%DD2_ln_wvfnq + jstqe_fct(iw)%DD2d_jstqe
      endif
      if (qdo_cusp_prs) then
         call qdo_jastrow(iw)%uDD2( obj%DD2_ln_wvfnq )
      endif
   end subroutine upd_DD2_wvfnq
   subroutine cmp_dp_wvfnq( obj, iw )
      class(wavfunq_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      integer(int32) :: ip
      ip = 1_int32
      select case ( abs(wf_type_d) )
       case(5)
         if ( ql_opt ) then
            call obj%prdq_d%cmp_dl( n_qdo, dorbs_mat(iw)%o, n_par_drd, &
            & obj%dp_ln_wvfnq(ip:ip+n_par_drd-1)  )
            ip = ip + n_par_drd
         endif
         if ( qc_opt ) then
            call obj%prdq_d%cmp_da( n_qdo, prdq_c_d, dorbs_mat(iw)%da_o, &
            & obj%dp_ln_wvfnq(ip:ip+n_par_dorbs-1) )
            ip = ip + n_par_dorbs
         endif
       case(6)
         if ( ql_opt ) then
            call obj%dipq_d%cmp_dl( n_qdo, n_par_drd, obj%dp_ln_wvfnq(1:n_par_drd) )
            ip = ip + n_par_drd
         endif
       case default
      end select
      if( jstqe_prs.and.jqe_opt ) then
         call jstqe_fct(iw)%cmp_dp( obj%dp_ln_wvfnq(ip:ip+n_par_jstqe-1 ) )
         ip = ip + n_par_jstqe
      endif
      if (qdo_cusp_opt) then
         call qdo_jastrow(iw)%dpar( obj%dp_ln_wvfnq(ip:ip+n_qdo_jst_par-1 ) )
         ip = ip + n_qdo_jst_par
      endif
   end subroutine cmp_dp_wvfnq
end module qdo_wavefunction_c
