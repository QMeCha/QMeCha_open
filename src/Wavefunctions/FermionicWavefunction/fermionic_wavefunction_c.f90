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
module fermionic_wavefunction_c
   use fortran_kinds_v, only: dp, int32
   use fermionic_wavefunction_v
   use jastrow_factors_v
   use quantum_monte_carlo_v, only: grad_updat, n_wlk_max
   use molecular_system_v, only: n_el, n_el_u, n_el_d, n_el_s, n_po, n_po_u, n_po_d, n_po_s, n_fe
   use fermionic_config_c, only: frm_cnf
   use slater_determinant_params_m, only: sld_c_e, sld_c_p
   use singlet_geminal_params_m, only: sgm_c_e, sgm_c_p
   use triplet_geminal_params_m, only: tgm_c_e, tgm_c_p
   use pfaffian_params_m, only: pff_c_e, pff_c_p
   use eppcfs_mod, only: epo_c
   use fermionic_orbitals_m, only: aoc_opt, forbs_mat, n_forbs, n_par_forbs, forbs_s
   use jastrow_orbitals_m, only: jorbs_mat, jorbs_s
   use feporb_mod, only: poc_opt, eporbs_mat, n_par_eporbs, n_eporbs
   use jeporb_mod, only: jepos_mat, n_par_jepos
   use jastrow_factors_m
   use slater_determinant_c, only: sldfun_t
   use singlet_geminal_c, only: sgmfun_t
   use triplet_geminal_c, only: tgmfun_t
   use pfaffian_c, only: pfffun_t
   use eppfun_cls, only: eppfun_t
   implicit none
   type, public :: wavfun_t
      real(dp), allocatable, dimension(:) :: wvfn, logwvfn, swvfn
      type(sldfun_t)           :: sld_e
      type(sldfun_t)           :: sld_p
      type(sgmfun_t)           :: sgm_e
      type(sgmfun_t)           :: sgm_p
      type(tgmfun_t)           :: tgm_e
      type(tgmfun_t)           :: tgm_p
      type(pfffun_t)           :: pff_e
      type(pfffun_t)           :: pff_p
      type(eppfun_t)           :: epp
      real(dp), allocatable, dimension(:) :: g_f, g_ep
      real(dp), allocatable, dimension(:) :: g_j1b, g_j2b
      real(dp), allocatable, dimension(:,:)   :: DD2_ln_wvfn
      real(dp), allocatable, dimension(:,:)   :: DD2_ln_wvfn_new
      real(dp), allocatable, dimension(:,:)   :: dp_ln_wvfn
   contains
      procedure :: ini    => ini_wvfn
      procedure :: cmp    => cmp_wvfn
      procedure :: ratio  => ratio_wvfn
      procedure :: ratde  => ratio_wvde
      procedure :: ratla  => ratio_wvla
      procedure :: rat1b  => ratio_wv1b
      procedure :: rat2b  => ratio_wv2b
      procedure :: drift  => drift_wvfn
      procedure :: upd    => upd_wvfn
      procedure :: cmp_D  => cmp_DD2_wvfn
      procedure :: upd_D  => upd_DD2_wvfn
      procedure :: cmp_dp => cmp_dp_wvfn
   end type wavfun_t
   private :: ini_wvfn, cmp_wvfn, ratio_wvfn, upd_wvfn, drift_wvfn, cmp_DD2_wvfn, &
   & upd_DD2_wvfn, cmp_dp_wvfn, ratio_wvde, ratio_wvla, ratio_wv1b, ratio_wv2b
contains
   subroutine ini_wvfn( obj )
      class(wavfun_t), intent(inout) :: obj
      allocate(obj%wvfn(1:n_wlk_max) ) ; obj%wvfn = 0.0_dp
      allocate(obj%logwvfn(1:n_wlk_max) ) ; obj%logwvfn = 0.0_dp
      allocate(obj%swvfn(1:n_wlk_max) ) ; obj%swvfn = 1.0_dp
      allocate( obj%g_f(1:n_wlk_max) ) ; obj%g_f = 1.0_dp
      allocate( obj%g_ep(1:n_wlk_max) ) ; obj%g_ep = 1.0_dp
      allocate( obj%g_j1b(1:n_wlk_max) ) ; obj%g_j1b = 1.0_dp
      allocate( obj%g_j2b(1:n_wlk_max) ) ; obj%g_j2b = 1.0_dp
      select case ( abs(wf_type_e) )
       case(1)
         call obj%sld_e%ini( n_el, n_el_u, n_el_d, edl_opt )
       case(2)
         call obj%sgm_e%ini( n_el, n_el_u, edl_opt )
       case(3)
         call obj%tgm_e%ini( n_el, n_el_u, n_el_d, edl_opt, tgm_c_e )
       case(4)
         call obj%pff_e%ini( n_el, n_el_u, n_el_d, edl_opt )
       case default
      end select
      if ( wf_type_p.ne.0 ) then
         select case ( abs(wf_type_p) )
          case(1)
            call obj%sld_p%ini( n_po, n_po_u, n_po_d, pdl_opt )
          case(2)
            call obj%sgm_p%ini( n_po, n_po_u, pdl_opt )
          case(3)
            call obj%tgm_p%ini( n_po, n_po_u, n_po_d, pdl_opt, tgm_c_p )
          case(4)
            call obj%pff_p%ini( n_po, n_po_u, n_po_d, pdl_opt )
          case default
         end select 
      endif
      if ( wf_type_ep.ne.0 ) then
         select case ( abs(wf_type_ep) )
          case(11)
            call obj%epp%ini(  )
          case default
         end select
      endif
      allocate( obj%DD2_ln_wvfn(1:4*n_fe,1:n_wlk_max) ) ; obj%DD2_ln_wvfn = 0.0_dp
      if ( grad_updat ) then
         allocate( obj%DD2_ln_wvfn_new(1:4,1:n_wlk_max) ) ; obj%DD2_ln_wvfn_new = 0.0_dp
      endif
      if ( n_par_frm_wvfn.gt.0_int32 ) then
         allocate( obj%dp_ln_wvfn(1:n_par_frm_wvfn,1:n_wlk_max)) ; obj%dp_ln_wvfn = 0.0_dp
      endif
   end subroutine ini_wvfn
   subroutine cmp_wvfn( obj, iw )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      if (n_forbs.ne.0) call forbs_mat(iw)%cmp( forbs_s)
      if ( jd1_prs .or. jd2_prs ) call jorbs_mat(iw)%cmp( jorbs_s)
      if (n_eporbs.ne.0) call eporbs_mat(iw)%cmp( )
      if (jdep_prs) call jepos_mat(iw)%cmp()
      obj%wvfn(iw) = 1.0_dp
      select case ( abs(wf_type_e) )
       case(1)
         call obj%sld_e%cmp( iw, n_el, n_el_u, n_el_d, sld_c_e, forbs_mat(iw)%o(:,1:n_el) )
         obj%wvfn(iw) = obj%wvfn(iw) * obj%sld_e%det_S_u(iw) * obj%sld_e%det_S_d(iw)
       case(2)
         call obj%sgm_e%cmp( iw, n_el, n_el_u, n_el_d, n_el_s, sgm_c_e, forbs_mat(iw)%o(:,1:n_el) )
         obj%wvfn(iw) = obj%wvfn(iw) * obj%sgm_e%det_G(iw)
       case(3)
         call obj%tgm_e%cmp( iw, n_el, n_el_u, n_el_d, tgm_c_e, forbs_mat(iw)%o(:,1:n_el) )
         obj%wvfn(iw) = obj%wvfn(iw) * obj%tgm_e%pff_T_u(iw) * obj%tgm_e%pff_T_d(iw)
       case(4)
         call obj%pff_e%cmp( iw, n_el, n_el_u, n_el_d, pff_c_e, forbs_mat(iw)%o(:,1:n_el) )
         obj%wvfn(iw) = obj%wvfn(iw) * obj%pff_e%pff_P(iw)
       case default
      end select 
      if (n_po.gt.0_int32) then
         select case ( abs(wf_type_p) )
          case(1)
            call obj%sld_p%cmp( iw, n_po, n_po_u, n_po_d, sld_c_p, forbs_mat(iw)%o(:,n_el+1:n_fe) )
            obj%wvfn(iw) = obj%wvfn(iw) * obj%sld_p%det_S_u(iw) * obj%sld_p%det_S_d(iw)
          case(2)
            call obj%sgm_p%cmp( iw, n_po, n_po_u, n_po_d, n_po_s, sgm_c_p, forbs_mat(iw)%o(:,n_el+1:n_fe) )
            obj%wvfn(iw) = obj%wvfn(iw) * obj%sgm_p%det_G(iw)
          case(3)
            call obj%tgm_p%cmp( iw, n_po, n_po_u, n_po_d, tgm_c_p, forbs_mat(iw)%o(:,n_el+1:n_fe) )
            obj%wvfn(iw) = obj%wvfn(iw) * obj%tgm_p%pff_T_u(iw) * obj%tgm_p%pff_T_d(iw)
          case(4)
            call obj%pff_p%cmp( iw, n_po, n_po_u, n_po_d, pff_c_p, forbs_mat(iw)%o(:,n_el+1:n_fe) )
            obj%wvfn(iw) = obj%wvfn(iw) * obj%pff_p%pff_P(iw)
          case default
         end select 
      endif
      if (wf_type_ep.ne.0) then
         select case ( abs(wf_type_ep) )
          case(11)
            call obj%epp%cmp( iw, n_el, n_po, epo_c, eporbs_mat(iw)%o )
            obj%wvfn(iw) = obj%wvfn(iw) * obj%epp%gem(iw)
          case default
         end select
      endif
      obj%swvfn(iw)   = sign(1.0_dp, obj%wvfn(iw))
      obj%logwvfn(iw) = log(abs(obj%wvfn(iw)))
      if(jst_prs) then
         call comp_jastrow_factors( iw )
         obj%logwvfn(iw) = obj%logwvfn(iw) + jst(iw)
      endif
   end subroutine cmp_wvfn
   subroutine ratio_wvfn( obj, iw, g )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(out)   :: g
      real(dp)                       :: g_1b, g_2b
      call obj%rat1b( iw, g_1b )
      call obj%rat2b( iw, g_2b )
      g = g_1b * g_2b
   end subroutine ratio_wvfn
   subroutine ratio_wv1b( obj, iw, g_1b )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(out)   :: g_1b
      real(dp)                       :: g_de
      call obj%ratde( iw, g_de )
      if(jst_prs) then
         call variation_1body_jastrow_factors( iw, obj%g_j1b(iw) )
      else
         obj%g_j1b(iw) = 1.0_dp
      endif
      g_1b = obj%g_j1b(iw) * g_de
   end subroutine ratio_wv1b
   subroutine ratio_wvde( obj, iw, g_de, pseudo_update )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(out)   :: g_de
      logical, optional, intent(in)  :: pseudo_update
      if (.not.present(pseudo_update)) then
         if ( grad_updat ) then
            if (n_forbs.ne.0) call forbs_mat(iw)%new_D(forbs_s )
            if( jd1_prs .or. jd2_prs ) call jorbs_mat(iw)%new_D(jorbs_s)
            if (n_eporbs.ne.0) call eporbs_mat(iw)%new_D()
            if (jdep_prs) call jepos_mat(iw)%new_D()
         else
            if (n_forbs.ne.0) call forbs_mat(iw)%new(forbs_s )
            if( jd1_prs .or. jd2_prs ) call jorbs_mat(iw)%new(jorbs_s)
            if (n_eporbs.ne.0) call eporbs_mat(iw)%new()
            if (jdep_prs) call jepos_mat(iw)%new()
         endif
      else
         if (pseudo_update) then
            if (n_forbs.ne.0) call forbs_mat(iw)%new(forbs_s)
            if( jd1_prs .or. jd2_prs ) call jorbs_mat(iw)%new(jorbs_s)
            if (n_eporbs.ne.0) call eporbs_mat(iw)%new()
            if (jdep_prs) call jepos_mat(iw)%new()
         endif
      endif
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         select case ( abs(wf_type_e) )
          case(1)
            call obj%sld_e%ratio( iw, frm_cnf(iw)%i_fe, n_el_u, n_el_d, sld_c_e, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(2)
            call obj%sgm_e%ratio( iw, frm_cnf(iw)%i_fe, n_el_u, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(3)
            call obj%tgm_e%ratio( iw, frm_cnf(iw)%i_fe, n_el_u, n_el_d, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(4)
            call obj%pff_e%ratio( iw, frm_cnf(iw)%i_fe, n_el_u, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case default
            obj%g_f(iw) = 1.0_dp
         end select 
      else
         select case ( abs(wf_type_p) )
          case(1)
            call obj%sld_p%ratio( iw, frm_cnf(iw)%i_fe-n_el, n_po_u, n_po_d, sld_c_p, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(2)
            call obj%sgm_p%ratio( iw, frm_cnf(iw)%i_fe-n_el, n_po_u, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(3)
            call obj%tgm_p%ratio( iw, frm_cnf(iw)%i_fe-n_el, n_po_u, n_po_d, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(4)
            call obj%pff_p%ratio( iw, frm_cnf(iw)%i_fe-n_el, n_po_u, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case default
            obj%g_f(iw) = 1.0_dp
         end select
      endif
      if (wf_type_ep.ne.0) then
         select case ( abs(wf_type_ep) )
          case(11)
            call obj%epp%ratio( iw, frm_cnf(iw)%i_fe, n_el, n_po, epo_c, eporbs_mat(iw)%o_dif, obj%g_ep(iw) )
          case default
         end select
      else
         obj%g_ep(iw) = 1.0_dp
      endif
      g_de = obj%g_f(iw) * obj%g_ep(iw)
   end subroutine ratio_wvde
   subroutine ratio_wvla( obj, iw, g, pseudo_update )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(out)   :: g
      logical, optional, intent(in)  :: pseudo_update
      if (.not.present(pseudo_update)) then
         if ( grad_updat ) then
            if (n_forbs.ne.0) call forbs_mat(iw)%new_D( forbs_s)
            if( jd1_prs .or. jd2_prs ) call jorbs_mat(iw)%new_D(jorbs_s)
            if (n_eporbs.ne.0) call eporbs_mat(iw)%new_D()
            if (jdep_prs) call jepos_mat(iw)%new_D()
         else
            if (n_forbs.ne.0) call forbs_mat(iw)%new(forbs_s )
            if( jd1_prs .or. jd2_prs ) call jorbs_mat(iw)%new(jorbs_s)
            if (n_eporbs.ne.0) call eporbs_mat(iw)%new()
            if (jdep_prs) call jepos_mat(iw)%new()
         endif
      else
         if (pseudo_update) then
            if (n_forbs.ne.0) call forbs_mat(iw)%new(forbs_s)
            if( jd1_prs .or. jd2_prs ) call jorbs_mat(iw)%new(jorbs_s)
            if (n_eporbs.ne.0) call eporbs_mat(iw)%new()
            if (jdep_prs) call jepos_mat(iw)%new()
         endif
      endif
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         select case ( abs(wf_type_e) )
          case(1)
            call obj%sld_e%ratio( iw, frm_cnf(iw)%i_fe, n_el_u, n_el_d, sld_c_e, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(2)
            call obj%sgm_e%ratio( iw, frm_cnf(iw)%i_fe, n_el_u, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(3)
            call obj%tgm_e%ratio( iw, frm_cnf(iw)%i_fe, n_el_u, n_el_d, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(4)
            call obj%pff_e%ratio( iw, frm_cnf(iw)%i_fe, n_el_u, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case default
            obj%g_f(iw) = 1.0_dp
         end select 
      else
         select case ( abs(wf_type_p) )
          case(1)
            call obj%sld_p%ratio( iw, frm_cnf(iw)%i_fe-n_el, n_po_u, n_po_d, sld_c_p, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(2)
            call obj%sgm_p%ratio( iw, frm_cnf(iw)%i_fe-n_el, n_po_u, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(3)
            call obj%tgm_p%ratio( iw, frm_cnf(iw)%i_fe-n_el, n_po_u, n_po_d, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case(4)
            call obj%pff_p%ratio( iw, frm_cnf(iw)%i_fe-n_el, n_po_u, forbs_mat(iw)%o_new(:), obj%g_f(iw) )
          case default
            obj%g_f(iw) = 1.0_dp
         end select 
      endif
      if (wf_type_ep.ne.0) then
         select case ( abs(wf_type_ep) )
          case(11)
            call obj%epp%ratio( iw, frm_cnf(iw)%i_fe, n_el, n_po, epo_c, eporbs_mat(iw)%o_dif, obj%g_ep(iw) )
          case default
         end select
      else
         obj%g_ep(iw) = 1.0_dp
      endif
      if(jst_prs) then
         call variation_1body_jastrow_factors( iw, obj%g_j1b(iw) )
         call variation_2body_jastrow_factors( iw, obj%g_j2b(iw) )
      else
         obj%g_j1b(iw) = 1.0_dp
         obj%g_j2b(iw) = 1.0_dp
      endif
      g = obj%g_j1b(iw) * obj%g_j2b(iw) *obj%g_f(iw) * obj%g_ep(iw)
   end subroutine ratio_wvla
   subroutine ratio_wv2b( obj, iw, g_2b )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(out)   :: g_2b
      if(jst_prs) then
         call variation_2body_jastrow_factors( iw, obj%g_j2b(iw) )
      else
         obj%g_j2b(iw) = 1.0_dp
      endif
      g_2b = obj%g_j2b(iw)
   end subroutine ratio_wv2b
   subroutine drift_wvfn( obj, iw )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         select case ( abs(wf_type_e) )
          case(1)
            call obj%sld_e%new_D( iw, frm_cnf(iw)%i_fe, forbs_mat(iw)%DD2_o_new(:,:) )
            obj%DD2_ln_wvfn_new(:,iw) = obj%sld_e%DD2_ln_det_S_new(:,iw)
          case(2)
            call obj%sgm_e%new_D( iw, frm_cnf(iw)%i_fe, forbs_mat(iw)%DD2_o_new(:,:)  )
            obj%DD2_ln_wvfn_new(:,iw) = obj%sgm_e%DD2_ln_det_G_new(:,iw)
          case(3)
            call obj%tgm_e%new_D( iw, frm_cnf(iw)%i_fe, n_el_u, forbs_mat(iw)%DD2_o_new(:,:) )
            obj%DD2_ln_wvfn_new(:,iw) = obj%tgm_e%DD2_ln_pff_T_new(:,iw)
          case(4)
            call obj%pff_e%new_D( iw, frm_cnf(iw)%i_fe, forbs_mat(iw)%DD2_o_new(:,:)  )
            obj%DD2_ln_wvfn_new(:,iw) = obj%pff_e%DD2_ln_pff_P_new(:,iw)
          case default
            obj%DD2_ln_wvfn_new(:,iw) = 0.0_dp
         end select 
      else
         select case ( abs(wf_type_p) )
          case(1)
            call obj%sld_p%new_D( iw, frm_cnf(iw)%i_fe-n_el, forbs_mat(iw)%DD2_o_new(:,:) )
            obj%DD2_ln_wvfn_new(:,iw) = obj%sld_p%DD2_ln_det_S_new(:,iw)
          case(2)
            call obj%sgm_p%new_D( iw, frm_cnf(iw)%i_fe-n_el, forbs_mat(iw)%DD2_o_new(:,:)  )
            obj%DD2_ln_wvfn_new(:,iw) = obj%sgm_p%DD2_ln_det_G_new(:,iw)
          case(3)
            call obj%tgm_p%new_D( iw, frm_cnf(iw)%i_fe-n_el, n_po_u, forbs_mat(iw)%DD2_o_new(:,:) )
            obj%DD2_ln_wvfn_new(:,iw) = obj%tgm_p%DD2_ln_pff_T_new(:,iw)
          case(4)
            call obj%pff_p%new_D(iw, frm_cnf(iw)%i_fe-n_el, forbs_mat(iw)%DD2_o_new(:,:)  )
            obj%DD2_ln_wvfn_new(:,iw) = obj%pff_p%DD2_ln_pff_P_new(:,iw)
          case default
            obj%DD2_ln_wvfn_new(:,iw) = 0.0_dp
         end select 
      endif
      if (wf_type_ep.ne.0) then
         select case ( abs(wf_type_ep) )
          case(11)
            call obj%epp%new_D( iw, frm_cnf(iw)%i_fe, n_el, n_po, epo_c, eporbs_mat(iw)%DD2_o_new )
            obj%DD2_ln_wvfn_new(:,iw) = obj%DD2_ln_wvfn_new(:,iw) + obj%epp%DD2_ln_gem_new(:,iw)
          case default
         end select
      endif
      if ( jst_prs ) then
         call new_DD2_jastrow_factors( iw )
         obj%DD2_ln_wvfn_new(:,iw) = obj%DD2_ln_wvfn_new(:,iw) + new_DD2_jst(:,iw)
      endif
   end subroutine drift_wvfn
   subroutine upd_wvfn( obj, iw )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      if (n_forbs.ne.0) call forbs_mat(iw)%upd( forbs_s )
      if( jd1_prs .or. jd2_prs ) call jorbs_mat(iw)%upd(jorbs_s )
      if (n_eporbs.ne.0) call eporbs_mat(iw)%upd( )
      if( jdep_prs ) call jepos_mat(iw)%upd( )
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         select case ( abs(wf_type_e) )
          case(1)
            call obj%sld_e%upd( iw, frm_cnf(iw)%i_fe, n_el_u, n_el_d )
          case(2)
            call obj%sgm_e%upd( iw, frm_cnf(iw)%i_fe, n_el, n_el_u, sgm_c_e, forbs_mat(iw)%o_new(:) )
          case(3)
            call obj%tgm_e%upd( iw, frm_cnf(iw)%i_fe, n_el, n_el_u, n_el_d, tgm_c_e, forbs_mat(iw)%o(:,1:n_el) )
          case(4)
            call obj%pff_e%upd( iw, frm_cnf(iw)%i_fe, n_el, n_el_u, n_el_d, pff_c_e, forbs_mat(iw)%o(:,1:n_el) )
         end select 
      else
         select case ( abs(wf_type_p) )
          case(1)
            call obj%sld_p%upd( iw, frm_cnf(iw)%i_fe-n_el, n_po_u, n_po_d )
          case(2)
            call obj%sgm_p%upd( iw, frm_cnf(iw)%i_fe-n_el, n_po, n_po_u, sgm_c_p, forbs_mat(iw)%o_new(:) )
          case(3)
            call obj%tgm_p%upd( iw, frm_cnf(iw)%i_fe-n_el, n_po, n_po_u, n_po_d, tgm_c_p, forbs_mat(iw)%o(:,n_el+1:n_fe) )
          case(4)
            call obj%pff_p%upd( iw, frm_cnf(iw)%i_fe-n_el, n_po, n_po_u, n_po_d, pff_c_p, forbs_mat(iw)%o(:,n_el+1:n_fe) )
          case default
         end select 
      endif
      if (wf_type_ep.ne.0) then
         select case ( abs(wf_type_ep) )
          case(11)
            call obj%epp%upd( iw, frm_cnf(iw)%i_fe, n_el, n_po )
          case default
         end select
      endif
      obj%wvfn(iw) = obj%wvfn(iw) * obj%g_f(iw) * obj%g_ep(iw)
      obj%swvfn(iw)   = sign(1.0_dp, obj%wvfn(iw))
      obj%logwvfn(iw) = log(abs(obj%wvfn(iw)))
      if(jst_prs) then
         call updt_jastrow_factors(iw)
         obj%logwvfn(iw) = obj%logwvfn(iw) + jst(iw)
      endif
   end subroutine upd_wvfn
   subroutine cmp_DD2_wvfn( obj, iw )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      select case ( abs(wf_type_e) )
       case(1)
         call obj%sld_e%cmp_D( iw, n_el, n_el_u, n_el_d, sld_c_e, forbs_mat(iw)%DD2_o(:,1:4*n_el) )
         obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%sld_e%DD2_ln_det_S(1:4*n_el,iw)
       case(2)
         call obj%sgm_e%cmp_D( iw, n_el, n_el_u, forbs_mat(iw)%DD2_o(:,1:4*n_el) )
         obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%sgm_e%DD2_ln_det_G(1:4*n_el,iw)
       case(3)
         call obj%tgm_e%cmp_D( iw, n_el, n_el_u, n_el_d, forbs_mat(iw)%DD2_o(:,1:4*n_el) )
         obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%tgm_e%DD2_ln_pff_T(1:4*n_el,iw)
       case(4)
         call obj%pff_e%cmp_D( iw, n_el, n_el_u, n_el_d, forbs_mat(iw)%DD2_o(:,1:4*n_el) )
         obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%pff_e%DD2_ln_pff_P(1:4*n_el,iw)
       case default
         obj%DD2_ln_wvfn(1:4*n_el,iw) = 0.0_dp
      end select 
      if (n_po.gt.0_int32) then
         select case ( abs(wf_type_p) )
          case(1)
            call obj%sld_p%cmp_D( iw, n_po, n_po_u, n_po_d, sld_c_p, forbs_mat(iw)%DD2_o(:,4*n_el+1:4*n_fe) )
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%sld_p%DD2_ln_det_S(1:4*n_po,iw)
          case(2)
            call obj%sgm_p%cmp_D( iw, n_po, n_po_u, forbs_mat(iw)%DD2_o(:,4*n_el+1:4*n_fe) )
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%sgm_p%DD2_ln_det_G(1:4*n_po,iw)
          case(3)
            call obj%tgm_p%cmp_D( iw, n_po, n_po_u, n_po_d, forbs_mat(iw)%DD2_o(:,4*n_el+1:4*n_fe) )
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%tgm_p%DD2_ln_pff_T(1:4*n_po,iw)
          case(4)
            call obj%pff_p%cmp_D( iw, n_po, n_po_u, n_po_d, forbs_mat(iw)%DD2_o(:,4*n_el+1:4*n_fe) )
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%pff_p%DD2_ln_pff_P(1:4*n_po,iw)
          case default
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = 0.0_dp
         end select 
      endif
      if (wf_type_ep.ne.0) then
         select case ( abs(wf_type_ep) )
          case(11)
            call obj%epp%cmp_D( iw, n_el, n_po, epo_c, eporbs_mat(iw)%DD2_o  )
            obj%DD2_ln_wvfn(:,iw) = obj%DD2_ln_wvfn(:,iw) + obj%epp%DD2_ln_gem(:,iw)
          case default
         end select
      endif
      if (jst_prs) then
         call DD2_jastrow_factors( iw )
         obj%DD2_ln_wvfn(:,iw) = obj%DD2_ln_wvfn(:,iw) + DD2_jst(:,iw)
      endif
   end subroutine cmp_DD2_wvfn
   subroutine upd_DD2_wvfn( obj, iw )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         select case ( abs(wf_type_e) )
          case(1)
            call obj%sld_e%upd_D( iw, frm_cnf(iw)%i_fe, n_el, n_el_u, n_el_d, sld_c_e, forbs_mat(iw)%DD2_o(:,1:4*n_el) )
            obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%sld_e%DD2_ln_det_S(1:4*n_el,iw)
          case(2)
            call obj%sgm_e%upd_D( iw, frm_cnf(iw)%i_fe, n_el, n_el_u, forbs_mat(iw)%DD2_o(:,1:4*n_el) )
            obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%sgm_e%DD2_ln_det_G(1:4*n_el,iw)
          case(3)
            call obj%tgm_e%upd_D( iw, frm_cnf(iw)%i_fe, n_el, n_el_u, n_el_d, forbs_mat(iw)%DD2_o(:,1:4*n_el) )
            obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%tgm_e%DD2_ln_pff_T(1:4*n_el,iw)
          case(4)
            call obj%pff_e%upd_D( iw, n_el, n_el_u, n_el_d, forbs_mat(iw)%DD2_o(:,1:4*n_el) )
            obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%pff_e%DD2_ln_pff_P(1:4*n_el,iw)
          case default
            obj%DD2_ln_wvfn(1:4*n_el,iw) = 0.0_dp
         end select 
         select case ( abs(wf_type_p) )
          case(1)
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%sld_p%DD2_ln_det_S(1:4*n_po,iw)
          case(2)
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%sgm_p%DD2_ln_det_G(1:4*n_po,iw)
          case(3)
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%tgm_p%DD2_ln_pff_T(1:4*n_po,iw)
          case(4)
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%pff_p%DD2_ln_pff_P(1:4*n_po,iw)
          case default
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = 0.0_dp
         end select 
      else
         select case ( abs(wf_type_e) )
          case(1)
            obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%sld_e%DD2_ln_det_S(1:4*n_el,iw)
          case(2)
            obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%sgm_e%DD2_ln_det_G(1:4*n_el,iw)
          case(3)
            obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%tgm_e%DD2_ln_pff_T(1:4*n_el,iw)
          case(4)
            obj%DD2_ln_wvfn(1:4*n_el,iw) = obj%pff_e%DD2_ln_pff_P(1:4*n_el,iw)
          case default
            obj%DD2_ln_wvfn(1:4*n_el,iw) = 0.0_dp
         end select 
         select case ( abs(wf_type_p) )
          case(1)
            call obj%sld_p%upd_D( iw, frm_cnf(iw)%i_fe-n_el, n_po, n_po_u, n_po_d, sld_c_p, forbs_mat(iw)%DD2_o(:,4*n_el+1:4*n_fe) )
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%sld_p%DD2_ln_det_S(1:4*n_po,iw)
          case(2)
            call obj%sgm_p%upd_D( iw, frm_cnf(iw)%i_fe-n_el, n_po, n_po_u, forbs_mat(iw)%DD2_o(:,4*n_el+1:4*n_fe) )
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%sgm_p%DD2_ln_det_G(1:4*n_po,iw)
          case(3)
            call obj%tgm_p%upd_D( iw, frm_cnf(iw)%i_fe-n_el, n_po, n_po_u, n_po_d, forbs_mat(iw)%DD2_o(:,4*n_el+1:4*n_fe) )
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%tgm_p%DD2_ln_pff_T(1:4*n_po,iw)
          case(4)
            call obj%pff_p%upd_D( iw, n_po, n_po_u, n_po_d, forbs_mat(iw)%DD2_o(:,4*n_el+1:4*n_fe) )
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = obj%pff_p%DD2_ln_pff_P(1:4*n_po,iw)
          case default
            obj%DD2_ln_wvfn(4*n_el+1:4*n_fe,iw) = 0.0_dp
         end select 
      endif
      if (wf_type_ep.ne.0) then
         select case ( abs(wf_type_ep) )
          case(11)
            call obj%epp%upd_D( iw, frm_cnf(iw)%i_fe, n_el, n_po, epo_c, eporbs_mat(iw)%DD2_o  )
            obj%DD2_ln_wvfn(:,iw) = obj%DD2_ln_wvfn(:,iw) + obj%epp%DD2_ln_gem(:,iw)
          case default
         end select
      endif
      if (jst_prs) then
         call updt_DD2_jastrow_factors( iw )
         obj%DD2_ln_wvfn(:,iw) = obj%DD2_ln_wvfn(:,iw) + DD2_jst(:,iw)
      endif
   end subroutine upd_DD2_wvfn
   subroutine cmp_dp_wvfn( obj, iw )
      class(wavfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      integer(int32) :: ip
      obj%dp_ln_wvfn(:,iw) = 0.0_dp ; ip = 0_int32
      if (edl_opt) then
         select case( abs(wf_type_e) )
          case(1)
            call obj%sld_e%cmp_dl( iw, n_el, n_el_u, n_el_d, sld_c_e, forbs_mat(iw)%o(:,1:n_el), &
            & n_par_det_e, obj%dp_ln_wvfn(ip+1:ip+n_par_det_e,iw) )
          case(2)
            call obj%sgm_e%cmp_dl( iw, n_el, n_el_u, n_el_d, n_el_s, sgm_c_e, forbs_mat(iw)%o(:,1:n_el), &
            & n_par_det_e, obj%dp_ln_wvfn(ip+1:ip+n_par_det_e,iw) )
          case(3)
            call obj%tgm_e%cmp_dl( iw, n_el, n_el_u, n_el_d, tgm_c_e, forbs_mat(iw)%o(:,1:n_el), &
            & n_par_det_e, obj%dp_ln_wvfn(ip+1:ip+n_par_det_e,iw) )
          case(4)
            call obj%pff_e%cmp_dl( iw, n_el, n_el_u, n_el_d, pff_c_e, forbs_mat(iw)%o(:,1:n_el), &
            & n_par_det_e, obj%dp_ln_wvfn(ip+1:ip+n_par_det_e,iw) )
          case default
         end select
         ip = ip + n_par_det_e
      endif
      if ( pdl_opt ) then
         if (wf_type_p.ne.0) then
            select case ( abs(wf_type_p) )
             case(1)
               call obj%sld_p%cmp_dl( iw, n_po, n_po_u, n_po_d, sld_c_p, forbs_mat(iw)%o(:,n_el+1:n_fe), &
               & n_par_det_p, obj%dp_ln_wvfn(ip+1:ip+n_par_det_p,iw) )
             case(2)
               call obj%sgm_p%cmp_dl( iw, n_po, n_po_u, n_po_d, n_po_s, sgm_c_p, forbs_mat(iw)%o(:,n_el+1:n_fe), &
               & n_par_det_p, obj%dp_ln_wvfn(ip+1:ip+n_par_det_p,iw) )
             case(3)
               call obj%tgm_p%cmp_dl( iw, n_po, n_po_u, n_po_d, tgm_c_p, forbs_mat(iw)%o(:,n_el+1:n_fe), &
               & n_par_det_p, obj%dp_ln_wvfn(ip+1:ip+n_par_det_p,iw) )
             case(4)
               call obj%pff_p%cmp_dl( iw, n_po, n_po_u, n_po_d, pff_c_p, forbs_mat(iw)%o(:,n_el+1:n_fe), &
               & n_par_det_p, obj%dp_ln_wvfn(ip+1:ip+n_par_det_p,iw) )
             case default
            end select
            ip = ip + n_par_det_p
         endif
         if (wf_type_ep.ne.0) then
            select case ( abs(wf_type_ep) )
             case(11)
               call obj%epp%cmp_dl( iw, n_el, n_po, eporbs_mat(iw)%o, n_par_det_ep, obj%dp_ln_wvfn(ip+1:ip+n_par_det_ep,iw) )
            end select
            ip = ip + n_par_det_ep
         endif
      endif 
      if (aoc_opt) then
         if (wf_type_e.ne.0) then
            select case( abs(wf_type_e) )
             case(1)
               call obj%sld_e%cmp_da( iw, n_el, n_el_u, n_el_d, forbs_mat(iw)%da_o(:,1:n_el), &
               & obj%dp_ln_wvfn(ip+1:ip+n_par_forbs,iw) )
             case(2)
               call obj%sgm_e%cmp_da( iw, n_el, forbs_mat(iw)%da_o(:,1:n_el), &
               & obj%dp_ln_wvfn(ip+1:ip+n_par_forbs,iw) )
             case(3)
               call obj%tgm_e%cmp_da( iw, n_el, n_el_u, n_el_d, forbs_mat(iw)%da_o(:,1:n_el), &
               & obj%dp_ln_wvfn(ip+1:ip+n_par_forbs,iw) )
             case(4)
               call obj%pff_e%cmp_da( iw, n_el, forbs_mat(iw)%da_o(:,1:n_el), &
               & obj%dp_ln_wvfn(ip+1:ip+n_par_forbs,iw) )
             case default
            end select 
         endif
         if (wf_type_p.ne.0) then
            select case ( abs(wf_type_p) )
             case(1)
               call obj%sld_p%cmp_da( iw, n_po, n_po_u, n_po_d, forbs_mat(iw)%da_o(:,n_el+1:n_fe), &
               & obj%dp_ln_wvfn(ip+1:ip+n_par_forbs,iw) )
             case(2)
               call obj%sgm_p%cmp_da( iw, n_po, forbs_mat(iw)%da_o(:,n_el+1:n_fe), &
               & obj%dp_ln_wvfn(ip+1:ip+n_par_forbs,iw) )
             case(3)
               call obj%tgm_p%cmp_da( iw, n_po, n_po_u, n_po_d, forbs_mat(iw)%da_o(:,n_el+1:n_fe), &
               & obj%dp_ln_wvfn(ip+1:ip+n_par_forbs,iw) )
             case(4)
               call obj%pff_p%cmp_da( iw, n_po, forbs_mat(iw)%da_o(:,n_el+1:n_fe), &
               & obj%dp_ln_wvfn(ip+1:ip+n_par_forbs,iw) )
             case default
            end select 
         endif
         if (wf_type_e.ne.0.or.wf_type_p.ne.0) ip = ip + n_par_forbs
      endif
      if (wf_type_ep.ne.0) then
         select case ( abs(wf_type_ep) )
          case(11)
            if (n_eporbs.gt.0_int32.and.poc_opt) then
               call obj%epp%cmp_da( iw, n_el, n_po, epo_c, eporbs_mat(iw)%da_o, &
               & obj%dp_ln_wvfn(ip+1:ip+n_par_eporbs,iw) )
               ip = ip + n_par_eporbs
            endif
          case default
         end select
      endif
      if ( n_par_jst.gt.0_int32 ) &
      & call dp_jastrow_factors( iw, obj%dp_ln_wvfn(n_par_frm+1:n_par_frm_wvfn,iw) )
   end subroutine cmp_dp_wvfn
end module fermionic_wavefunction_c
