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
module local_energy_m
   use fortran_kinds_v,          only: int32, dp
   use quantum_monte_carlo_v,    only: n_wlk_max, qmc_mthd
   use wavefunction_optimization_v, only: opt_mthd
   use qdo_system_v,             only: n_qdo, n_pchrg, qdoham, qdos
   use molecular_system_v,       only: n_fe, n_el, n_po, n_at, m_e, m_p
   use pseudopotentials_v,       only: psd_prs
   use fermionic_potentials_m,   only: cmp_v_nn, cmp_v_fn, cmp_v_ee, cmp_v_pp, cmp_v_ep, &
   & upd_v_fn, upd_v_ee, upd_v_pp, upd_v_ep
   use qdo_potentials_m,         only: cmp_v_qq, cmp_v_pchpch, cmp_v_qpch, cmp_v_dq, cmp_v_dd,&
   & cmp_v_dpch, upd_v_dq, cmp_v_dd, upd_v_dpch
   use embedding_potential_m,    only: cmp_v_qn, cmp_v_npch, cmp_v_dn, cmp_v_fqd, cmp_v_fpch, &
   & upd_v_fqd, upd_v_fpch, upd_v_dn
   use qdo_wavefunction_m,       only: wvfnq
   use fermionic_wavefunction_m, only: wvfn
   use jstqe_mod,                only: jstqe_fct
   use jstqepar_var,             only: jstqe_prs
   use external_field_m,         only: ext_E_field, cmp_extfld_frm, cmp_extfld_drd
   use pseudopotentials_m,       only: eval_pseudo_ene
   use kinetic_energy_m,         only: cmp_kinene, upd_kinene, upd_kinene_emb
   implicit none
   real(dp), public, save, allocatable, dimension(:)   :: e_l
   real(dp), public, save, allocatable, dimension(:)   :: k_l
   real(dp), public, save, allocatable, dimension(:)   :: v_l
   real(dp), public, save, allocatable, dimension(:)   :: v_psd
   real(dp), public, save, allocatable, dimension(:,:) :: e_l_pp
   real(dp), public, save, allocatable, dimension(:,:) :: v_psd_pp
   real(dp), public, save, allocatable, dimension(:,:) :: v_fn
   real(dp), public, save, allocatable, dimension(:,:) :: v_ee
   real(dp), public, save, allocatable, dimension(:,:)   :: v_pp
   real(dp), public, save, allocatable, dimension(:,:,:) :: v_ep
   real(dp), public, save, allocatable, dimension(:,:) :: v_dq
   real(dp), public, save, allocatable, dimension(:,:) :: v_dd
   real(dp), public, save, allocatable, dimension(:,:) :: v_dn
   real(dp), public, save, allocatable, dimension(:,:) :: v_fq
   real(dp), public, save, allocatable, dimension(:,:,:) :: v_df
   real(dp), public, save, allocatable, dimension(:,:) :: v_fpch
   real(dp), public, save, allocatable, dimension(:,:) :: v_dpch
   real(dp), public, save :: v_ext
   real(dp), public, save :: v_nn
   real(dp), public, save :: v_qq
   real(dp), public, save :: v_qn
   real(dp), public, save :: v_pchpch
   real(dp), public, save :: v_npch
   real(dp), public, save :: v_qpch
   public :: ini_locene, cmp_locene
contains
   subroutine ini_locene()
      call cmp_v_nn( v_nn )
      if (qdoham.eq.'cou') then
         call cmp_v_qq( v_qq )
      else
         v_qq = 0.0_dp
      endif
      call cmp_v_qn( v_qn )
      call cmp_v_pchpch( v_pchpch )
      call cmp_v_qpch( v_qpch )
      call cmp_v_npch( v_npch )
      v_ext = v_nn + v_qq + v_qn + v_pchpch + v_qpch + v_npch
      allocate( e_l(1:n_wlk_max) ) ; e_l = 0.0_dp
      allocate( k_l(1:n_wlk_max) ) ; k_l = 0.0_dp
      allocate( v_l(1:n_wlk_max) ) ; v_l = 0.0_dp
      if ( n_fe+n_qdo.ge.1_int32 ) then
         allocate( e_l_pp(1:n_fe+n_qdo,1:n_wlk_max) ) ; e_l_pp = 0.0_dp
         if (psd_prs) then
            allocate( v_psd(1:n_wlk_max) ) ; v_psd = 0.0_dp
            allocate( v_psd_pp(1:n_fe+n_qdo,1:n_wlk_max) ) ; v_psd_pp = 0.0_dp
         endif
      endif
      if ( n_at.ge.1_int32 ) then
         allocate( v_fn(1:n_fe,1:n_wlk_max) ) ; v_fn = 0.0_dp
      endif
      if ( n_el.ge.2_int32 ) then
         allocate( v_ee(1:n_el*(n_el-1)/2,1:n_wlk_max) ) ; v_ee = 0.0_dp
      endif
      if ( n_po.ge.2_int32 ) then
         allocate( v_pp(1:n_po*(n_po-1)/2,1:n_wlk_max) ) ; v_pp = 0.0_dp
      endif
      if ( n_el.ge.1_int32 .and. n_po.ge.1_int32 ) then
         allocate( v_ep(1:n_po,1:n_el,1:n_wlk_max) ) ; v_ep = 0.0_dp
      endif
      if ( n_qdo.gt.0_int32 ) then
         allocate( v_dq(1:n_qdo,1:n_wlk_max) )   ; v_dq = 0.0_dp
         if ( n_qdo.ge.2_int32 ) then
            allocate( v_dd(1:n_qdo*(n_qdo-1)/2,1:n_wlk_max) ) ; v_dd = 0.0_dp
         endif
         if ( n_at.gt.0_int32 ) then
            allocate( v_dn(1:n_qdo,1:n_wlk_max) ) ; v_dn = 0.0_dp
         endif
         if ( n_fe.gt.0_int32 ) then
            allocate( v_fq(1:n_fe,1:n_wlk_max) )              ; v_fq = 0.0_dp
            allocate( v_df(1:n_qdo,1:n_fe,1:n_wlk_max) )      ; v_df = 0.0_dp
         endif
         if ( n_pchrg.gt.0_int32 ) then
            allocate( v_dpch(1:n_qdo,1:n_wlk_max) )           ; v_dpch = 0.0_dp
         endif
         if( (n_pchrg.gt.0_int32).and.(n_at.gt.0_int32) ) then
            allocate( v_fpch(1:n_fe,1:n_wlk_max) )            ; v_fpch = 0.0_dp
         endif
      endif
   end subroutine ini_locene
   subroutine cmp_locene( iw, K_cmp )
      integer(int32), intent(in) :: iw
      logical,        intent(in) :: K_cmp
      real(dp)                   :: DD2_tmp(1:4)
      integer(int32) :: i1, i2, i3
      if ( K_cmp ) call cmp_kinene( iw )
      if ( n_fe.gt.0_int32 ) then
         do i1 = 1_int32, n_fe
            if ( jstqe_prs ) then
               DD2_tmp(1:4) =  wvfn%DD2_ln_wvfn(4*i1-3:4*i1,iw) + jstqe_fct(iw)%DD2e_jstqe(4*i1-3:4*i1)
            else
               DD2_tmp(1:4) =  wvfn%DD2_ln_wvfn(4*i1-3:4*i1,iw)
            endif
            if ( i1.le.n_el ) then
               e_l_pp(i1,iw) = - 0.5_dp/m_e*( DD2_tmp(4) + sum(DD2_tmp(1:3)**2) )
            else
               e_l_pp(i1,iw) = - 0.5_dp/m_p*( DD2_tmp(4) + sum(DD2_tmp(1:3)**2) )
            endif
         enddo ! i1 (n_fe)
      endif
      if ( n_qdo.gt.0_int32 ) then
         do i1 = 1_int32, n_qdo
            e_l_pp(i1+n_fe,iw) = - 0.5_dp/qdos(i1)%qdo_m &
               & *(wvfnq(iw)%DD2_ln_wvfnq(4*i1)+sum(wvfnq(iw)%DD2_ln_wvfnq(4*i1-3:4*i1-1)**2) )
         enddo ! i1 (n_qdo)
      endif
      k_l(iw) = sum(e_l_pp(:,iw))
      v_l(iw) = 0.0_dp
      if ( n_fe.gt.0_int32 ) then
         if ( n_at.gt.0_int32  ) then
            call cmp_v_fn( iw, v_fn(:,iw) )
            v_l(iw) = v_l(iw) + sum( v_fn(:,iw) )
            e_l_pp(1:n_fe,iw) = e_l_pp(1:n_fe,iw) + v_fn(1:n_fe,iw)
            if ( psd_prs ) then
               call eval_pseudo_ene( iw, v_psd_pp(1:n_fe,iw) )
               v_psd(iw) = sum(v_psd_pp(1:n_fe,iw))
               v_l(iw) = v_l(iw) + v_psd(iw)
               e_l_pp(1:n_fe,iw) = e_l_pp(1:n_fe,iw) + v_psd_pp(1:n_fe,iw)
            endif
         endif
         if ( n_el.ge.2_int32 ) then
            call cmp_v_ee( iw, v_ee(:,iw) )
            i3 = 0_int32
            do i1 = 1_int32, n_el ; do i2 = i1 + 1_int32, n_el
                  i3 = i3 + 1_int32
                  e_l_pp(i1,iw) = e_l_pp(i1,iw) + 0.5 * v_ee(i3,iw)
                  e_l_pp(i2,iw) = e_l_pp(i2,iw) + 0.5 * v_ee(i3,iw)
               enddo ; enddo
            v_l(iw) = v_l(iw) + sum(v_ee(:,iw))
         endif
         if ( n_po.ge.2_int32 ) then
            call cmp_v_pp( iw, v_pp(:,iw) )
            i3 = 0_int32
            do i1 = 1_int32, n_po ; do i2 = i1 + 1_int32, n_po
                  i3 = i3 + 1_int32
                  e_l_pp(i1+n_el,iw) = e_l_pp(i1+n_el,iw) + 0.5 * v_pp(i3,iw)
                  e_l_pp(i2+n_el,iw) = e_l_pp(i2+n_el,iw) + 0.5 * v_pp(i3,iw)
               enddo ; enddo
            v_l(iw) = v_l(iw) + sum(v_pp(:,iw))
         endif
         if ( n_el.gt.0_int32 .and. n_po.gt.0_int32 ) then
            call cmp_v_ep( iw, v_ep(:,:,iw) )
            do i1 = 1_int32, n_el ; do i2 = 1_int32, n_po
                  e_l_pp(i1,iw) = e_l_pp(i1,iw) + 0.5 * v_ep(i2,i1,iw)
                  e_l_pp(i2+n_el,iw) = e_l_pp(i2+n_el,iw) + 0.5 * v_ep(i2,i1,iw)
               enddo ; enddo
            v_l(iw) = v_l(iw) + sum(v_ep(:,:,iw))
         endif
         if ( (sum(abs(ext_E_field)).gt.0.0_dp)) call cmp_extfld_frm( iw, e_l_pp(1:n_fe,iw), v_l(iw) )
      endif 
      if ( n_qdo.gt.0_int32 ) then
         call cmp_v_dq( iw, v_dq(:,iw) )
         v_l(iw) = v_l(iw) + sum(v_dq(:,iw))
         e_l_pp(n_fe+1:n_fe+n_qdo,iw) = e_l_pp(n_fe+1:n_fe+n_qdo,iw) + v_dq(1:n_qdo,iw)
         if ( n_qdo.ge.2_int32 ) then
            call cmp_v_dd( iw, v_dd(:,iw) )
            i3 = 0_int32
            do i1 = 1_int32, n_qdo ; do i2 = i1 + 1_int32, n_qdo
                  i3 = i3 + 1_int32
                  e_l_pp(n_fe+i1,iw) = e_l_pp(n_fe+i1,iw) + 0.5 * v_dd(i3,iw)
                  e_l_pp(n_fe+i2,iw) = e_l_pp(n_fe+i2,iw) + 0.5 * v_dd(i3,iw)
               enddo ; enddo
            v_l(iw) = v_l(iw) + sum(v_dd(:,iw))
         endif
         if ( n_at.gt.0_int32 ) then
            call cmp_v_dn( iw, v_dn(:,iw) )
            v_l(iw) = v_l(iw) + sum(v_dn(:,iw))
            e_l_pp(n_fe+1:n_fe+n_qdo,iw) = e_l_pp(n_fe+1:n_fe+n_qdo,iw) + v_dn(1:n_qdo,iw)
         endif
         if ( n_fe.gt.0_int32 ) then
            call cmp_v_fqd( iw, v_fq(:,iw), v_df(:,:,iw) )
            v_l(iw) = v_l(iw) + sum(v_fq(:,iw)) + sum(v_df(:,:,iw))
            e_l_pp(1:n_fe,iw) = e_l_pp(1:n_fe,iw) + v_fq(1:n_fe,iw)
            do i1 = 1_int32, n_fe ; do i2 = 1_int32, n_qdo
                  e_l_pp(i1,iw)      = e_l_pp(i1,iw)      + 0.5_dp * v_df(i2,i1,iw)
                  e_l_pp(n_fe+i2,iw) = e_l_pp(n_fe+i2,iw) + 0.5_dp * v_df(i2,i1,iw)
               enddo ; enddo
         endif
         if ( n_pchrg.gt.0_int32 ) then
            call cmp_v_dpch( iw, v_dpch(:,iw) )
            v_l(iw) = v_l(iw) + sum(v_dpch(:,iw))
            e_l_pp(n_fe+1:n_fe+n_qdo,iw) = e_l_pp(n_fe+1:n_fe+n_qdo,iw) + v_dpch(1:n_qdo,iw)
         endif
         if ( (n_pchrg.gt.0_int32).and.(n_fe.gt.0_int32) ) then
            call cmp_v_fpch( iw, v_fpch(:,iw) )
            v_l(iw) = v_l(iw) + sum(v_fpch(:,iw))
            e_l_pp(1:n_fe,iw) = e_l_pp(1:n_fe,iw) + v_fpch(1:n_fe,iw)
         endif
         if (sum(abs(ext_E_field)).gt.0.0_dp ) call cmp_extfld_drd( iw, e_l_pp(n_fe+1:n_fe+n_qdo,iw), v_l(iw) )
      endif 
      e_l(iw) = k_l(iw) + v_l(iw) + v_ext
   end subroutine cmp_locene
end module local_energy_m
