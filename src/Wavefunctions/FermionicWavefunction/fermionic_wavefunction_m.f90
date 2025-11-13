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
module fermionic_wavefunction_m
   use fortran_kinds_v, only: dp, int32, stdout
   use openmp_mpi_m,    only: mpi_rank
   use write_lines_m,   only: write_separator_line, write_simple_line, write_empty_line, &
   & write_variable_line
   use molecular_system_v, only: n_fe, n_el, n_el_u, n_el_d, n_po, n_po_u, n_po_d, n_at
   use jastrow_factors_v
   use fermionic_wavefunction_v
   use quantum_monte_carlo_v, only: restart, n_wlk_max
   use fermionic_wavefunction_c, only: wavfun_t
   use fermionic_orbitals_m, only: ini_forbs_mat
   use feporb_mod, only: ini_eporbs_mat
   use jeporb_mod, only: ini_jepos_mat
   use jastrow_orbitals_m, only: ini_jorbs_mat
   use slater_determinant_params_m, only: ini_sldcfs, sld_c_e, sld_c_p
   use singlet_geminal_params_m, only: ini_sgmcfs, sgm_c_e, sgm_c_p
   use triplet_geminal_params_m, only: ini_tgmcfs, tgm_c_e, tgm_c_p
   use pfaffian_params_m, only: ini_pffcfs, pff_c_e, pff_c_p
   use eppcfs_mod, only: ini_eppcfs
   use jastrow_factors_m, only: init_jastrow_factors
   use jastrow_factor_params_m, only: read_jastrow_factor_params
   use wavefunction_optimization_v, only: n_par_opt, n_par_opt_s
   implicit none
   type(wavfun_t), public, save, target :: wvfn
   public  :: init_fermionic_wavefunction_v, upd_fermionic_wavefunction_v, save_fermionic_wavefunction_v
contains
   subroutine init_fermionic_wavefunction_v()
      character(3), external :: wvfn_name
      n_par_frm_wvfn   = 0_int32 ; n_par_frm_wvfn_s   = 0_int32
      n_par_frm    = 0_int32 ; n_par_frm_s    = 0_int32
      n_par_det_e  = 0_int32 ; n_par_det_e_s  = 0_int32
      n_par_det_p  = 0_int32 ; n_par_det_p_s  = 0_int32
      n_par_det_ep = 0_int32 ; n_par_det_ep_s = 0_int32
      n_par_det    = 0_int32 ; n_par_det_s    = 0_int32
      if ( n_at.eq.0_int32 ) then
         wf_type_e = 0
         wf_type_p = 0
      endif
      if (n_el.eq.0_int32) then
         wf_type_e = 0
         edl_opt = .false.
      endif
      if (n_po.eq.0_int32) then
         wf_type_p  = 0
         wf_type_ep = 0
         pdl_opt = .false.
      endif
      if ( abs(wf_type_e).eq.4_int32 ) then
         if ( n_el_u.lt.2_int32 .and. n_el_d.lt.2_int32 ) wf_type_e = sign(1_int32, wf_type_e) * 2
      endif
      if ( wf_type_e.eq.-2_int32 ) then
         if ( n_el_u.eq.1_int32 .and. n_el_d.eq.1_int32 ) wf_type_e = 2
      endif
      if ( abs(wf_type_e).eq.4_int32 .or. abs(wf_type_e).eq.3_int32 ) then
         if ( n_el_u.ne.n_el_d) wf_type_e = -abs(wf_type_e)
      endif
      if ( abs(wf_type_e).eq.2_int32 .or. abs(wf_type_e).eq.4_int32 ) then
         if ( n_el_u.eq.0_int32 .or. n_el_d.eq.0_int32 ) wf_type_e = -3
      endif
      if ( abs(wf_type_e).eq.3_int32 ) then
         if ( n_el_u.eq.1_int32 .and. n_el_d.eq.1_int32 ) wf_type_e = sign(1_int32, wf_type_e)
      endif
      if ( abs(wf_type_p).eq.4_int32 ) then
         if ( n_po_u.lt.2_int32 .and. n_po_d.lt.2_int32 ) wf_type_p = sign(1_int32, wf_type_p) * 2
      endif
      if ( wf_type_p.eq.-2_int32 ) then
         if ( n_po_u.eq.1_int32 .and. n_po_d.eq.1_int32 ) wf_type_p = 2
      endif
      if ( abs(wf_type_p).eq.4_int32 .or. abs(wf_type_p).eq.3_int32 ) then
         if ( n_po_u.ne.n_po_d) wf_type_p = -abs(wf_type_p)
      endif
      if ( abs(wf_type_p).eq.2_int32 .or. abs(wf_type_p).eq.4_int32 ) then
         if ( n_po_u.eq.0_int32 .or. n_po_d.eq.0_int32 ) wf_type_p = -3
      endif
      if ( abs(wf_type_p).eq.3_int32 ) then
         if ( n_po_u.eq.1_int32 .and. n_po_d.eq.1_int32 ) wf_type_p = sign(1_int32, wf_type_p)
      endif
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","INITIALIZING WAVE FUNCTION")
      call write_simple_line(stdout,0,mpi_rank,2,"c","________________________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","Wave function types")
      call write_empty_line(stdout,0,mpi_rank)
      if (wf_type_e.ne.0) then
         call write_variable_line(stdout,0,mpi_rank,2,"Wave function type for electrons", wvfn_name( wf_type_e ),var_name="wf_type_e")
      endif
      if (wf_type_p.ne.0) then
         call write_variable_line(stdout,0,mpi_rank,2,"Wave function type for positrons", wvfn_name( wf_type_p ),var_name="wf_type_p")
      endif
      if ( wf_type_ep.ne.0 ) then
         call write_variable_line(stdout,0,mpi_rank,2,"Wave function for ele.-pos. coupling", wvfn_name( wf_type_ep ),var_name="wf_type_ep")
      endif
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"l","Initializing atomic orbitals.")
      call write_empty_line(stdout,0,mpi_rank)
      call ini_forbs_mat()
      call ini_jorbs_mat()
      call ini_eporbs_mat()
      call ini_jepos_mat()
      if (restart) then
         call read_forbs_par()
         call read_eporbs_par()
      endif
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"l","Initializing determinantal part of the wave function.")
      call write_empty_line(stdout,0,mpi_rank)
      if ( wf_type_e.lt.0_int32 ) then
         spin = 'U'
      else
         spin = 'R'
      endif
      select case ( abs(wf_type_e) )
       case(1)
         call ini_sldcfs( sld_c_e , n_el_u, n_el_d, edl_opt, n_par_det_e, n_par_det_e_s, head_fle='slde' )
       case(2)
         call ini_sgmcfs( sgm_c_e, n_el_u, n_el_d, edl_opt, n_par_det_e, n_par_det_e_s, head_fle='sgme' )
       case(3)
         call ini_tgmcfs( tgm_c_e, n_el_u, n_el_d, edl_opt, n_par_det_e, n_par_det_e_s, head_fle='tgme' )
       case(4)
         call ini_pffcfs( pff_c_e, n_el, n_el_d, edl_opt, n_par_det_e, n_par_det_e_s, fle_head='pffe' )
       case default
      end select 
      select case ( abs(wf_type_p) )
       case(1)
         call ini_sldcfs( sld_c_p, n_po_u, n_po_d, pdl_opt, n_par_det_p, n_par_det_p_s, head_fle='sldp' )
       case(2)
         call ini_sgmcfs( sgm_c_p, n_po_u, n_po_d, pdl_opt, n_par_det_p, n_par_det_p_s, head_fle='sgmp' )
       case(3)
         call ini_tgmcfs( tgm_c_p, n_po_u, n_po_d, pdl_opt, n_par_det_p, n_par_det_p_s, head_fle='tgmp' )
       case(4)
         call ini_pffcfs( pff_c_p, n_po, n_po_d, pdl_opt, n_par_det_p, n_par_det_p_s, fle_head='pffp' )
       case default
      end select 
      select case ( abs(wf_type_ep) )
       case(11)
         call ini_eppcfs( pdl_opt, n_par_det_ep, n_par_det_ep_s, head_fle='epog' )
       case default
      end select
      n_par_det   = n_par_det + n_par_det_e + n_par_det_p + n_par_det_ep
      n_par_det_s = n_par_det_s + n_par_det_e_s + n_par_det_p_s + n_par_det_ep_s
      n_par_frm   = n_par_frm   + n_par_det
      n_par_frm_s = n_par_frm_s + n_par_det_s
      call init_jastrow_factors()
      n_par_frm_wvfn   = n_par_frm    + n_par_jst
      n_par_frm_wvfn_s = n_par_frm_s  + n_par_jst_s
      n_par_opt   = n_par_opt + n_par_frm_wvfn
      n_par_opt_s = n_par_opt_s + n_par_frm_wvfn_s
      call wvfn%ini( )
   end subroutine init_fermionic_wavefunction_v
   subroutine upd_fermionic_wavefunction_v( vec_wvfn_var )
      use molecular_system_v, only: n_el, n_el_d, n_el_u, n_el_s, n_po, n_po_d, n_po_s, n_po_u
      use fermionic_orbitals_m, only: forbs_s, n_par_forbs_s
      use jastrow_orbitals_m, only: jc_opt, jorbs_s, n_par_jorbs_s
      use feporb_mod, only: eporbs_s, n_par_eporbs_s
      use jeporb_mod, only: jpoc_opt, jepos_s, n_par_jepos_s
      use jastrow_factors_v, only: n_par_jst, n_par_jstd_s, n_par_jstc, n_par_jepd_s
      use fermionic_wavefunction_v, only: wf_type_e, wf_type_p, wf_type_ep
      use slater_determinant_params_m, only: upd_sldcfs, sld_c_e, sld_c_p
      use singlet_geminal_params_m, only: upd_sgmcfs, sgm_c_e, sgm_c_p
      use triplet_geminal_params_m, only: upd_tgmcfs, tgm_c_e, tgm_c_p
      use pfaffian_params_m, only: upd_pffcfs, pff_c_e, pff_c_p
      use eppcfs_mod, only: upd_eppcfs
      use jastrow_factor_params_m, only: updt_jastrow_factor_params
      real(dp), dimension(n_par_frm_wvfn_s), intent(in) :: vec_wvfn_var
      integer(int32) :: ip
      ip = 1_int32
      if ( edl_opt ) then
         select case ( abs(wf_type_e) )
          case(1)
            call upd_sldcfs( sld_c_e, n_el_u ,n_el_d, n_par_det_e_s, vec_wvfn_var(1:n_par_det_e_s) )
          case(2)
            call upd_sgmcfs( sgm_c_e, n_el_s, n_par_det_e_s, vec_wvfn_var(1:n_par_det_e_s) )
          case(3)
            call upd_tgmcfs( tgm_c_e, n_par_det_e_s, vec_wvfn_var(1:n_par_det_e_s) )
          case(4)
            call upd_pffcfs( pff_c_e, n_el, n_el_d, n_par_det_e_s, vec_wvfn_var(1:n_par_det_e_s) )
          case default
         end select 
         ip = ip + n_par_det_e_s
      endif
      if ( pdl_opt ) then
         select case ( abs(wf_type_p) )
          case(1)
            call upd_sldcfs( sld_c_p, n_po_u, n_po_d, n_par_det_p_s, vec_wvfn_var(ip:ip+n_par_det_p_s-1) )
          case(2)
            call upd_sgmcfs( sgm_c_p, n_po_s, n_par_det_p_s, vec_wvfn_var(ip:ip+n_par_det_p_s-1) )
          case(3)
            call upd_tgmcfs( tgm_c_p, n_par_det_p_s, vec_wvfn_var(ip:ip+n_par_det_p_s-1) )
          case(4)
            call upd_pffcfs( pff_c_p, n_po, n_po_d, n_par_det_p_s, vec_wvfn_var(ip:ip+n_par_det_p_s-1) )
          case default
         end select 
         ip = ip + n_par_det_p_s
      endif
      if ( wf_type_ep.ne.0 ) then
         select case ( abs(wf_type_ep) )
          case(11)
            if ( pdl_opt ) then
               call upd_eppcfs( n_par_det_ep_s, vec_wvfn_var(ip:ip+n_par_det_ep_s-1) )
               ip = ip + n_par_det_ep_s
            endif
          case default
         end select 
      endif
      if (n_par_forbs_s.gt.0_int32) then
         call forbs_s%upd( vec_wvfn_var(ip:ip+n_par_forbs_s-1), 0.0001_dp )
         ip = ip + n_par_forbs_s
      endif
      if (n_par_eporbs_s.gt.0_int32) then
         call eporbs_s%upd( vec_wvfn_var(ip:ip+n_par_eporbs_s-1), 0.0001_dp)
         ip = ip + n_par_eporbs_s
      endif
      if ( n_par_jst.gt.0_int32 ) then
         call updt_jastrow_factor_params( vec_wvfn_var(ip:ip+n_par_jstc+n_par_jstd_s+n_par_jepd_s+n_par_jepa_s-1) )
         ip = ip + n_par_jstc + n_par_jstd_s + n_par_jepd_s + n_par_jepa_s
         if ( jc_opt ) then
            call jorbs_s%upd( vec_wvfn_var(ip:ip+n_par_jorbs_s-1), 0.0001_dp )
            ip = ip + n_par_jorbs_s
         endif
         if ( jpoc_opt ) then
            call jepos_s%upd( vec_wvfn_var(ip:ip+n_par_jepos_s-1), 0.0001_dp )
         endif
      endif
   end subroutine upd_fermionic_wavefunction_v
   subroutine save_fermionic_wavefunction_v()
      use molecular_system_v, only: n_po, n_el, n_el_u, n_el_d, n_el_s, n_po_d, n_po_u, n_po_s
      use fermionic_wavefunction_v, only: wf_type_e, wf_type_p
      use slater_determinant_params_m, only: save_sldcfs, save_sldcfs_st, sld_c_e, sld_c_p
      use singlet_geminal_params_m, only: save_sgmcfs, save_sgmcfs_st, sgm_c_e, sgm_c_p
      use triplet_geminal_params_m, only: save_tgmcfs, save_tgmcfs_st, tgm_c_e, tgm_c_p
      use pfaffian_params_m, only: save_pffcfs, save_pffcfs_st, pff_c_e, pff_c_p
      use eppcfs_mod, only: save_eppcfs, save_eppcfs_st
      use jastrow_factor_params_m, only: save_jastrow_factor_params
      call system('mkdir -p wvfn.save')
      if (wf_type_e.ne.0) then
         select case ( abs(wf_type_e) )
          case(1)
            call save_sldcfs( sld_c_e, n_el_u, n_el_d, 'slde' )
            call save_sldcfs_st( sld_c_e, n_el_d, 'slde' )
          case(2)
            call save_sgmcfs( sgm_c_e, n_el_s, 'sgme' )
            call save_sgmcfs_st( sgm_c_e, n_el_s, 'sgme' )
          case(3)
            call save_tgmcfs( tgm_c_e, 'tgme' )
            call save_tgmcfs_st( tgm_c_e, 'tgme' )
          case(4)
            call save_pffcfs( pff_c_e, n_el, n_el_d, 'pffe' )
            call save_pffcfs_st( pff_c_e, 'pffe' )
          case default
         end select
      endif
      if ( wf_type_p.ne.0_int32 ) then
         select case ( abs(wf_type_p) )
          case(1)
            call save_sldcfs( sld_c_p, n_po_u, n_po_d, 'sldp' )
            call save_sldcfs_st( sld_c_p, n_po_d, 'sldp' )
          case(2)
            call save_sgmcfs( sgm_c_p, n_po_s, 'sgmp' )
            call save_sgmcfs_st( sgm_c_p, n_po_s, 'sgmp' )
          case(3)
            call save_tgmcfs( tgm_c_p, 'tgmp' )
            call save_tgmcfs_st( tgm_c_p, 'tgmp' )
          case(4)
            call save_pffcfs( pff_c_p, n_po, n_po_d, 'pffp' )
            call save_pffcfs_st( pff_c_p, 'pffp' )
          case default
         end select 
      endif
      if ( wf_type_ep.ne.0_int32 ) then
         select case ( abs(wf_type_ep) )
          case(11)
            call save_eppcfs( 'epog' )
            call save_eppcfs_st( 'epog' )
          case default
         end select 
      endif
      call save_jastrow_factor_params()
      call save_forbs_par()
      call save_eporbs_par( )
   end subroutine save_fermionic_wavefunction_v
   subroutine sym_wvfn_par( O, O_s )
      use molecular_system_v, only: n_el, n_el_d, n_el_s, n_po, n_po_d, n_po_s
      use fermionic_wavefunction_v, only: edl_opt, pdl_opt, wf_type_e, wf_type_p
      use slater_determinant_params_m, only: sym_sldcfs, sld_c_e, sld_c_p
      use singlet_geminal_params_m, only: sym_sgmcfs, sgm_c_e, sgm_c_p
      use triplet_geminal_params_m, only: sym_tgmcfs, tgm_c_e, tgm_c_p
      use pfaffian_params_m, only: sym_pffcfs, pff_c_e, pff_c_p
      use eppcfs_mod, only: sym_eppcfs
      use fermionic_orbitals_m, only: forbs_s, n_par_forbs_s, n_par_forbs
      use feporb_mod, only: eporbs_s, n_par_eporbs_s, n_par_eporbs
      use jeporb_mod, only: jpoc_opt, jepos_s, n_par_jepos_s, n_par_jepos
      use jastrow_orbitals_m, only: jc_opt, jorbs_s, n_par_jorbs_s, n_par_jorbs
      use jastrow_factor_params_m, only: applysym_jastrow_factor_params
      real(dp), dimension(n_par_frm_wvfn),   intent(in) :: O
      real(dp), dimension(n_par_frm_wvfn_s), intent(inout) :: O_s
      integer(int32)                                    :: ip, ip_s
      ip = 0 ; ip_s = 0
      if ( edl_opt ) then
         select case( abs(wf_type_e) )
          case(1)
            call sym_sldcfs( sld_c_e, n_el_d, n_par_det_e, n_par_det_e_s, O(ip+1:ip+n_par_det_e), O_s(ip_s+1:ip_s+n_par_det_e_s) )
          case(2)
            call sym_sgmcfs( sgm_c_e, n_el_s, n_par_det_e, n_par_det_e_s, O(ip+1:ip+n_par_det_e), O_s(ip_s+1:ip_s+n_par_det_e_s) )
          case(3)
            call sym_tgmcfs( tgm_c_e, n_par_det_e, n_par_det_e_s, O(ip+1:ip+n_par_det_e), O_s(ip_s+1:ip_s+n_par_det_e_s) )
          case(4)
            call sym_pffcfs( pff_c_e, n_el, n_el_d, n_par_det_e, n_par_det_e_s, O(ip+1:ip+n_par_det_e), O_s(ip_s+1:ip_s+n_par_det_e_s) )
          case default
         end select 
         ip   = ip + n_par_det_e
         ip_s = ip_s + n_par_det_e_s
      endif
      if (pdl_opt) then
         if (wf_type_p.ne.0) then
            select case( abs(wf_type_p) )
             case(1)
               call sym_sldcfs( sld_c_p, n_po_d, n_par_det_p, n_par_det_p_s, O(ip+1:ip+n_par_det_p), O_s(ip_s+1:ip_s+n_par_det_p_s) )
             case(2)
               call sym_sgmcfs( sgm_c_p, n_po_s, n_par_det_p, n_par_det_p_s,  O(ip+1:ip+n_par_det_p), O_s(ip_s+1:ip_s+n_par_det_p_s) )
             case(3)
               call sym_tgmcfs( tgm_c_p, n_par_det_p, n_par_det_p_s, O(ip+1:ip+n_par_det_p), O_s(ip_s+1:ip_s+n_par_det_p_s) )
             case(4)
               call sym_pffcfs( pff_c_p, n_po, n_po_d, n_par_det_p, n_par_det_p_s, O(ip+1:ip+n_par_det_p), O_s(ip_s+1:ip_s+n_par_det_p_s) )
             case default
            end select 
            ip   = ip + n_par_det_p
            ip_s = ip_s + n_par_det_p_s
         endif
      endif
      if (wf_type_ep.ne.0) then
         select case( abs(wf_type_ep) )
          case(11)
            if (pdl_opt) call sym_eppcfs( n_par_det_ep, n_par_det_ep_s, O(ip+1:ip+n_par_det_ep), O_s(ip_s+1:ip_s+n_par_det_ep_s) )
          case default
         end select 
         ip   = ip + n_par_det_ep
         ip_s = ip_s + n_par_det_ep_s
      endif
      if (n_par_forbs_s.gt.0_int32) then
         call forbs_s%sym(O(ip+1:ip+n_par_forbs),O_s(ip_s+1:ip_s+n_par_forbs_s) )
         ip   = ip + n_par_forbs
         ip_s = ip_s + n_par_forbs_s
      endif
      if (n_par_eporbs_s.gt.0_int32) then
         call eporbs_s%sym(O(ip+1:ip+n_par_eporbs),O_s(ip_s+1:ip_s+n_par_eporbs_s) )
         ip   = ip + n_par_eporbs
         ip_s = ip_s + n_par_eporbs_s
      endif
      if ( n_par_jst.gt.0_int32 ) then
         call applysym_jastrow_factor_params( O(ip+1:ip+n_par_jstc+n_par_jstd+n_par_jepd+n_par_jepa), &
         & O_s(ip_s+1:ip_s+n_par_jstc+n_par_jstd_s+n_par_jepd_s+n_par_jepa_s) )
         ip   = ip + n_par_jstc+n_par_jstd+n_par_jepd+ n_par_jepa
         ip_s = ip_s + n_par_jstc+n_par_jstd_s+n_par_jepd_s+ n_par_jepa_s
         if( jc_opt ) then
            call jorbs_s%sym( O(ip+1:ip+n_par_jorbs), O_s(ip_s+1:ip_s+n_par_jorbs_s) )
            ip   = ip + n_par_jorbs
            ip_s = ip_s + n_par_jorbs_s
         endif
         if( jpoc_opt ) then
            call jepos_s%sym( O(ip+1:ip+n_par_jepos), O_s(ip_s+1:ip_s+n_par_jepos_s))
         endif
      endif
   end subroutine sym_wvfn_par
   subroutine ini_frmwvf_red( wf_type_e_new, wf_type_p_new )
      use openmp_mpi_m, only: mpi_rank
      use quantum_monte_carlo_v, only: restart
      use molecular_system_v, only: n_fe, n_el, n_el_u, n_el_d, n_po, n_po_u, n_po_d, n_at
      use slater_determinant_params_m, only: ini_sldcfs, sld_c_e, sld_c_p
      use singlet_geminal_params_m, only: ini_sgmcfs, sgm_c_e, sgm_c_p
      use triplet_geminal_params_m, only: ini_tgmcfs, tgm_c_e, tgm_c_p
      use pfaffian_params_m, only: ini_pffcfs, pff_c_e, pff_c_p
      use eppcfs_mod, only: ini_eppcfs
      integer(int32), intent(inout) :: wf_type_e_new, wf_type_p_new
      character(3), external :: wvfn_name
      restart = .false.
      if (n_el.eq.0_int32) then
         wf_type_e_new = 0
      endif
      if (n_po.eq.0_int32) then
         wf_type_p_new = 0
      endif
      if (n_fe.gt.0_int32) then
         if ( n_el.eq.1_int32 ) then
            wf_type_e_new = -1
         endif
         if ( n_el_d.eq.0_int32 ) then
            if ( wf_type_e_new.eq.4 ) wf_type_e_new = -3
         endif
         if ( n_po.gt.0_int32 .and. wf_type_p_new.ne.0_int32 ) then
            if ( n_po.eq.1_int32 ) then
               wf_type_p_new = -1
            endif
            if ( n_po_d.eq.0_int32 ) then
               if ( wf_type_p_new.eq.4 ) wf_type_p_new = -3
            endif
         else 
            wf_type_p_new = 0
         endif 
      endif
      if ( n_at.eq.0_int32) then
         wf_type_e_new = 0
         wf_type_p_new = 0
      endif
      if( mpi_rank.eq.0_int32 ) then
         call write_empty_line(stdout,0,mpi_rank)
         write(*,'(3X,"================================================================")')
         write(*,'(3X,"                 WAVE FUNCTION INITIALIZATION                   ")')
         write(*,'(3X,"               ________________________________                 ")')
         write(*,'(3X,"                      Wave function types                       ")')
         call write_empty_line(stdout,0,mpi_rank)
         if (wf_type_e_new.ne.0) then
            write(*,'(3X,"Wave function type for electrons          (wf_type_e) : ",5X,A3)') wvfn_name( wf_type_e_new )
         else
            write(*,'(3X,"No fermionic wave function for electrons  (wf_type_e) : ",5X,A3)') wvfn_name( wf_type_e_new )
         endif
         if (wf_type_p_new.ne.0) then
            write(*,'(3X,"Wave function type for positrons          (wf_type_p) : ",5X,A3)') wvfn_name( wf_type_p_new )
         else
            write(*,'(3X,"No fermionic wave function for positrons  (wf_type_p) : ",5X,A3)') wvfn_name( wf_type_p_new )
         endif
         call write_empty_line(stdout,0,mpi_rank)
         write(*,'(3X,"               ________________________________                 ")')
         write(*,'(3X,"        Initialize fermionic part of the wave function          ")')
         call write_empty_line(stdout,0,mpi_rank)
      endif
      if ( wf_type_e.lt.0_int32 ) then
         spin = 'U'
      else
         spin = 'R'
      endif
      if (wf_type_e_new.ne.0_int32) then
         select case ( abs(wf_type_e_new) )
          case(1)
            if( mpi_rank.eq.0_int32 ) &
            & write(*,'(3X,"Initialize Slater det. coefficients for electrons ...     ")',advance='no')
            call ini_sldcfs( sld_c_e , n_el_u, n_el_d, edl_opt, n_par_det_e, n_par_det_e_s )
          case(2)
            if( mpi_rank.eq.0_int32 ) &
            & write(*,'(3X,"Initialize Singlet Geminal coefficients for electrons ... ")',advance='no')
            call ini_sgmcfs( sgm_c_e, n_el_u, n_el_d, edl_opt, n_par_det_e, n_par_det_e_s )
          case(3)
            if( mpi_rank.eq.0_int32 ) &
            & write(*,'(3X,"Initialize Triplet Geminal coefficients for electrons ... ")',advance='no')
            call ini_tgmcfs( tgm_c_e, n_el_u, n_el_d, edl_opt, n_par_det_e, n_par_det_e_s )
          case(4)
            if( mpi_rank.eq.0_int32 ) &
            & write(*,'(3X,"Initialize Pfaffian coefficients for electrons ...        ")',advance='no')
            call ini_pffcfs( pff_c_e, n_el, n_el_d, edl_opt, n_par_det_e, n_par_det_e_s )
          case default
         end select 
         if ( mpi_rank.eq.0_int32 )  write(*,'("(done)")')
      endif
      if ( wf_type_p_new.ne.0 ) then
         select case ( abs(wf_type_p_new) )
          case(1)
            if( mpi_rank.eq.0_int32 ) &
            & write(*,'(3X,"Initialize Slater det. coefficients for positrons ...     ")',advance='no')
            call ini_sldcfs( sld_c_p, n_po_u, n_po_d, pdl_opt, n_par_det_p, n_par_det_p_s )
          case(2)
            if( mpi_rank.eq.0_int32 ) &
            & write(*,'(3X,"Initialize Geminal coefficients for positrons ...         ")',advance='no')
            call ini_sgmcfs( sgm_c_p, n_po_u, n_po_d, pdl_opt, n_par_det_p, n_par_det_p_s )
          case(3)
            if( mpi_rank.eq.0_int32 ) &
            & write(*,'(3X,"Initialize Triplet Geminal coefficients for positrons ... ")',advance='no')
            call ini_tgmcfs( tgm_c_p, n_po_u, n_po_d, pdl_opt, n_par_det_p, n_par_det_p_s )
          case(4)
            if( mpi_rank.eq.0_int32 ) &
            & write(*,'(3X,"Initialize Pfaffian coefficients for positrons ...        ")',advance='no')
            call ini_pffcfs( pff_c_p, n_po, n_po_d, pdl_opt, n_par_det_p, n_par_det_p_s )
          case default
         end select
         if ( mpi_rank.eq.0_int32 )  write(*,'("(done)")')
      endif
   end subroutine ini_frmwvf_red
   subroutine save_frmwvf_red(wf_type_e_new, wf_type_p_new)
      use molecular_system_v, only: n_po, n_el, n_el_u, n_el_d, n_el_s, n_po_d, n_po_u, n_po_s
      use slater_determinant_params_m, only: save_sldcfs, save_sldcfs_st, sld_c_e, sld_c_p
      use singlet_geminal_params_m, only: save_sgmcfs, save_sgmcfs_st, sgm_c_e, sgm_c_p
      use triplet_geminal_params_m, only: save_tgmcfs, save_tgmcfs_st, tgm_c_e, tgm_c_p
      use pfaffian_params_m, only: save_pffcfs, save_pffcfs_st, pff_c_e, pff_c_p
      integer(int32), intent(inout) :: wf_type_e_new, wf_type_p_new
      call system('mkdir -p wvfn.save')
      if (wf_type_e_new.ne.0) then
         select case ( abs(wf_type_e_new) )
          case(1)
            call save_sldcfs( sld_c_e, n_el_u, n_el_d, 'slde' )
            call save_sldcfs_st( sld_c_e, n_el_d, 'slde' )
          case(2)
            call save_sgmcfs( sgm_c_e, n_el_s, 'sgme' )
            call save_sgmcfs_st( sgm_c_e, n_el_s, 'sgme' )
          case(3)
            call save_tgmcfs( tgm_c_e, 'tgme' )
            call save_tgmcfs_st( tgm_c_e, 'tgme' )
          case(4)
            call save_pffcfs( pff_c_e, n_el, n_el_d, 'pffe' )
            call save_pffcfs_st( pff_c_e, 'pffe' )
          case default
         end select 
      endif
      if ( wf_type_p_new.ne.0_int32 ) then
         select case ( abs(wf_type_p_new) )
          case(1)
            call save_sldcfs( sld_c_p, n_po_u, n_po_d, 'sldp' )
            call save_sldcfs_st( sld_c_p, n_po_d, 'sldp' )
          case(2)
            call save_sgmcfs( sgm_c_p, n_po_s, 'sgmp' )
            call save_sgmcfs_st( sgm_c_p, n_po_s, 'sgmp' )
          case(3)
            call save_tgmcfs( tgm_c_p, 'tgmp' )
            call save_tgmcfs_st( tgm_c_p, 'tgmp' )
          case(4)
            call save_pffcfs( pff_c_p, n_po, n_po_d, 'pffp' )
            call save_pffcfs_st( pff_c_p, 'pffp' )
          case default
         end select 
      endif
   end subroutine save_frmwvf_red
end module fermionic_wavefunction_m
