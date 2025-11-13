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
module io_datasheet_m
   use fortran_kinds_v,             only: dp, int32
   use openmp_mpi_m
   use write_lines_m,               only: write_separator_line, write_simple_line, &
   & write_empty_line, write_simple_strn_line
   use wavefunction_optimization_v, only: da, eps, eps_auto_tune, n_opt_step, opt_mthd,&
   & corr_samp, n_crs_step, norm_cut, n_par_opt, n_par_opt_s, force_cut, lin_solv_mthd,&
   & eps_force_reg, cg_acc, n_avg_step, err_fct, norm_min, alpha_norm
   use quantum_monte_carlo_v,       only: restart, save_step, save_status, qmc_mthd, &
   & n_bra, n_bra_vmc, n_bin, bin_l, n_trm, n_scr, dt, dt_e, dt_p, dt_d, var_dt, smp_type,  &
   & arate, n_wlk, n_wlk_max, nuc_corr, a_drft, grad_updat, mapping_type, sysname, &
   & spc_map, smp_mthd, two_step_sampling
   use diffusion_monte_carlo_v,     only: e_r, g_r, g_e, brnch_type, encrr_type,&
   & max_wlk_life, w_chk, w_max, n_trm_dmc, trotter_order
   use jastrow_factors_v,                  only: jc1_opt, jc2_opt, jd_opt, jd1_opt, jd2_opt,&
   & jce_opt, ejc1_opt, ejc2_opt, jd_opt, ejd1_opt, ejd2_opt, jaep_opt, &
   & pjc1_opt, pjc2_opt, pjd1_opt, pjd2_opt, jdep_opt, jdc_opt
   use orbitalsbased_jastrow_factors_params_m, only: cpl_j, cpl_m, spin_j, corr_type
   use jstqepar_var,                only: jqe_opt, jstqe_prs
   use fermionic_wavefunction_v,    only: dl_opt, pdl_opt, edl_opt, n_par_frm_wvfn, n_par_frm_wvfn_s
   use jastrow_orbitals_m,          only: je_opt, jc_opt
   use fermionic_orbitals_m,        only: aoe_opt, aoc_opt
   use feporb_mod,                  only: poe_opt, poc_opt
   use jeporb_mod,                  only: jpoe_opt, jpoc_opt
   use molec_symmetry_m,            only: atm_rot, sys_rot, atm_eqv, atm_trs
   use pseudopotentials_v,          only: fle_pspot
   use molecular_system_v,          only: fle_coord, molsys_prs
   use mersenne_twister19937_m,     only: rndm_seed
   use electronic_properties_m,     only: dip_mom, quad_mom, expec_dist, expec_pote
   use qdo_wavefunction_v,          only: ql_opt, n_par_drd_wvfn, n_par_drd_wvfn_s
   use drudonic_orbitals_m,         only: qc_opt, qe_opt
   use qdo_system_v,                only: fle_coord_qdo, qdosys_prs
   use basisset_v,                  only: fle_basis_qdo, fle_basis
   implicit none
   logical :: de_opt, dc_opt, el_opt, po_opt
   character(len=100), public :: fle_dataio
   namelist /qmcsmp/ qmc_mthd, bin_l, n_trm, n_bra, n_bin, dt, dt_e, dt_p, dt_d, var_dt, rndm_seed, &
   & n_scr, smp_mthd, spc_map, arate, n_wlk, nuc_corr, a_drft, mapping_type, sysname, two_step_sampling
   namelist /dmc/ e_r, g_r, g_e, brnch_type, encrr_type, max_wlk_life, w_chk, n_trm_dmc, w_max, &
   & trotter_order
   namelist /optmet/ da, eps, n_opt_step, jc1_opt, jc2_opt, n_crs_step, n_avg_step, &
   & jd_opt, jd1_opt, jd2_opt, je_opt, jc_opt, jdc_opt, de_opt, dc_opt, dl_opt, &
   & ql_opt, qc_opt, qe_opt, jqe_opt, restart, jce_opt, norm_cut, err_fct, &
   & el_opt, po_opt, force_cut, lin_solv_mthd, eps_force_reg, cg_acc, eps_auto_tune, &
   & norm_min, alpha_norm
   namelist /basset/ fle_basis, fle_basis_qdo, fle_pspot
   namelist /molsys/ fle_coord, fle_coord_qdo
   namelist /elcprp/ dip_mom, quad_mom, expec_dist, expec_pote
   public :: read_data_fle, write_data_fle, dflt_data_var, bcst_data_var
contains
   subroutine dflt_data_var
      restart     = .false.
      save_step   = 0_int32
      save_status = 0_int32
      qmc_mthd   = 'vmc'
      n_wlk      = 1_int32
      n_wlk_max  = 1_int32
      n_trm  = 10_int32
      n_bin  = 100_int32
      bin_l  =   0_int32
      n_bra  =   0_int32
      n_scr  =   0_int32
      rndm_seed = 142857369_int32
      force_cut = 1.0d-8
      a_drft = 0.5_dp
      var_dt = .true.
      spc_map = .false.
      mapping_type = 2_int32
      smp_mthd = 'g3d'
      smp_type = 4_int32
      two_step_sampling = .true.
      grad_updat = .false.
      arate = 0.5_dp
      nuc_corr = .false.
      dt     = 1.000_dp
      dt_e   = 1.000_dp
      dt_p   = 1.000_dp
      dt_d   = 1.000_dp
      molsys_prs = .false.
      qdosys_prs = .false.
      sysname = 'non'
      opt_mthd   = 'src'
      corr_samp  = .false.
      da         = 0.005_dp
      n_opt_step = 0_int32
      n_crs_step = 0_int32
      n_avg_step = 1_int32
      eps        = 0.001_dp
      eps_force_reg = v
      eps_auto_tune   = .false.
      err_fct = 4.0_dp
      norm_cut   = 0.1_dp
      norm_min   = 0.001_dp
      alpha_norm   = 0.0_dp
      el_opt     = .true.
      po_opt     = .true.
      jc1_opt    = .false. ; ejc1_opt = .false. ; pjc1_opt = .false.
      jc2_opt    = .false. ; ejc2_opt = .false. ; pjc2_opt = .false.
      jce_opt    = .false.
      jd_opt     = .false. ; jdc_opt = .false.
      jd1_opt    = .false. ; ejd1_opt = .false. ; pjd1_opt = .false.
      jd2_opt    = .false. ; ejd2_opt = .false. ; pjd2_opt = .false.
      jdep_opt   = .false.
      jaep_opt   = .false.
      jc_opt     = .false. ; je_opt   = .false.
      jpoe_opt   = .false. ; jpoc_opt = .false.
      de_opt     = .false. ; aoe_opt  = .false. ; poe_opt  = .false.
      dc_opt     = .false. ; aoc_opt  = .false. ; poc_opt  = .false.
      dl_opt     = .false. ; pdl_opt  = .false. ; edl_opt  = .false.
      ql_opt     = .false.
      qc_opt     = .false.
      qe_opt     = .false.
      jqe_opt    = .false.
      lin_solv_mthd = 'lu'
      cg_acc = 0.000000001_dp
      n_par_opt   = 0_int32
      n_par_opt_s = 0_int32
      n_par_frm_wvfn = 0_int32
      n_par_frm_wvfn_s = 0_int32
      n_par_drd_wvfn = 0_int32
      n_par_drd_wvfn_s = 0_int32
      jstqe_prs   = .false.
      e_r = 0.0_dp
      g_r = 1.0_dp
      g_e = 0.2_dp
      brnch_type = 2_int32
      encrr_type = 1_int32
      max_wlk_life = 0_int32
      w_chk = 0.5_dp
      w_max = 10.0_dp
      n_trm_dmc = 10_int32
      trotter_order = 2_int32
      atm_rot = .true.
      sys_rot = .true.
      atm_eqv = .false.
      atm_trs = .false.
      fle_basis     = 'none'
      fle_basis_qdo = 'none'
      fle_pspot = 'none'
      fle_coord     = 'none'
      fle_coord_qdo = 'none'
      dip_mom    = .false.
      quad_mom   = .false.
      expec_dist = .false.
      expec_pote = .false.
      cpl_j = 'f'
      cpl_m = 'm'
      spin_j = 'R'
   end subroutine dflt_data_var
   subroutine read_data_fle( fle_name, ierr )
      character(*),   intent(in)    :: fle_name
      integer(int32), intent(inout) :: ierr
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Reading input datafile", fle_name,var_name="file_data")
      if( mpi_rank.eq.0_int32 ) then
         open(unit=1, file=fle_name, action='read', form='formatted', status='old', iostat=ierr )
         if(ierr.ne.0_int32) then
            call write_variable_line(stdout,0,mpi_rank,2,"ERROR!!! Cannot open input file", trim(fle_name))
         endif
      endif
      if( mpi_rank.eq.0_int32 .and. ierr.eq.0_int32 ) then
         read(1,nml=qmcsmp, iostat=ierr ) ; rewind(1)
         call check_namelist_mandatory(ierr, "&qmcsmp")
         if ( qmc_mthd.eq.'dmc' ) then
            read(1,nml=dmc, iostat=ierr ) ; rewind(1)
            call check_namelist_mandatory(ierr, "&dmc")
            select case ( brnch_type )
             case(1)
               n_wlk_max = n_wlk + int(n_wlk / 2_int32)
             case(2)
               n_wlk_max = n_wlk
             case default
               n_wlk_max = n_wlk
            end select
            if (encrr_type.ge.20_int32) then
               spc_map = .true.
               mapping_type = -1
            else
               spc_map = .false.
               mapping_type = 0
            endif
            smp_mthd = 'l3d'
            var_dt = .false.
         else
            n_wlk_max = n_wlk
         endif
         if (dt.gt.0.0_dp) then
            dt_e = dt
            dt_p = dt
            dt_d = dt
         endif
         if ( n_bra.eq.0_int32 .and. qmc_mthd.eq.'vmc' ) n_bra = 2_int32
         select case (smp_mthd)
          case('s1d')
            smp_type = 0
          case('g1d')
            smp_type = 1
          case('s3d')
            smp_type = 2
          case('g3d')
            smp_type = 3
          case('l3d')
            arate = 0.99_dp
            smp_type = 10
          case default
         end select
         if ( smp_type.eq.10_int32.or.qmc_mthd.eq.'dmc') then
            grad_updat = .true.
         else
            grad_updat = .false.
         endif
      endif
      if( mpi_rank.eq.0_int32 .and. ierr.eq.0_int32 ) then
         read(1,nml=basset, iostat=ierr) ; rewind(1)
         call check_namelist_optional( ierr, "&basset")
      endif
      if(mpi_rank.eq.0_int32 .and. ierr.eq.0_int32) then
         read(1,nml=molsys, iostat=ierr) ; rewind(1)
         call check_namelist_optional( ierr, "&molsys")
      endif
      if( mpi_rank.eq.0_int32 .and. ierr.eq.0_int32) then
         read(1,nml=optmet, iostat=ierr) ; rewind(1)
         if ( ierr .lt. 0_int32 ) then
            opt_mthd = 'non'
            if ( .not.restart ) restart = .true.
         endif
         call check_namelist_optional( ierr, "&optmet")
         if (n_crs_step.gt.0_int32) corr_samp = .true.
         if (opt_mthd .eq. 'ams') force_cut = 0.0d0 ! No cuts on forces (Order is important)
         if ( jd_opt ) then
            jd1_opt = .true.
            jd2_opt = .true.
         endif
         if ( qe_opt ) qc_opt = .true.
         if ( de_opt ) dc_opt = .true.
         if ( je_opt ) jc_opt = .true.
         if ( po_opt ) then
            if ( jd1_opt ) pjd1_opt = .true.
            if ( jd2_opt ) then
               pjd2_opt = .true.
               jdep_opt = .true.
               jaep_opt =  .true.
            endif
            if ( jc1_opt ) pjc1_opt = .true.
            if ( jc2_opt ) pjc2_opt = .true.
            if ( dc_opt )  poc_opt = .true.
            if ( de_opt )  poe_opt = .true.
            if ( dl_opt )  pdl_opt = .true.
            if ( jc_opt )  jpoc_opt = .true.
            if ( je_opt )  jpoe_opt = .true.
         endif
         if (el_opt) then
            if ( jd1_opt ) ejd1_opt = .true.
            if ( jd2_opt ) ejd2_opt = .true.
            if ( jc1_opt ) ejc1_opt = .true.
            if ( jc2_opt ) ejc2_opt = .true.
            if ( dc_opt )  aoc_opt = .true.
            if ( de_opt )  aoe_opt = .true.
            if ( dl_opt ) edl_opt = .true.
         endif
      endif
      if( mpi_rank.eq.0_int32  .and. ierr.eq.0_int32) then
         read(1,nml=elcprp, iostat=ierr) ; rewind(1)
         call check_namelist_optional( ierr, "&elcprp")
      endif
      if( mpi_rank.eq.0_int32 ) close(1)
   end subroutine read_data_fle
   subroutine write_data_fle( dflt, inp_type, fle_name )
      logical       , intent(in) :: dflt
      character(*),   intent(in) :: inp_type
      character(*),   intent(in) :: fle_name
      call write_variable_line(stdout,0,mpi_rank,2,"Writing input file",trim(fle_name),var_name=trim(inp_type))
      call write_empty_line(stdout,0,mpi_rank)
      if( dflt ) call dflt_data_var()
      if( mpi_rank.eq.0_int32 ) then
         open(unit=1,file=fle_name,action='write',form='formatted')
         write(1,nml=qmcsmp)
         write(1,nml=basset)
         write(1,nml=molsys)
         if(inp_type.eq.'opt') then
            write(1,nml=optmet)
         endif
         if(inp_type.eq.'dmc') then
            write(1,nml=dmc)
         endif
         write(1,nml=elcprp)
         close(1)
      endif
   end subroutine write_data_fle
   subroutine bcst_data_var
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(cpl_j, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(cpl_m, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(corr_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(spin_j, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(fle_pspot,         100, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(fle_coord,         100, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(fle_coord_qdo,     100, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(fle_basis,         100, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(fle_basis_qdo,     100, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(sysname,           100, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(lin_solv_mthd,       3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(opt_mthd, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jc1_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jc2_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(ejc1_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(ejc2_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(pjc1_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(pjc2_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jce_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jd_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jd1_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jd2_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jdc_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(ejd1_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(ejd2_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(pjd1_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(pjd2_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jdep_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jaep_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jc_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(je_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(grad_updat ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(dl_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(edl_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(pdl_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(aoe_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(aoc_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(poc_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(poe_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jpoc_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jpoe_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(ql_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(qc_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(qe_opt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jqe_opt,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(da   ,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(eps  ,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(eps_force_reg,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(cg_acc,           1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(eps_auto_tune,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_opt_step,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_crs_step,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_avg_step,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(force_cut,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(corr_samp ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(restart,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(qmc_mthd,     3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(smp_mthd,     3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_wlk,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_wlk_max,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_trm,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_bin,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(bin_l,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_bra,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_scr,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(rndm_seed,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(smp_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(two_step_sampling,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(e_r,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(g_r,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(g_e,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(a_drft,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_trm_dmc,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(w_chk ,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(w_max ,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(encrr_type,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(brnch_type,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(max_wlk_life,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(trotter_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(norm_cut, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(norm_min, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(alpha_norm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(err_fct, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(dt,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(dt_e,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(dt_p,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(dt_d,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(arate    ,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(var_dt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(spc_map , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(mapping_type,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(nuc_corr,      1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(atm_rot,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(sys_rot,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(atm_eqv,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(atm_trs,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(dip_mom,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(quad_mom,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(expec_dist, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(expec_pote, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      if ( mpierr .ne. 0 ) then
         call write_simple_line(stdout,0,mpi_rank,2,'l',"ERROR!!! While bcasting input variables")
         stop
      endif
#endif
   end subroutine bcst_data_var
end module io_datasheet_m
