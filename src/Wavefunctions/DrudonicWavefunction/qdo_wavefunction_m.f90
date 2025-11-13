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
module qdo_wavefunction_m
   use fortran_kinds_v,       only: dp, int8, int32, stdout
   use openmp_mpi_m,          only: mpi_rank
   use quantum_monte_carlo_v, only: restart, n_wlk_max
   use write_lines_m
   use wavefunction_optimization_v, only: n_par_opt, n_par_opt_s
   use fermionic_wavefunction_v,    only: n_par_frm_wvfn, n_par_frm_wvfn_s
   use fermionic_wavefunction_m,    only: upd_fermionic_wavefunction_v
   use qdo_wavefunction_c,  only: wavfunq_t
   use qdo_jastrow_m, only: init_qdocusps_par, n_qdo_jst_par, qdo_jastrow, &
   & save_jstcsp, qdo_cusp_prs, upd_par_jstcsp, qdo_cusp_opt
   use drudonic_orbitals_m,   only: qc_opt, n_par_dorbs, n_par_dorbs_s, dorbs_s
   use qdo_wavefunction_v,   only: ql_opt, wf_type_d, n_par_drd, n_par_drd_s, n_par_drd_wvfn, n_par_drd_wvfn_s
   use jstqepar_var, only: jstqe_prs, jqe_opt, n_par_jstqe
   use jstqe_mod,    only: ini_jstqe_fct, upd_jstqe_fct_par, save_jstqe_fct_par
   use qdo_system_v,   only: n_qdo, qdosys_prs
   use drudonic_orbitals_m,   only: ini_dorbs_mat, n_par_dorbs, n_par_dorbs_s
   use prdqcfs_mod,  only: ini_prdqcfs, save_prdqcfs, upd_prdqcfs, prdq_c_d
   use dipqcfs_mod,  only: ini_dipqcfs, save_dipqcfs, upd_dipqcfs, dipq_c_d
   implicit none
   type(wavfunq_t), public, save, allocatable, dimension(:), target :: wvfnq
   public  :: ini_qdowvf_par, upd_qdowvf_par, save_qdowvf_par
contains
   subroutine ini_qdowvf_par()
      character(3), external :: wvfn_name
      integer(int32) :: iw
      n_par_drd_wvfn   = 0_int32 ; n_par_drd_wvfn_s = 0_int32
      n_par_drd     = 0_int32 ; n_par_drd_s     = 0_int32
      if( .not.qdosys_prs ) wf_type_d = 0_int32
      if ( abs(wf_type_d).eq.5_int32) call ini_dorbs_mat()
      if (restart) call read_dorbs_par()
      if ( abs(wf_type_d).eq.5_int32) then
         n_par_drd_wvfn   = n_par_dorbs
         n_par_drd_wvfn_s = n_par_dorbs_s
      endif
      if( mpi_rank.eq.0_int32 ) then
         call write_empty_line(stdout,0,mpi_rank)
         write(*,'(3X,"================================================================")')
         write(*,'(3X,"               QDO WAVE FUNCTION INITIALIZATION                 ")')
         write(*,'(3X,"               ________________________________                 ")')
         write(*,'(3X,"                       Input parameters                         ")')
         call write_empty_line(stdout,0,mpi_rank)
         write(*,'(3X,"Wave type for drudons                      (wf_type_d) : ",5X,A3)') wvfn_name( wf_type_d )
         call write_empty_line(stdout,0,mpi_rank)
      endif
      select case ( abs(wf_type_d) )
       case(5)
         if (restart) then
            if( mpi_rank.eq.0_int32 ) write(*,'(3X,"Reading Product coefficients for drudons ... ")')
            call ini_prdqcfs( prdq_c_d, n_qdo, n_par_drd, n_par_drd_s, head_fle='prdq' )
         else
            if( mpi_rank.eq.0_int32 ) write(*,'(3X,"Initialize Product coefficients for drudons ... ")')
            call ini_prdqcfs( prdq_c_d, n_qdo, n_par_drd, n_par_drd_s)
         endif
       case(6)
         if (restart) then
            if( mpi_rank.eq.0_int32 ) write(*,'(3X,"Reading Dipole-dipole coefficients for drudons ... ")')
            call ini_dipqcfs( dipq_c_d, n_qdo, n_par_drd, n_par_drd_s, head_fle='dipq' )
         else
            if( mpi_rank.eq.0_int32 ) write(*,'(3X,"Initialize Dipole-dipole coefficients for drudons ... ")')
            call ini_dipqcfs( dipq_c_d, n_qdo, n_par_drd, n_par_drd_s)
         endif
       case default
      end select 
      call init_qdocusps_par()
      n_par_drd_wvfn   = n_par_drd_wvfn + n_par_drd + n_qdo_jst_par
      n_par_drd_wvfn_s = n_par_drd_wvfn_s + n_par_drd_s + n_qdo_jst_par
      if (jstqe_prs) then
         call ini_jstqe_fct()
         n_par_drd_wvfn   = n_par_drd_wvfn   + n_par_jstqe
         n_par_drd_wvfn_s = n_par_drd_wvfn_s + n_par_jstqe
      endif
      n_par_opt   = n_par_opt + n_par_drd_wvfn
      n_par_opt_s = n_par_opt_s + n_par_drd_wvfn_s
      allocate( wvfnq(1:n_wlk_max) )
      if (qdo_cusp_prs) allocate( qdo_jastrow(1:n_wlk_max) )
      do iw = 1_int32, n_wlk_max
         call wvfnq(iw)%ini( )
         if (qdo_cusp_prs) call qdo_jastrow(iw)%ini(iw)
      enddo
   end subroutine ini_qdowvf_par
   subroutine upd_qdowvf_par( vec_wvfnq_var )
      real(dp), dimension(n_par_drd_wvfn_s), intent(in) :: vec_wvfnq_var
      integer(int32) :: ip
      ip = 1_int32
      if ( ql_opt ) then
         select case ( abs(wf_type_d) )
          case(5)
            call upd_prdqcfs( prdq_c_d, n_qdo, n_par_drd_s, vec_wvfnq_var(1:n_par_drd_s) )
            ip = ip + n_par_drd_s
          case(6)
            call upd_dipqcfs( dipq_c_d, n_qdo, n_par_drd_s, vec_wvfnq_var(1:n_par_drd_s) )
            ip = ip + n_par_drd_s
          case default
         end select
      endif
      if ( (abs(wf_type_d).eq.5).and.qc_opt ) then
         call dorbs_s%upd( vec_wvfnq_var(ip:ip+n_par_dorbs_s-1), 0.001_dp )
         ip = ip + n_par_dorbs_s
      endif
      if( jstqe_prs.and.jqe_opt ) then
         call upd_jstqe_fct_par( vec_wvfnq_var(ip:ip+n_par_jstqe-1) )
         ip = ip + n_par_jstqe
      endif
      if (qdo_cusp_opt) then
         call upd_par_jstcsp( vec_wvfnq_var(ip:ip+n_qdo_jst_par-1) )
         ip = ip + n_qdo_jst_par
      endif
   end subroutine upd_qdowvf_par
   subroutine sym_qdowvf_par( O, O_s )
      real(dp), dimension(n_par_drd_wvfn),   intent(in) :: O
      real(dp), dimension(n_par_drd_wvfn_s), intent(inout) :: O_s
      integer(int32) :: ip, ip_s
      if ( ql_opt ) then
         select case ( abs(wf_type_d) )
          case(5)
            O_s(1:n_par_drd_s) = O(1:n_par_drd)
          case(6)
            O_s(1:n_par_drd_s) = O(1:n_par_drd)
          case default
         end select 
      endif
      ip = n_par_drd ; ip_s = n_par_drd_s
      if ( (abs(wf_type_d).eq.5).and.qc_opt ) then
         if ( n_par_dorbs_s.gt.0_int32) then
            O_s(ip_s+1:ip_s+n_par_dorbs_s) = O(ip+1:ip+n_par_dorbs)
            ip   = ip + n_par_dorbs
            ip_s = ip_s + n_par_dorbs_s
         endif ; endif
      if( jstqe_prs.and.jqe_opt ) then
         O_s(ip_s+1:ip_s+n_par_jstqe) = O(ip+1:ip+n_par_jstqe)
         ip = ip + n_par_jstqe
         ip_s = ip_s + n_par_jstqe
      endif
      if (qdo_cusp_opt) then
         O_s(ip_s+1:ip_s+n_qdo_jst_par) = O(ip+1:ip+n_qdo_jst_par)
         ip = ip + n_qdo_jst_par
      endif
   end subroutine sym_qdowvf_par
   subroutine save_qdowvf_par()
      call system('mkdir -p wvfn.save')
      select case ( abs(wf_type_d) )
       case(5)
         call save_prdqcfs( prdq_c_d, n_qdo, 'prdq' )
       case(6)
         call save_dipqcfs( dipq_c_d, n_qdo, 'dipq' )
       case default
      end select
      if ( (abs(wf_type_d).eq.5) ) call save_dorbs_par()
      if (jstqe_prs) call save_jstqe_fct_par()
      if (qdo_cusp_prs) call save_jstcsp( )
   end subroutine save_qdowvf_par
end module qdo_wavefunction_m
