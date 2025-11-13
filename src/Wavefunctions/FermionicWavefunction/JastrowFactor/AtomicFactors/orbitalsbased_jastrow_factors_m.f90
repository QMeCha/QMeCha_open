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
module orbitalsbased_jastrow_factors_m
   use fortran_kinds_v, only: dp, int32, stdout
   use openmp_mpi_m, only: mpi_rank
   use write_lines_m
   use molecular_system_v, only: n_fe, n_po, n_el
   use jastrow_factors_v, only: jd1_prs, ejd1_prs, pjd1_prs, jd2_prs, ejd2_prs, pjd2_prs, &
   & jd_opt,  ejd_opt,  pjd_opt, jd1_opt, ejd1_opt, pjd1_opt , n_par_jst, &
   & jd2_opt, ejd2_opt, pjd2_opt, n_par_jstd, n_par_jstd_s, n_par_jstc, &
   & n_par_jepd, n_par_jepd_s, jdep_opt, jdep_prs, jdte_prs, jdtp_prs, &
   & jqep_prs, n_par_jepa, n_par_jepa_s, jaep_opt, jaep_prs, jdc_opt
   use orbitalsbased_jastrow_factors_params_m, only: n_par_dyn2b,n_par_dyn2b_s, n_par_dyn1b, &
   & n_par_dyn1b_s, spin_j, cpl_m, cpl_j, corr_type, order_coupling, comp_fermions_fermions_atomic_factors_symtable
   use fermionic_nuclei_atomic_factors_m, only: init_fermionic_nuclei_atomic_factors, &
   & comp_fermionic_nuclei_atomic_factors, variation_fermionic_nuclei_atomic_factors, &
   & DD2_fermionic_nuclei_atomic_factors, new_DD2_fermionic_nuclei_atomic_factors, &
   & updt_DD2_fermionic_nuclei_atomic_factors, dg_fermionic_nuclei_atomic_factors, &
   & da_fermionic_nuclei_atomic_factors, da_orbitals_sums, comp_orbitals_sums, & 
   & updt_orbitals_sums
   use fermions_fermions_atomic_factors_m, only: init_fermions_fermions_atomic_factors, &
   & comp_fermions_fermions_atomic_factors, variation_fermions_fermions_atomic_factors, &
   & updt_fermions_fermions_atomic_factors, DD2_fermions_fermions_atomic_factors, &
   & new_DD2_fermions_fermions_atomic_factors, updt_DD2_fermions_fermions_atomic_factors, &
   & dg_fermions_fermions_atomic_factors, da_fermions_fermions_atomic_factors
   use jastrow_orbitals_m, only: n_jorbs, n_par_jorbs, n_par_jorbs_s, jc_opt, je_opt
   use jeporb_mod, only: n_jepos, n_par_jepos, n_par_jepos_s, jpoc_opt, jpoe_opt
   implicit none
   public :: init_orbitalsbased_jastrow_factors, comp_orbitalsbased_jastrow_factors, &
   & variation_orbitalsbased_1body_factors, variation_orbitalsbased_2body_factors, &
   & updt_orbitalsbased_jastrow_factors, DD2_orbitalsbased_jastrow_factors, &
   & new_DD2_orbitalsbased_jastrow_factors, updt_DD2_orbitalsbased_jastrow_factors, &
   & dp_orbitalsbased_jastrow_factors
contains
   subroutine init_orbitalsbased_jastrow_factors()
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Two-body dynamical factor", cpl_m,var_name="cpl_m")
      call write_variable_line(stdout,0,mpi_rank,2,"Localization of two-body coupling", cpl_j,var_name="cpl_j")
      call write_variable_line(stdout,0,mpi_rank,2,"Spin of the two-body coupling", spin_j,var_name="spin_j")
      n_par_jstd = 0_int32 ; n_par_jstd_s = 0_int32
      if (n_par_jorbs_s.eq.0_int32) then
         jc_opt = .false. ; je_opt = .false.
      endif
      if (n_par_jepos_s.eq.0_int32) then
         jpoc_opt = .false. ; jpoe_opt = .false.
      endif
      if ( n_jorbs.gt.0_int32 ) then
         jd1_prs = .true. 
         if ( n_fe.ge.2_int32 ) jd2_prs = .true.
         call init_fermionic_nuclei_atomic_factors( )
         if (jd2_prs) call comp_fermions_fermions_atomic_factors_symtable()
         if (jd2_opt) then 
            if (n_el.gt.1_int32) ejd2_opt = .true. 
            if (n_po.gt.1_int32) pjd2_opt = .true.
         endif
         if( cpl_m.eq.'m' ) then
            call init_fermions_fermions_atomic_factors( )
         endif
      else
         ejd1_opt = .false. ; ejd2_opt = .false.
         pjd1_opt = .false. ; pjd2_opt = .false.
         jd1_prs = .false. ;  jd2_prs = .false.
         jdc_opt = .false. ; jd2_opt = .false. 
         jd1_opt = .false.
      endif
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"fe-nu dynamical part", jd1_prs,var_name="jd1_prs")
      call write_variable_line(stdout,0,mpi_rank,2,"fe-fe dynamical part", jd2_prs,var_name="jd2_prs")
   end subroutine init_orbitalsbased_jastrow_factors
   subroutine comp_orbitalsbased_jastrow_factors( iw , jst )
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(inout) :: jst
      if ( jd1_prs.or.(jd2_prs.and.cpl_m.ne.'l') ) call comp_orbitals_sums( iw )
      if ( jd1_prs ) call comp_fermionic_nuclei_atomic_factors( iw, jst )
      if ( jd2_prs ) then
         if( cpl_m.eq.'m' ) then
            call comp_fermions_fermions_atomic_factors( iw, jst )
         endif
      endif
   end subroutine comp_orbitalsbased_jastrow_factors
   subroutine variation_orbitalsbased_1body_factors( iw, jst_var_1b )
      integer(int32), intent(in)    :: iw
      real(dp),       intent(inout) :: jst_var_1b
      if ( jd1_prs ) call variation_fermionic_nuclei_atomic_factors( iw, jst_var_1b )
   end subroutine variation_orbitalsbased_1body_factors
   subroutine variation_orbitalsbased_2body_factors( iw, jst_var_2b )
      integer(int32), intent(in)    :: iw
      real(dp),       intent(inout) :: jst_var_2b
      if ( jd2_prs ) then
         if( cpl_m.eq.'m' ) then
            call variation_fermions_fermions_atomic_factors( iw, jst_var_2b )
         endif
      endif ! jd_prs
   end subroutine variation_orbitalsbased_2body_factors
   subroutine updt_orbitalsbased_jastrow_factors( iw )
      integer(int32),  intent(in)    :: iw
      if ( jd1_prs .or. (jd2_prs.and.cpl_m.ne.'l') ) call updt_orbitals_sums( iw )
      if ( jd2_prs ) then
         if( cpl_m.eq.'m' ) then
            call updt_fermions_fermions_atomic_factors( iw )
         endif
      endif
   end subroutine updt_orbitalsbased_jastrow_factors
   subroutine DD2_orbitalsbased_jastrow_factors( iw, DD2_jst )
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(inout) :: DD2_jst(1:4*n_fe)
      if ( jd1_prs ) call DD2_fermionic_nuclei_atomic_factors( iw, DD2_jst(:) )
      if ( jd2_prs ) then
         if (cpl_m.eq.'m') then
            call DD2_fermions_fermions_atomic_factors( iw, DD2_jst(:) )
         endif
      endif 
   end subroutine DD2_orbitalsbased_jastrow_factors
   subroutine new_DD2_orbitalsbased_jastrow_factors( iw, new_DD2_jst )
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(inout) :: new_DD2_jst(1:4)
      if ( jd1_prs ) call new_DD2_fermionic_nuclei_atomic_factors( iw, new_DD2_jst(:) )
      if ( jd2_prs ) then
         if (cpl_m.eq.'m') then
            call new_DD2_fermions_fermions_atomic_factors( iw, new_DD2_jst(:) )
         endif
      endif 
   end subroutine new_DD2_orbitalsbased_jastrow_factors
   subroutine updt_DD2_orbitalsbased_jastrow_factors( iw, DD2_jst )
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(inout) :: DD2_jst(1:4*n_fe)
      if ( jd1_prs ) call updt_DD2_fermionic_nuclei_atomic_factors( iw, DD2_jst(:) )
      if ( jd2_prs ) then
         if (cpl_m.eq.'m') then
            call updt_DD2_fermions_fermions_atomic_factors( iw, DD2_jst(:) )
         endif
      endif ! jd2_prs
   end subroutine updt_DD2_orbitalsbased_jastrow_factors
   subroutine dp_orbitalsbased_jastrow_factors( iw, dp_ln_jst )
      integer(int32),                              intent(in) :: iw
      real(dp), dimension(n_par_jst-n_par_jstc), intent(inout) :: dp_ln_jst
      integer(int32) :: ip
      ip = 1_int32
      if ( jd1_prs) then
         call dg_fermionic_nuclei_atomic_factors( iw, dp_ln_jst(ip:ip+n_par_dyn1b-1) )
         ip = ip + n_par_dyn1b
      endif
      if ( jd2_prs) then
         if (cpl_m.eq.'m') then
            call dg_fermions_fermions_atomic_factors( iw, dp_ln_jst(ip:ip+n_par_dyn2b-1) )
         endif
         ip = ip + n_par_dyn2b
      endif
      if ( jc_opt ) then
         dp_ln_jst(ip:ip+n_par_jorbs-1) = 0.0_dp
         if ( jd1_prs .or. jd2_prs) call da_orbitals_sums( iw )
         if ( jd1_prs ) call da_fermionic_nuclei_atomic_factors( iw, dp_ln_jst(ip:ip+n_par_jorbs-1) )
         if ( jd2_prs ) then
            if (cpl_m.eq.'m')then
               call da_fermions_fermions_atomic_factors( iw, dp_ln_jst(ip:ip+n_par_jorbs-1) )
            endif
         endif 
         ip = ip + n_par_jorbs
      endif 
   end subroutine dp_orbitalsbased_jastrow_factors
end module orbitalsbased_jastrow_factors_m
