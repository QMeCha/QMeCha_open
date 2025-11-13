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
module jastrow_factor_params_m
   use fortran_kinds_v,          only: dp, int32, stdout
   use openmp_mpi_m,             only: mpi_rank
   use write_lines_m
   use fermionic_wavefunction_v, only: spin
   use molecular_system_v,       only: n_el, n_po, n_el_u, n_el_d
   use jastrow_factors_v
   use orbitalsbased_jastrow_factors_params_m
   use fermionic_nuclei_atomic_factors_m,  only: g_en_u, g_en_d, g_pn
   use fermions_fermions_atomic_factors_m, only: updt_fermions_fermions_atomic_factors_params, &
   & applysym_fermions_fermions_atomic_factors_params, n_par_dyn2b_s, n_par_dyn2b, &
   & save_fermions_fermions_atomic_factors_params, read_fermions_fermions_atomic_factors_params
   use cusps_functions_params_m, only: b_en, b_pn, b_ee, b_pp, b_ep, jcen_prs, jcpn_prs, jcee_prs, jcpp_prs, &
   & jcep_prs, n_par_jcen, n_par_jcpn, n_par_jcee, n_par_jcpp, n_par_jcep, &
   & save_cusps_functions_params, read_cusps_functions_params
   use jastrow_orbitals_m,       only: n_jorbs, jorbs_s, n_par_jorbs, jorbs_mat, n_par_jorbs_s, jc_opt, je_opt
   use jeporb_mod,               only: n_jepos, n_par_jepos, jepos_mat, n_par_jepos_s, jpoc_opt, jpoe_opt
   implicit none
   public  ::  updt_jastrow_factor_params, save_jastrow_factor_params, applysym_jastrow_factor_params, read_jastrow_factor_params
contains
   subroutine updt_jastrow_factor_params( vec_jst_fct_var )
      use fermions_nuclei_cusps_m, only: updt_fermions_nuclei_cusps_params
      use fermions_fermions_cusps_m, only: updt_fermions_fermions_cusps_params
      use electrons_positrons_cusps_m, only: updt_electrons_positrons_cusps_params
      real(dp), dimension(n_par_jstc+n_par_jstd_s+n_par_jepd_s+n_par_jepa_s), intent(in) :: vec_jst_fct_var
      integer(int32) :: i1, i2, i3
      integer(int32) :: i_par
      i_par = 0_int32
      if ( ejc1_opt ) then
         call updt_fermions_nuclei_cusps_params( b_en, jce_opt, n_par_jcen, vec_jst_fct_var(i_par+1:i_par+n_par_jcen) )
         i_par = i_par + n_par_jcen
      endif
      if ( pjc1_opt ) then
         call updt_electrons_positrons_cusps_params( b_ep, jce_opt, n_par_jcep, vec_jst_fct_var(i_par+1:i_par+n_par_jcep) )
         i_par = i_par + n_par_jcep
      endif
      if ( ejc2_opt ) then
         if (jcee_prs) then
            call updt_fermions_fermions_cusps_params( b_ee, jce_opt, n_par_jcee, vec_jst_fct_var(i_par+1:i_par+n_par_jcee) )
            i_par = i_par + n_par_jcee
         endif
      endif
      if ( pjc2_opt ) then
         if (jcpn_prs) then
            call updt_fermions_nuclei_cusps_params( b_pn, jce_opt, n_par_jcpn, vec_jst_fct_var(i_par+1:i_par+n_par_jcpn) )
            i_par = i_par + n_par_jcpn
         endif
         if (jcpp_prs) then
            call updt_fermions_fermions_cusps_params( b_pp, jce_opt, n_par_jcpp, vec_jst_fct_var(i_par+1:i_par+n_par_jcpp) )
            i_par = i_par + n_par_jcpp
         endif
      endif
      if ( ejd1_opt ) then
         do i1 = 1_int32, jorbs_s%orbs_sym_tab%n_indp ; do i2 = 1_int32, jorbs_s%orbs_sym_tab%indp(i1)%n_cnct
               i3 = abs(jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2))
               g_en_u(i3) = g_en_u(i3) + vec_jst_fct_var(i_par+i1) * dble(sign(1,jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2)))
            enddo ; enddo
         i_par = i_par + jorbs_s%orbs_sym_tab%n_indp
         if (spin_j.eq.'U'.and.n_el_d.gt.0_int32) then
            do i1 = 1_int32, jorbs_s%orbs_sym_tab%n_indp ; do i2 = 1_int32, jorbs_s%orbs_sym_tab%indp(i1)%n_cnct
                  i3 = abs(jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2))
                  g_en_d(i3) = g_en_d(i3) + vec_jst_fct_var(i_par+i1) * dble(sign(1,jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2)))
               enddo ; enddo
            i_par = i_par + jorbs_s%orbs_sym_tab%n_indp
         endif
      endif
      if ( pjd1_opt ) then
         do i1 = 1_int32, jorbs_s%orbs_sym_tab%n_indp ; do i2 = 1_int32, jorbs_s%orbs_sym_tab%indp(i1)%n_cnct
               i3 = abs(jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2))
               g_pn(i3) = g_pn(i3) + vec_jst_fct_var(i_par+i1) * dble(sign(1,jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2)))
            enddo ; enddo
         i_par = i_par + jorbs_s%orbs_sym_tab%n_indp
      endif
      if ( jd2_opt.or.jdc_opt ) then
         if (cpl_m.eq.'m') then
            if (jd2_opt) call updt_fermions_fermions_atomic_factors_params( vec_jst_fct_var(i_par+1_int32:i_par+n_par_dyn2b_s)  )
         endif
         i_par = i_par + n_par_dyn2b_s
      endif
   end subroutine updt_jastrow_factor_params
   subroutine applysym_jastrow_factor_params( dp_ln_jst, dp_ln_jst_s )
      real(dp), dimension(n_par_jstc+n_par_jstd+n_par_jepd+n_par_jepa),   intent(in)    :: dp_ln_jst
      real(dp), dimension(n_par_jstc+n_par_jstd_s+n_par_jepd_s+n_par_jepa_s), intent(inout) :: dp_ln_jst_s
      integer(int32) :: i1, i2, i3, i4
      if ( n_par_jstc.gt.0_int32 ) then
         dp_ln_jst_s(1:n_par_jstc) = dp_ln_jst_s(1:n_par_jstc) + dp_ln_jst(1:n_par_jstc)
         i3 = n_par_jstc
         i4 = n_par_jstc
      else
         i3 = 0_int32
         i4 = 0_int32
      endif
      if ( ejd1_opt ) then
         do i1 = 1_int32, jorbs_s%orbs_sym_tab%n_indp ; do i2 = 1_int32, jorbs_s%orbs_sym_tab%indp(i1)%n_cnct
               dp_ln_jst_s(i3+i1) = dp_ln_jst_s(i3+i1) + dble(sign(1,jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2))) * &
               & dp_ln_jst(i4+abs(jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2)))
            enddo ; enddo
         i4 = i4 + n_jorbs
         i3 = i3 + jorbs_s%orbs_sym_tab%n_indp
         if (spin_j.eq.'U' .and. n_el_d.gt.0_int32) then
            do i1 = 1_int32, jorbs_s%orbs_sym_tab%n_indp ; do i2 = 1_int32, jorbs_s%orbs_sym_tab%indp(i1)%n_cnct
                  dp_ln_jst_s(i3+i1) = dp_ln_jst_s(i3+i1) + dble(sign(1,jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2))) * &
                  & dp_ln_jst(i4+abs(jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2)))
               enddo ; enddo
            i4 = i4 + n_jorbs
            i3 = i3 + jorbs_s%orbs_sym_tab%n_indp
         endif
      endif
      if ( pjd1_opt ) then
         if (n_po.gt.0_int32) then
            do i1 = 1_int32, jorbs_s%orbs_sym_tab%n_indp ; do i2 = 1_int32, jorbs_s%orbs_sym_tab%indp(i1)%n_cnct
                  dp_ln_jst_s(i3+i1) = dp_ln_jst_s(i3+i1) + dble(sign(1,jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2))) * &
                  & dp_ln_jst(i4+abs(jorbs_s%orbs_sym_tab%indp(i1)%cnct(i2)))
               enddo ; enddo
            i4 = i4 + n_jorbs
            i3 = i3 + jorbs_s%orbs_sym_tab%n_indp
         endif
      endif
      if ( jd2_opt .or. jdc_opt ) then
         if (cpl_m.eq.'m') then
            if (jd2_opt) call applysym_fermions_fermions_atomic_factors_params( dp_ln_jst(i4+1:i4+n_par_dyn2b),  dp_ln_jst_s(i3+1:i3+n_par_dyn2b_s) )
         endif
         i3 = i3 + n_par_dyn2b_s
         i4 = i4 + n_par_dyn2b
      endif
   end subroutine applysym_jastrow_factor_params
   subroutine save_jastrow_factor_params( )
      use fermionic_nuclei_atomic_factors_m, only: save_fermionic_nuclei_atomic_factors_params
      use fermions_fermions_atomic_factors_m, only: save_fermions_fermions_atomic_factors_params
      open(unit=1,file='wvfn.save/jst.sav',action='write',form='formatted',status='unknown')
      write(1,'("# One body Cusp parameters")')
      call save_cusps_functions_params( 1, b_en )
      call save_cusps_functions_params( 1, b_pn )
      write(1,'("# Two body Cusp parameters")')
      call save_cusps_functions_params( 1, b_ee )
      call save_cusps_functions_params( 1, b_pp )
      call save_cusps_functions_params( 1, b_ep )
      call save_fermionic_nuclei_atomic_factors_params(1)
      if (cpl_m.eq.'m') then
         call save_fermions_fermions_atomic_factors_params(1)
      endif
      close(1)
   end subroutine save_jastrow_factor_params
   subroutine read_jastrow_factor_params( )
      use fermionic_nuclei_atomic_factors_m, only: read_fermionic_nuclei_atomic_factors_params
      use fermions_fermions_atomic_factors_m, only: read_fermions_fermions_atomic_factors_params
      if ( mpi_rank.eq.0_int32 ) open(unit=1,file='wvfn.save/jst.sav',action='read',form='formatted',status='old')
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"l","Reading Jastrow factor save file.")
      if ( mpi_rank.eq.0_int32 ) read(1,*) ! '("# One body Cusp parameters")'
      call read_cusps_functions_params( 1, b_en )
      call read_cusps_functions_params( 1, b_pn )
      if ( mpi_rank.eq.0_int32 ) read(1,*) ! '("# Two body Cusp parameters")'
      call read_cusps_functions_params( 1, b_ee )
      call read_cusps_functions_params( 1, b_pp )
      call read_cusps_functions_params( 1, b_ep )
      call read_fermionic_nuclei_atomic_factors_params( 1 )
      if (cpl_m.eq.'m') then
         call read_fermions_fermions_atomic_factors_params(1)
      endif
      if (mpi_rank.eq.0_int32) close(1)
   end subroutine read_jastrow_factor_params
end module jastrow_factor_params_m
