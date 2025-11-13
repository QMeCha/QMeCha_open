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
module cusps_factors_m
   use fortran_kinds_v, only: dp, int32, stdout
   use openmp_mpi_m, only: mpi_rank
   use write_lines_m
   use quantum_monte_carlo_v, only: n_wlk_max
   use molecular_system_v, only: n_el, n_el_u, n_el_d, n_po, n_po_u, n_po_d, n_fe, f_spin
   use fermionic_config_c, only: frm_cnf
   use jastrow_factors_v, only: jc1_opt, ejc1_opt, pjc1_opt, jc2_opt, ejc2_opt, pjc2_opt, jce_opt, &
   & n_par_jstc, jc2_prs
   use cusps_functions_params_m, only: b_en, b_pn, b_ee, b_pp, b_ep, jcen_prs, jcpn_prs, jcee_prs, jcpp_prs, &
   & jcep_prs, n_par_jcen, n_par_jcpn, n_par_jcee, n_par_jcpp, n_par_jcep
   use fermions_nuclei_cusps_m, only: init_fermions_nuclei_cusps, comp_fermions_nuclei_cusps, &
   & variation_fermions_nuclei_cusps, DD2_fermions_nuclei_cusps, new_DD2_fermions_nuclei_cusps, &
   & updt_DD2_fermions_nuclei_cusps, db_fermions_nuclei_cusps
   use fermions_fermions_cusps_m, only: init_fermions_fermions_cusps, comp_fermions_fermions_cusps, &
   & variation_fermions_fermions_cusps, DD2_fermions_fermions_cusps, new_DD2_fermions_fermions_cusps, &
   & updt_DD2_fermions_fermions_cusps, db_fermions_fermions_cusps
   use electrons_positrons_cusps_m, only: init_electrons_positrons_cusps, comp_electrons_positrons_cusps, &
   & variation_electrons_positrons_cusps, DD2_electrons_positrons_cusps, new_DD2_electrons_positrons_cusps, &
   & updt_DD2_electrons_positrons_cusps, db_electrons_positrons_cusps
   implicit none
   real(dp), public, save, allocatable, dimension(:,:)   :: vec_jcee, vec_jcpp
   real(dp), public, save, allocatable, dimension(:,:,:) :: vec_jcep
   real(dp), public, save, allocatable, dimension(:,:)   :: vec_jcff_new
   real(dp), public, save, allocatable, dimension(:,:)   :: DD2_en, DD2_pn
   real(dp), public, save, allocatable, dimension(:,:)   :: DD2_fn_new
   real(dp), public, save, allocatable, dimension(:,:,:) :: DD2_ee, DD2_pp, DD2_ep
   real(dp), public, save, allocatable, dimension(:,:)   :: DD2_ff_new
   real(dp), public, save, allocatable, dimension(:,:,:) :: db_en_c, db_pn_c
   real(dp), public, save, allocatable, dimension(:,:,:) :: db_ee_c, db_pp_c, db_ep_c
   public :: init_cusps, comp_cusps, variation_1body_cusps, variation_2body_cusps, &
   & updt_cusps, DD2_cusps, new_DD2_cusps, updt_DD2_cusps, dp_cusps
contains
   subroutine init_cusps()
      call init_fermions_nuclei_cusps( b_en, jcen_prs, n_par_jcen, ejc1_opt, jce_opt )
      call init_fermions_nuclei_cusps( b_pn, jcpn_prs, n_par_jcpn, pjc2_opt, jce_opt )
      call init_fermions_fermions_cusps( n_el, n_el_u, n_el_d, b_ee, jcee_prs, n_par_jcee, ejc2_opt, jce_opt )
      call init_fermions_fermions_cusps( n_po, n_po_u, n_po_d, b_pp, jcpp_prs, n_par_jcpp, pjc2_opt, jce_opt )
      call init_electrons_positrons_cusps( n_el, n_po, b_ep, jcep_prs, n_par_jcep, pjc1_opt, jce_opt )
      if ( n_par_jcen .eq. 0_int32) ejc1_opt = .false.
      if ( n_par_jcee .eq. 0_int32) ejc2_opt = .false.
      if ( n_par_jcep .eq. 0_int32) pjc1_opt = .false.
      if ( (n_par_jcpn+n_par_jcpp) .eq. 0_int32) pjc2_opt = .false.
      if ( .not.ejc1_opt .and. .not.ejc2_opt .and. .not.pjc1_opt .and. .not.pjc2_opt) jce_opt = .false.
      n_par_jstc  = n_par_jcee + n_par_jcpp + n_par_jcep + n_par_jcen + n_par_jcpn
      if ( jcee_prs.or.jcep_prs.or.jcpp_prs ) then
         jc2_prs = .true.
      else
         jc2_prs = .false.
      endif
      call write_variable_line(stdout,0,mpi_rank,2,"el-nu cusp present", jcen_prs,var_name="jcen_prs")
      call write_variable_line(stdout,0,mpi_rank,2,"po-nu cusp present", jcpn_prs,var_name="jcpn_prs")
      call write_variable_line(stdout,0,mpi_rank,2,"el-el cusp present", jcee_prs,var_name="jcee_prs")
      call write_variable_line(stdout,0,mpi_rank,2,"po-po cusp present", jcpp_prs,var_name="jcpp_prs")
      call write_variable_line(stdout,0,mpi_rank,2,"el-po cusp present ", jcep_prs,var_name="jcep_prs")
      if ( jcen_prs ) then
         allocate( DD2_en(1:4*n_el,1:n_wlk_max) ) ; DD2_en = 0.0_dp
      endif
      if ( jcpn_prs ) then
         allocate( DD2_pn(1:4*n_po,1:n_wlk_max) ) ; DD2_pn = 0.0_dp
      endif
      if ( jcen_prs.or.jcpn_prs ) then
         allocate( DD2_fn_new(1:4,1:n_wlk_max) )      ; DD2_fn_new = 0.0_dp
      endif
      if ( jcee_prs ) then
         allocate( vec_jcee(1:n_el*(n_el-1)/2,1:n_wlk_max) ) ; vec_jcee = 0.0_dp
         allocate( DD2_ee(1:4*n_el,1:n_el,1:n_wlk_max) ) ; DD2_ee = 0.0_dp
      endif
      if ( jcpp_prs ) then
         allocate( vec_jcpp(1:n_po*(n_po-1)/2,1:n_wlk_max) ) ; vec_jcpp = 0.0_dp
         allocate( DD2_pp(1:4*n_po,1:n_po,1:n_wlk_max) ) ; DD2_pp = 0.0_dp
      endif
      if ( jcep_prs ) then
         allocate( vec_jcep(1:n_po,1:n_el,1:n_wlk_max) ) ; vec_jcep = 0.0_dp
         allocate( DD2_ep(1:4*n_el,1:n_po,1:n_wlk_max) ) ; DD2_ep = 0.0_dp
      endif
      if ( jc2_prs ) then
         allocate( vec_jcff_new(1:n_fe,1:n_wlk_max) )    ; vec_jcff_new = 0.0_dp
         allocate( DD2_ff_new(1:4*n_fe,1:n_wlk_max) )      ; DD2_ff_new = 0.0_dp
      endif
      if ( ejc1_opt) then
         allocate( db_en_c(1:n_par_jcen,1:n_el,n_wlk_max) ) ; db_en_c = 0.0_dp
      endif
      if ( ejc2_opt ) then
         allocate( db_ee_c(1:n_par_jcee,1:n_el*(n_el-1)/2,n_wlk_max) ) ; db_ee_c = 0.0_dp
      endif
      if ( pjc1_opt ) then
         allocate( db_ep_c(1:n_par_jcep,1:n_el*n_po,n_wlk_max) ) ; db_ep_c = 0.0_dp
      endif
      if (pjc2_opt) then
         if( jcpn_prs ) then
            allocate( db_pn_c(1:n_par_jcpn,1:n_po,n_wlk_max) ) ; db_pn_c = 0.0_dp
         endif
         if ( jcpp_prs ) then
            allocate( db_pp_c(1:n_par_jcpp,1:n_po*(n_po-1)/2,n_wlk_max) ) ; db_pp_c = 0.0_dp
         endif
      endif
   end subroutine init_cusps
   subroutine comp_cusps( iw, vec_1b, jst )
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(inout) :: vec_1b(1:n_fe)
      real(dp),        intent(inout) :: jst
      if( jcen_prs ) then
         call comp_fermions_nuclei_cusps( -1.0_dp, n_el, frm_cnf(iw)%d_fn(:,1:n_el), b_en, vec_1b(1:n_el) )
      endif
      if ( jcpn_prs ) then
         call comp_fermions_nuclei_cusps(  1.0_dp, n_po, frm_cnf(iw)%d_fn(:,n_el+1:n_fe), b_pn, vec_1b(n_el+1:n_fe) )
      endif
      if (jcpn_prs.or.jcen_prs) jst = jst + sum(vec_1b(:))
      if ( jcee_prs ) then
         call comp_fermions_fermions_cusps( n_el, frm_cnf(iw)%d_ee, f_spin(1:n_el), b_ee, vec_jcee(:,iw) )
         jst = jst + sum(vec_jcee(:,iw))
      endif
      if ( jcpp_prs ) then
         call comp_fermions_fermions_cusps( n_po, frm_cnf(iw)%d_pp, f_spin(n_el+1:n_fe), b_pp, vec_jcpp(:,iw) )
         jst = jst + sum(vec_jcpp(:,iw))
      endif
      if ( jcep_prs ) then
         call comp_electrons_positrons_cusps( n_el, n_po, frm_cnf(iw)%d_ep, b_ep, vec_jcep(:,:,iw) )
         jst = jst + sum(vec_jcep(:,:,iw))
      endif
   end subroutine comp_cusps
   subroutine variation_1body_cusps( iw, vec_1b, new_1b, jst_var_1b )
      integer(int32), intent(in)    :: iw
      real(dp),       intent(inout) :: new_1b, vec_1b(1:n_fe)
      real(dp),       intent(inout) :: jst_var_1b
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         if ( jcen_prs ) then
            call variation_fermions_nuclei_cusps( -1.0_dp, frm_cnf(iw)%d_fn_new, b_en, new_1b )
            jst_var_1b = jst_var_1b + new_1b - vec_1b(frm_cnf(iw)%i_fe)
         endif
         if ( jcep_prs ) then
            call variation_electrons_positrons_cusps( n_po, frm_cnf(iw)%d_ff_new(n_el+1:n_fe), b_ep, vec_jcff_new(n_el+1:n_fe,iw) )
            jst_var_1b = jst_var_1b + sum(vec_jcff_new(n_el+1:n_fe,iw)) - sum(vec_jcep(:,frm_cnf(iw)%i_fe,iw))
         endif
      else
         if ( jcep_prs ) then
            call variation_electrons_positrons_cusps( n_el, frm_cnf(iw)%d_ff_new(1:n_el), b_ep, vec_jcff_new(1:n_el,iw) )
            jst_var_1b = jst_var_1b + sum(vec_jcff_new(1:n_el,iw)) - sum(vec_jcep(frm_cnf(iw)%i_fe-n_el,:,iw))
         endif
      endif
   end subroutine variation_1body_cusps
   subroutine variation_2body_cusps( iw, vec_1b, new_1b, jst_var_2b )
      integer(int32), intent(in)    :: iw
      real(dp),       intent(inout) :: new_1b, vec_1b(1:n_fe)
      real(dp),       intent(inout) :: jst_var_2b
      integer(int32) :: i1, i2, ip
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         if ( jcee_prs ) then
            call variation_fermions_fermions_cusps( frm_cnf(iw)%i_fe, n_el, frm_cnf(iw)%d_ff_new(1:n_el), &
            & f_spin(1:n_el), b_ee, vec_jcff_new(1:n_el,iw)  )
            jst_var_2b = jst_var_2b + sum(vec_jcff_new(1:n_el,iw))
            do i1 = 1, frm_cnf(iw)%i_fe - 1
               i2 = n_el*(i1-1) - i1*(i1+1)/2 + frm_cnf(iw)%i_fe
               jst_var_2b = jst_var_2b - vec_jcee(i2,iw)
            enddo
            i2 = n_el*(frm_cnf(iw)%i_fe-1)-frm_cnf(iw)%i_fe*(frm_cnf(iw)%i_fe+1)/2 + frm_cnf(iw)%i_fe
            do i1 = frm_cnf(iw)%i_fe + 1, n_el
               i2 = i2 + 1_int32
               jst_var_2b = jst_var_2b - vec_jcee(i2,iw)
            enddo
         endif
      else
         if ( jcpn_prs ) then
            new_1b = 0.0_dp
            call variation_fermions_nuclei_cusps( 1.0_dp, frm_cnf(iw)%d_fn_new, b_pn, new_1b )
            jst_var_2b = jst_var_2b + new_1b - vec_1b(frm_cnf(iw)%i_fe)
         endif
         if ( jcpp_prs ) then
            ip = frm_cnf(iw)%i_fe-n_el
            call variation_fermions_fermions_cusps( ip, n_po, frm_cnf(iw)%d_ff_new(n_el+1:n_fe), &
            & f_spin(n_el+1:n_fe), b_pp, vec_jcff_new(n_el+1:n_fe,iw)  )
            jst_var_2b = jst_var_2b + sum(vec_jcff_new(n_el+1:n_fe,iw))
            do i1 = 1, ip - 1
               i2 = n_po*(i1-1) - i1*(i1+1)/2 + ip
               jst_var_2b = jst_var_2b - vec_jcpp(i2,iw)
            enddo
            i2 = n_po*(ip-1)-ip*(ip+1)/2 + ip
            do i1 = ip + 1, n_po
               i2 = i2 + 1_int32
               jst_var_2b = jst_var_2b - vec_jcpp(i2,iw)
            enddo
         endif
      endif
   end subroutine variation_2body_cusps
   subroutine updt_cusps( iw )
      integer(int32),  intent(in)    :: iw
      integer(int32) :: i1, i2, ip
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         if (jcee_prs) then
            do i1 = 1, frm_cnf(iw)%i_fe - 1
               i2 = n_el*(i1-1) - i1*(i1+1)/2 + frm_cnf(iw)%i_fe
               vec_jcee(i2,iw) = vec_jcff_new(i1,iw)
            enddo
            i2 = n_el*(frm_cnf(iw)%i_fe-1)-frm_cnf(iw)%i_fe*(frm_cnf(iw)%i_fe+1)/2 + frm_cnf(iw)%i_fe
            do i1 = frm_cnf(iw)%i_fe + 1, n_el
               i2 = i2 + 1_int32
               vec_jcee(i2,iw) = vec_jcff_new(i1,iw)
            enddo
         endif
         if (jcep_prs) then
            vec_jcep(:,frm_cnf(iw)%i_fe,iw) = vec_jcff_new(n_el+1:n_fe,iw)
         endif
      else
         ip = frm_cnf(iw)%i_fe-n_el
         if (jcep_prs) then
            vec_jcep(ip,:,iw) = vec_jcff_new(1:n_el,iw)
         endif
         if (jcpp_prs) then
            do i1 = 1, ip - 1
               i2 = n_po*(i1-1) - i1*(i1+1)/2 + ip
               vec_jcpp(i2,iw) = vec_jcff_new(i1+n_el,iw)
            enddo
            i2 = n_po*(ip-1)-ip*(ip+1)/2 + ip
            do i1 = ip + 1, n_po
               i2 = i2 + 1_int32
               vec_jcpp(i2,iw) = vec_jcff_new(i1+n_el,iw)
            enddo
         endif
      endif
   end subroutine updt_cusps
   subroutine DD2_cusps( iw, DD2_jst )
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(inout) :: DD2_jst(1:4*n_fe)
      integer(int32) :: i1, i2
      if( jcen_prs ) then
         DD2_en(1:4*n_el,iw) = 0.0_dp
         call DD2_fermions_nuclei_cusps( -1.0_dp, n_el, frm_cnf(iw)%d_fn(:,1:n_el), b_en, DD2_en(1:4*n_el,iw) )
         DD2_jst(1:4*n_el) = DD2_en(1:4*n_el,iw)
      endif 
      if( jcpn_prs ) then
         DD2_pn(1:4*n_po,iw) = 0.0_dp
         call DD2_fermions_nuclei_cusps(  1.0_dp, n_po, frm_cnf(iw)%d_fn(:,n_el+1:n_fe), b_pn, DD2_pn(1:4*n_po,iw) )
         DD2_jst(4*n_el+1:4*n_fe) = DD2_pn(1:4*n_po,iw)
      endif 
      if( jcee_prs ) then
         call DD2_fermions_fermions_cusps( n_el, frm_cnf(iw)%d_ee, f_spin(1:n_el), b_ee, DD2_ee(:,:,iw) )
         do i1 = 1_int32, n_el
            DD2_jst(1:4*n_el) = DD2_jst(1:4*n_el) + DD2_ee(:,i1,iw)
         enddo 
      endif
      if( jcpp_prs ) then
         call DD2_fermions_fermions_cusps( n_po, frm_cnf(iw)%d_pp, f_spin(n_el+1:n_fe), b_pp, DD2_pp(:,:,iw) )
         do i1 = 1_int32, n_po
            DD2_jst(4*n_el+1:4*n_fe) = DD2_jst(4*n_el+1:4*n_fe) + DD2_pp(:,i1,iw)
         enddo 
      endif
      if( jcep_prs ) then
         call DD2_electrons_positrons_cusps( n_el, n_po, frm_cnf(iw)%d_ep, b_ep, DD2_ep(:,:,iw) )
         do i1 = 1_int32, n_po
            DD2_jst(1:4*n_el) = DD2_jst(1:4*n_el) + DD2_ep(:,i1,iw)
            do i2 = 1_int32, n_el
               DD2_jst(4*(i1+n_el)-3:4*(i1+n_el)-1) = DD2_jst(4*(i1+n_el)-3:4*(i1+n_el)-1) - DD2_ep(4*i2-3:4*i2-1,i1,iw)
               DD2_jst(4*(i1+n_el))                 = DD2_jst(4*(i1+n_el))                 + DD2_ep(4*i2,i1,iw)
            enddo 
         enddo 
      endif 
   end subroutine DD2_cusps
   subroutine new_DD2_cusps( iw, new_DD2_jst )
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(inout) :: new_DD2_jst(1:4)
      integer(int32) :: i1
      if (allocated(DD2_ff_new) ) DD2_ff_new(:,iw) = 0.0_dp
      if (allocated(DD2_fn_new) ) DD2_fn_new(:,iw) = 0.0_dp
      if (frm_cnf(iw)%i_fe.le.n_el) then
         if( jcen_prs ) then
            call new_DD2_fermions_nuclei_cusps ( -1.0_dp, frm_cnf(iw)%d_fn_new, b_en, DD2_fn_new(:,iw) )
            new_DD2_jst(:) = DD2_fn_new(:,iw)
         endif
         if( jcee_prs ) then
            call new_DD2_fermions_fermions_cusps( frm_cnf(iw)%i_fe, n_el, frm_cnf(iw)%d_ff_new(1:n_el), &
            & f_spin(1:n_el), b_ee, DD2_ff_new(1:4*n_el,iw) )
         endif 
         if ( jcep_prs ) then
            call new_DD2_electrons_positrons_cusps( n_po, frm_cnf(iw)%d_ff_new(n_el+1:n_fe), b_ep, DD2_ff_new(4*n_el+1:4*n_fe,iw) )
         endif 
         if (jcee_prs.or.jcep_prs) then
            do i1 = 1_int32, n_fe
               new_DD2_jst(1:4) = new_DD2_jst(1:4) + DD2_ff_new(4*i1-3:4*i1,iw)
            enddo
         endif
      else
         if( jcpn_prs ) then
            call new_DD2_fermions_nuclei_cusps ( 1.0_dp, frm_cnf(iw)%d_fn_new, b_pn, DD2_fn_new(:,iw) )
            new_DD2_jst(:) = DD2_fn_new(:,iw)
         endif 
         if ( jcep_prs ) then
            call new_DD2_electrons_positrons_cusps( n_el, frm_cnf(iw)%d_ff_new(1:n_el), b_ep, DD2_ff_new(1:4*n_el,iw) )
         endif 
         if( jcpp_prs ) then
            call new_DD2_fermions_fermions_cusps( frm_cnf(iw)%i_fe-n_el, n_po, frm_cnf(iw)%d_ff_new(n_el+1:n_fe),&
            & f_spin(n_el+1:n_fe), b_pp, DD2_ff_new(4*n_el+1:4*n_fe,iw) )
         endif 
         if ( jcpp_prs.or.jcep_prs ) then
            do i1 = 1_int32, n_fe
               new_DD2_jst(1:4) = new_DD2_jst(1:4) + DD2_ff_new(4*i1-3:4*i1,iw)
            enddo
         endif
      endif
   end subroutine new_DD2_cusps
   subroutine updt_DD2_cusps( iw, DD2_jst )
      integer(int32),  intent(in)    :: iw
      real(dp),        intent(inout) :: DD2_jst(1:4*n_fe)
      integer(int32) :: i1, i2, ie
      if (frm_cnf(iw)%i_fe.le.n_el) then
         ie = frm_cnf(iw)%i_fe
         if( jcen_prs ) DD2_en(4*frm_cnf(iw)%i_fe-3:4*frm_cnf(iw)%i_fe,iw) = DD2_fn_new(:,iw)
         if ( jcee_prs ) call updt_DD2_fermions_fermions_cusps( frm_cnf(iw)%i_fe, n_el, DD2_ff_new(1:4*n_el,iw), DD2_ee(:,:,iw) )
         if ( jcep_prs ) call updt_DD2_electrons_positrons_cusps( frm_cnf(iw)%i_fe, n_po, n_el, n_po, DD2_ff_new(4*n_el+1:4*n_fe,iw), DD2_ep(:,:,iw) )
         if( jcen_prs ) then
            DD2_jst(1:4*(ie-1)) = DD2_en(1:4*(ie-1),iw)
            DD2_jst(4*ie+1:4*n_el) = DD2_en(4*ie+1:4*n_el,iw)
         endif 
         if( jcpn_prs ) DD2_jst(4*n_el+1:4*n_fe) = DD2_pn(1:4*n_po,iw)
         if( jcee_prs ) then
            do i1 = 1_int32, n_el
               DD2_jst(1:4*(ie-1)) = DD2_jst(1:4*(ie-1)) + DD2_ee(1:4*(ie-1),i1,iw)
               DD2_jst(4*ie+1:4*n_el) = DD2_jst(4*ie+1:4*n_el) + DD2_ee(4*ie+1:4*n_el,i1,iw)
            enddo 
         endif
         if( jcpp_prs ) then
            do i1 = 1_int32, n_po
               DD2_jst(4*n_el+1:4*n_fe) = DD2_jst(4*n_el+1:4*n_fe) + DD2_pp(:,i1,iw)
            enddo 
         endif
         if( jcep_prs ) then
            do i1 = 1_int32, n_po
               DD2_jst(1:4*(ie-1)) = DD2_jst(1:4*(ie-1)) + DD2_ep(1:4*(ie-1),i1,iw)
               DD2_jst(4*ie+1:4*n_el) = DD2_jst(4*ie+1:4*n_el) + DD2_ep(4*ie+1:4*n_el,i1,iw)
               do i2 = 1_int32, n_el
                  DD2_jst(4*(i1+n_el)-3:4*(i1+n_el)-1) = DD2_jst(4*(i1+n_el)-3:4*(i1+n_el)-1) - DD2_ep(4*i2-3:4*i2-1,i1,iw)
                  DD2_jst(4*(i1+n_el))                 = DD2_jst(4*(i1+n_el))                 + DD2_ep(4*i2,i1,iw)
               enddo 
            enddo 
         endif 
      else
         ie = frm_cnf(iw)%i_fe - n_el
         if( jcpn_prs ) DD2_pn(4*ie-3:4*ie,iw) = DD2_fn_new(:,iw)
         if ( jcep_prs ) call updt_DD2_electrons_positrons_cusps( frm_cnf(iw)%i_fe, n_el, n_el, n_po, DD2_ff_new(1:4*n_el,iw), DD2_ep(:,:,iw) )
         if ( jcpp_prs ) call updt_DD2_fermions_fermions_cusps( ie, n_po, DD2_ff_new(4*n_el+1:4*n_fe,iw), DD2_pp(:,:,iw) )
         if( jcen_prs ) then
            DD2_jst(1:4*n_el) = DD2_en(1:4*n_el,iw)
         endif 
         if( jcpn_prs ) then
            DD2_jst(4*n_el+1:4*(frm_cnf(iw)%i_fe-1)) = DD2_pn(1:4*(ie-1),iw)
            DD2_jst(4*frm_cnf(iw)%i_fe+1:4*n_fe) = DD2_pn(4*ie+1:4*n_po,iw)
         endif 
         if( jcee_prs ) then
            do i1 = 1_int32, n_el
               DD2_jst(1:4*n_el) = DD2_jst(1:4*n_el) + DD2_ee(:,i1,iw)
            enddo 
         endif
         if( jcpp_prs ) then
            do i1 = 1_int32, n_po
               DD2_jst(4*n_el+1:4*(frm_cnf(iw)%i_fe-1)) = DD2_jst(4*n_el+1:4*(frm_cnf(iw)%i_fe-1)) + DD2_pp(1:4*(ie-1),i1,iw)
               DD2_jst(4*frm_cnf(iw)%i_fe+1:4*n_fe) = DD2_jst(4*frm_cnf(iw)%i_fe+1:4*n_fe) + DD2_pp(4*ie+1:4*n_po,i1,iw)
            enddo 
         endif
         if( jcep_prs ) then
            do i1 = 1_int32, n_po
               DD2_jst(1:4*n_el) = DD2_jst(1:4*n_el) + DD2_ep(:,i1,iw)
            enddo
            if (ie-1.gt.0_int32) then
               do i1 = 1_int32, ie-1
                  do i2 = 1_int32, n_el
                     DD2_jst(4*(i1+n_el)-3:4*(i1+n_el)-1) = DD2_jst(4*(i1+n_el)-3:4*(i1+n_el)-1) - DD2_ep(4*i2-3:4*i2-1,i1,iw)
                     DD2_jst(4*(i1+n_el))                 = DD2_jst(4*(i1+n_el))                 + DD2_ep(4*i2,i1,iw)
                  enddo 
               enddo 
            endif
            if (ie+1.le.n_po) then
               do i1 = ie+1, n_po
                  do i2 = 1_int32, n_el
                     DD2_jst(4*(i1+n_el)-3:4*(i1+n_el)-1) = DD2_jst(4*(i1+n_el)-3:4*(i1+n_el)-1) - DD2_ep(4*i2-3:4*i2-1,i1,iw)
                     DD2_jst(4*(i1+n_el))                 = DD2_jst(4*(i1+n_el))                 + DD2_ep(4*i2,i1,iw)
                  enddo 
               enddo 
            endif 
         endif
      endif
   end subroutine updt_DD2_cusps
   subroutine dp_cusps( iw, dp_ln_jst )
      integer(int32),                    intent(in) :: iw
      real(dp), dimension(n_par_jstc), intent(inout) :: dp_ln_jst
      integer(int32) :: ip
      ip = 1_int32
      if ( ejc1_opt ) then
         call db_fermions_nuclei_cusps( -1.0_dp, n_el, frm_cnf(iw)%d_fn(:,1:n_el), b_en, jce_opt, &
         & n_par_jcen, db_en_c(:,:,iw), dp_ln_jst(ip:ip+n_par_jcen-1) )
         ip = ip + n_par_jcen
      endif
      if ( pjc1_opt ) then
         call db_electrons_positrons_cusps( n_el, n_po, frm_cnf(iw)%d_ep, b_ep, jce_opt, &
         & n_par_jcep, db_ep_c(:,:,iw), dp_ln_jst(ip:ip+n_par_jcep-1) )
         ip = ip + n_par_jcep
      endif
      if (ejc2_opt) then
         call db_fermions_fermions_cusps(n_el, f_spin(1:n_el), frm_cnf(iw)%d_ee, b_ee, jce_opt, &
         & n_par_jcee, db_ee_c(:,:,iw), dp_ln_jst(ip:ip+n_par_jcee-1)  )
         ip = ip + n_par_jcee
      endif
      if ( pjc2_opt ) then
         if (jcpn_prs) then
            call db_fermions_nuclei_cusps( 1.0_dp, n_po, frm_cnf(iw)%d_fn(:,n_el+1:n_fe), b_pn, jce_opt, &
            & n_par_jcpn, db_pn_c(:,:,iw), dp_ln_jst(ip:ip+n_par_jcpn-1) )
            ip = ip + n_par_jcpn
         endif
         if (jcpp_prs) then
            call db_fermions_fermions_cusps(n_po, f_spin(n_el+1:n_fe), frm_cnf(iw)%d_pp, b_pp, jce_opt, &
            & n_par_jcpp, db_pp_c(:,:,iw), dp_ln_jst(ip:ip+n_par_jcpp-1)  )
            ip = ip + n_par_jcpp
         endif
      endif 
   end subroutine dp_cusps
end module cusps_factors_m
