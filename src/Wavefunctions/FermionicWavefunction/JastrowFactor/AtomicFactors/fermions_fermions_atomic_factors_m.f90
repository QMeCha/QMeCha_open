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
module fermions_fermions_atomic_factors_m
   use fortran_kinds_v,       only: dp, int32, stdout
   use openmp_mpi_m
   use write_lines_m
   use quantum_monte_carlo_v, only: n_wlk_max
   use orbitalsbased_jastrow_factors_params_m
   use molecular_system_v,    only: n_at, n_el, n_el_d, n_el_u, n_po, n_fe, n_sys_sym
   use fermionic_config_c,    only: frm_cnf
   use geminal_symmetries_m,  only: cmp_nzr_fee, cmp_sgm_sym, lmbd_ele_t, lmbd_sym_t
   use jastrow_factors_v,     only: jd2_prs, jd2_opt, ejd2_opt, pjd2_opt, n_par_jstd, n_par_jstd_s
   use orbbss_cls,            only: orbspa_t
   use jastrow_orbitals_m,    only: jorbs, n_jorbs, n_par_jorbs, jorbs_mat, jc_opt, je_opt
   implicit none
   real(dp), save, public, allocatable, dimension(:,:,:) :: mat_GX
   real(dp), save, public, allocatable, dimension(:,:)   :: vec_GX_s
   real(dp), save, public, allocatable, dimension(:,:)   :: vec_GX_d
   real(dp), save, public, allocatable, dimension(:,:,:) :: dg_2b
   public  :: init_fermions_fermions_atomic_factors, comp_fermions_fermions_atomic_factors, &
   & variation_fermions_fermions_atomic_factors, updt_fermions_fermions_atomic_factors, &
   & DD2_fermions_fermions_atomic_factors, new_DD2_fermions_fermions_atomic_factors, &
   & updt_DD2_fermions_fermions_atomic_factors, dg_fermions_fermions_atomic_factors, &
   & da_fermions_fermions_atomic_factors, read_fermions_fermions_atomic_factors_params, &
   & save_fermions_fermions_atomic_factors_params, applysym_fermions_fermions_atomic_factors_params, &
   & updt_fermions_fermions_atomic_factors_params
   private :: comp_fermions_fermions_atomic_factors_symtable
contains
   subroutine init_fermions_fermions_atomic_factors( )
      g_ee_u_prs = .false.
      g_ee_d_prs = .false.
      g_ee_a_prs = .false.
      g_pe_u_prs = .false.
      g_pe_d_prs = .false.
      g_pp_prs = .false.
      if ( n_el.ge.2_int32 ) then
         if (spin_j.eq.'U') then
            if (n_el_u.ge.2_int32) then
               g_ee_u_prs = .true.
               allocate( G_ee_u(1:n_jorbs,1:n_jorbs) ) ; G_ee_u = 0.0_dp
            endif
            if (n_el_d.gt.0_int32) then
               g_ee_a_prs = .true.
               allocate( G_ee_a(1:n_jorbs,1:n_jorbs) ) ; G_ee_a = 0.0_dp
               if (n_el_d.ge.2_int32) then
                  g_ee_d_prs = .true.
                  allocate( G_ee_d(1:n_jorbs,1:n_jorbs) ) ; G_ee_d = 0.0_dp
               endif
            endif
         else
            g_ee_u_prs = .true.
            allocate( G_ee_u(1:n_jorbs,1:n_jorbs) ) ; G_ee_u = 0.0_dp
         endif
      else
         ejd2_opt = .false.
      endif
      if (n_po.gt.0_int32) then
         g_pe_u_prs = .true.
         allocate( G_pe_u(1:n_jorbs,1:n_jorbs) ) ; G_pe_u = 0.0_dp
         if (spin_j.eq.'U'.and.n_el_d.gt.0_int32) then
            g_pe_d_prs = .true.
            allocate( G_pe_d(1:n_jorbs,1:n_jorbs) ) ; G_pe_d = 0.0_dp
         endif
         if (n_po.ge.2_int32) then
            g_pp_prs = .true.
            allocate( G_pp(1:n_jorbs,1:n_jorbs) ) ; G_pp = 0.0_dp
         endif
      else
         pjd2_opt = .false.
      endif
      if ( ejd2_opt ) then
         if ( g_ee_u_prs ) then
            n_par_dyn2b   = n_par_dyn2b  + n_nz_G
            n_par_dyn2b_s = n_par_dyn2b_s + n_sym_G_ee
         endif
         if ( g_ee_a_prs ) then
            n_par_dyn2b   = n_par_dyn2b   + n_nz_G
            n_par_dyn2b_s = n_par_dyn2b_s + n_sym_G_pe
         endif
         if ( g_ee_d_prs ) then
            n_par_dyn2b   = n_par_dyn2b  + n_nz_G
            n_par_dyn2b_s = n_par_dyn2b_s + n_sym_G_ee
         endif
      endif  
      if ( pjd2_opt ) then
         if ( g_pe_u_prs ) then
            n_par_dyn2b   = n_par_dyn2b   + n_nz_G
            n_par_dyn2b_s = n_par_dyn2b_s + n_sym_G_pe
         endif
         if ( g_pe_d_prs ) then
            n_par_dyn2b   = n_par_dyn2b   + n_nz_G
            n_par_dyn2b_s = n_par_dyn2b_s + n_sym_G_pe
         endif
         if ( g_pp_prs  ) then
            n_par_dyn2b   = n_par_dyn2b  + n_nz_G
            n_par_dyn2b_s = n_par_dyn2b_s + n_sym_G_ee
         endif
      endif 
      n_par_jstd   = n_par_jstd + n_par_dyn2b
      n_par_jstd_s = n_par_jstd_s + n_par_dyn2b_s
      if (.not.allocated(vec_oeu_s)) then
         allocate( vec_oeu_s(1:n_jorbs,1:n_wlk_max) ) ; vec_oeu_s = 0.0_dp
      endif
      if (spin_j.eq.'U') then
         if (.not.allocated(vec_oed_s)) then
            allocate( vec_oed_s(1:n_jorbs,1:n_wlk_max) ) ; vec_oed_s = 0.0_dp
         endif
      endif
      if ( n_po.gt.0_int32 ) then
         if (.not.allocated(vec_op_s)) then
            allocate( vec_op_s(1:n_jorbs,1:n_wlk_max) )  ; vec_op_s = 0.0_dp
         endif
      endif
      allocate( mat_GX(1:n_jorbs,1:n_fe,1:n_wlk_max) ) ; mat_GX = 0.0_dp
      allocate( vec_GX_s(1:n_jorbs,1:n_wlk_max) )      ; vec_GX_s = 0.0_dp
      allocate( vec_GX_d(1:n_jorbs,1:n_wlk_max) )    ; vec_GX_d = 0.0_dp
      if ( ejd2_opt .or. pjd2_opt )then
         allocate( dg_2b(1:n_jorbs,1:n_jorbs,1:n_wlk_max) ) ; dg_2b = 0.0_dp
      endif
      if( jc_opt ) then
         if (.not.allocated(vec_da_oeu_s)) then
            allocate( vec_da_oeu_s(1:n_par_jorbs,1:n_wlk_max) ) ; vec_da_oeu_s = 0.0_dp
         endif
         if (spin_j.eq.'U') then
            if (.not.allocated(vec_da_oed_s)) then
               allocate( vec_da_oed_s(1:n_par_jorbs,1:n_wlk_max) ) ; vec_da_oed_s = 0.0_dp
            endif
         endif
         if(n_po.gt.0_int32) then
            if (.not.allocated(vec_da_op_s)) then
               allocate( vec_da_op_s(1:n_par_jorbs,1:n_wlk_max) ) ; vec_da_op_s = 0.0_dp
            endif
         endif
      endif
   end subroutine init_fermions_fermions_atomic_factors
   subroutine save_fermions_fermions_atomic_factors_params_st( head_fle )
      character(4), intent(in) :: head_fle
      character(100) :: svf_fle
      integer(int32), pointer :: o1
      integer(int32)          :: i1, i2
      svf_fle = 'wvfn.save/'//trim(head_fle)//'_s.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      if ( n_el.ge.2_int32 .or. n_po.ge.2_int32 ) then
      endif
      if ( n_el.ge.1_int32 .and. n_po.ge.1_int32 ) then
      endif
      close(1)
   end subroutine save_fermions_fermions_atomic_factors_params_st
   subroutine comp_fermions_fermions_atomic_factors( iw, jst_2b )
      integer(int32), intent(in)    :: iw
      real(dp),       intent(inout) :: jst_2b
      integer(int32) :: i1
      if (spin_j.eq.'U') then
         vec_GX_s(:,iw) = 0.0_dp
         if (n_el_u.ge.2_int32) then
            call dsymv( 'L', n_jorbs, 1.0_dp, G_ee_u, n_jorbs, vec_oeu_s(:,iw), &
            & 1_int32, 0.0_dp, vec_GX_s(:,iw), 1_int32 )
         endif
         if (n_el_d.gt.0_int32) then
            call dgemv( 'T', n_jorbs, n_jorbs, 1.0_dp, G_ee_a, n_jorbs, &
            & vec_oed_s(:,iw), 1_int32, 1.0_dp, vec_GX_s(:,iw), 1_int32 )
         endif
         if (n_po.gt.0_int32) then
            call dgemv( 'T', n_jorbs, n_jorbs, 1.0_dp, G_pe_u, n_jorbs, &
            & vec_op_s(:,iw), 1_int32, 1.0_dp, vec_GX_s(:,iw), 1_int32 )
         endif
         do i1 = 1_int32, n_el_u
            mat_GX(:,i1,iw) = vec_GX_s(:,iw)
         enddo
         if (n_el_d.gt.0_int32) then
            call dgemv( 'N', n_jorbs, n_jorbs, 1.0_dp, G_ee_a, n_jorbs, &
            & vec_oeu_s(:,iw), 1_int32, 0.0_dp, vec_GX_s(:,iw), 1_int32 )
            if (n_el_d.ge.2_int32) then
               call dsymv( 'L', n_jorbs, 1.0_dp, G_ee_d, n_jorbs, vec_oed_s(:,iw), &
               & 1_int32, 1.0_dp, vec_GX_s(:,iw), 1_int32 )
            endif
            if (n_po.gt.0_int32) then
               call dgemv( 'T', n_jorbs, n_jorbs, 1.0_dp, G_pe_d, n_jorbs, &
               & vec_op_s(:,iw), 1_int32, 1.0_dp, vec_GX_s(:,iw), 1_int32 )
            endif
            do i1 = n_el_u+1_int32, n_el
               mat_GX(:,i1,iw) = vec_GX_s(:,iw)
            enddo
         endif
         if (n_po.gt.0_int32) then
            call dgemv( 'N', n_jorbs, n_jorbs, 1.0_dp, G_pe_u, n_jorbs, &
            & vec_oeu_s(:,iw), 1_int32, 0.0_dp, vec_GX_s(:,iw), 1_int32 )
            if (n_el_d.gt.0_int32) then
               call dgemv( 'N', n_jorbs, n_jorbs, 1.0_dp, G_pe_d, n_jorbs, &
               & vec_oed_s(:,iw), 1_int32, 1.0_dp, vec_GX_s(:,iw), 1_int32 )
            endif
            if (n_po.gt.1_int32) then
               call dsymv( 'L', n_jorbs, 1.0_dp, G_pp, n_jorbs, vec_op_s(:,iw), &
               & 1_int32, 1.0_dp, vec_GX_s(:,iw), 1_int32 )
            endif
         endif
         do i1 = n_el+1_int32, n_fe
            mat_GX(:,i1,iw) = vec_GX_s(:,iw)
         enddo
         if (n_el_u.ge.2_int32) then
            call dsymm( 'L', 'L', n_jorbs, n_el_u, -1.0_dp, G_ee_u, n_jorbs, jorbs_mat(iw)%o(:,1:n_el_u), &
            & n_jorbs, 1.0_dp, mat_GX(:,1:n_el_u,iw), n_jorbs )
         endif
         if (n_el_d.ge.2_int32) then
            call dsymm( 'L', 'L', n_jorbs, n_el_d, -1.0_dp, G_ee_d, n_jorbs, jorbs_mat(iw)%o(:,n_el_u+1:n_el), &
            & n_jorbs, 1.0_dp, mat_GX(:,n_el_u+1:n_el,iw), n_jorbs )
         endif
         if (n_po.ge.2_int32) then
            call dsymm( 'L', 'L', n_jorbs, n_po, -1.0_dp, G_pp, n_jorbs, jorbs_mat(iw)%o(:,n_el+1:n_fe), &
            & n_jorbs, 1.0_dp, mat_GX(:,n_el+1:n_fe,iw), n_jorbs )
         endif
      else
         call dsymv( 'L', n_jorbs, 1.0_dp, G_ee_u, n_jorbs, vec_oeu_s(:,iw), &
         & 1_int32, 0.0_dp, vec_GX_s(:,iw), 1_int32 )
         if (n_po.gt.0_int32) then
            call dgemv( 'T', n_jorbs, n_jorbs, 1.0_dp, G_pe_u, n_jorbs, &
            & vec_op_s(:,iw), 1_int32, 1.0_dp, vec_GX_s(:,iw), 1_int32 )
         endif
         do i1 = 1_int32, n_el
            mat_GX(:,i1,iw) = vec_GX_s(:,iw)
         enddo
         if (n_po.gt.0_int32) then
            call dgemv( 'N', n_jorbs, n_jorbs, 1.0_dp, G_pe_u, n_jorbs, &
            & vec_oeu_s(:,iw), 1_int32, 0.0_dp, vec_GX_s(:,iw), 1_int32 )
            if (n_po.ge.2_int32) then
               call dsymv( 'L', n_jorbs, 1.0_dp, G_pp, n_jorbs, vec_op_s(:,iw), &
               & 1_int32, 1.0_dp, vec_GX_s(:,iw), 1_int32 )
            endif
            do i1 = n_el+1_int32, n_fe
               mat_GX(:,i1,iw) = vec_GX_s(:,iw)
            enddo
         endif
         call dsymm( 'L', 'L', n_jorbs, n_el, -1.0_dp, G_ee_u, n_jorbs, jorbs_mat(iw)%o(:,1:n_el), &
         & n_jorbs, 1.0_dp, mat_GX(:,1:n_el,iw), n_jorbs )
         if (n_po.ge.2_int32) then
            call dsymm( 'L', 'L', n_jorbs, n_po, -1.0_dp, G_pp, n_jorbs, jorbs_mat(iw)%o(:,n_el+1:n_fe), &
            & n_jorbs, 1.0_dp, mat_GX(:,n_el+1:n_fe,iw), n_jorbs )
         endif
      endif
      do i1 = 1_int32, n_fe
         jst_2b = jst_2b + 0.5_dp * dot_product( jorbs_mat(iw)%o(:,i1), mat_GX(:,i1,iw) )
      enddo
   end subroutine comp_fermions_fermions_atomic_factors
   subroutine variation_fermions_fermions_atomic_factors( iw, var_2b )
      integer(int32), intent(in)    :: iw
      real(dp),       intent(inout) :: var_2b
      var_2b = var_2b + dot_product( jorbs_mat(iw)%o_dif(:), mat_GX(:,frm_cnf(iw)%i_fe,iw) )
   end subroutine variation_fermions_fermions_atomic_factors
   subroutine updt_fermions_fermions_atomic_factors( iw )
      integer(int32), intent(in)    :: iw
      integer(int32) :: i1
      if (spin_j.eq.'U') then
         if ( frm_cnf(iw)%i_fe.le.n_el) then
            if ( frm_cnf(iw)%i_fe.le.n_el_u) then
               if ( n_el_u.ge.2_int32 ) then
                  call dsymv( 'L', n_jorbs, 1.0_dp, G_ee_u, n_jorbs, jorbs_mat(iw)%o_dif(:), &
                  & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
                  do i1 = 1_int32, n_el_u
                     if (i1.ne. frm_cnf(iw)%i_fe) mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
                  enddo
               endif
               if ( n_el_d.ge.1_int32 ) then
                  call dgemv( 'N', n_jorbs, n_jorbs, 1.0_dp, G_ee_a, n_jorbs, jorbs_mat(iw)%o_dif(:), &
                  & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
                  do i1 = n_el_u + 1_int32, n_el
                     mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
                  enddo
               endif
               if ( n_po.gt.0_int32 ) then
                  call dgemv( 'N', n_jorbs, n_jorbs, 1.0_dp, G_pe_u, n_jorbs, jorbs_mat(iw)%o_dif(:), &
                  & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
                  do i1 = n_el + 1_int32, n_fe
                     mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
                  enddo
               endif
            else
               if ( n_el_u.ge.1_int32 ) then
                  call dgemv( 'T', n_jorbs, n_jorbs, 1.0_dp, G_ee_a, n_jorbs, jorbs_mat(iw)%o_dif(:), &
                  & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
                  do i1 = 1_int32, n_el_u
                     mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
                  enddo
               endif
               if ( n_el_d.ge.2_int32 ) then
                  call dsymv( 'L', n_jorbs, 1.0_dp, G_ee_d, n_jorbs, jorbs_mat(iw)%o_dif(:), &
                  & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
                  do i1 = n_el_u + 1_int32, n_el
                     if (i1.ne.frm_cnf(iw)%i_fe) mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
                  enddo
               endif
               if ( n_po.gt.0_int32 ) then
                  call dgemv( 'N', n_jorbs, n_jorbs, 1.0_dp, G_pe_d, n_jorbs, jorbs_mat(iw)%o_dif(:), &
                  & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
                  do i1 = n_el + 1_int32, n_fe
                     mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
                  enddo
               endif
            endif
         else
            if ( n_el_u.gt.0_int32 ) then
               call dgemv( 'T', n_jorbs, n_jorbs, 1.0_dp, G_pe_u, n_jorbs, jorbs_mat(iw)%o_dif(:), &
               & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
               do i1 = 1_int32, n_el_u
                  mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
               enddo
            endif
            if ( n_el_d.gt.0_int32 ) then
               call dgemv( 'T', n_jorbs, n_jorbs, 1.0_dp, G_pe_u, n_jorbs, jorbs_mat(iw)%o_dif(:), &
               & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
               do i1 = n_el_u + 1_int32, n_el
                  mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
               enddo
            endif
            if ( n_po.gt.0_int32 ) then
               call dsymv( 'L', n_jorbs, 1.0_dp, G_pp, n_jorbs, jorbs_mat(iw)%o_dif(:), &
               & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
               do i1 = n_el + 1_int32, n_fe
                  if (i1.ne.frm_cnf(iw)%i_fe) mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
               enddo
            endif
         endif
      else
         if ( frm_cnf(iw)%i_fe.le.n_el) then
            if ( n_el.ge.2_int32 ) then
               call dsymv( 'L', n_jorbs, 1.0_dp, G_ee_u, n_jorbs, jorbs_mat(iw)%o_dif(:), &
               & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
               do i1 = 1_int32, n_el
                  if (i1.ne. frm_cnf(iw)%i_fe) mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
               enddo
            endif
            if (n_po.gt.0_int32) then
               call dgemv( 'N', n_jorbs, n_jorbs, 1.0_dp, G_pe_u, n_jorbs, jorbs_mat(iw)%o_dif(:), &
               & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
               do i1 = n_el + 1_int32, n_fe
                  mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
               enddo
            endif
         else
            call dgemv( 'T', n_jorbs, n_jorbs, 1.0_dp, G_pe_u, n_jorbs, jorbs_mat(iw)%o_dif(:), &
            & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
            do i1 = 1_int32, n_el
               mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
            enddo
            if ( n_po.ge.2_int32 ) then
               call dsymv( 'L', n_jorbs, 1.0_dp, G_pp, n_jorbs, jorbs_mat(iw)%o_dif(:), &
               & 1_int32, 0.0_dp, vec_GX_d(:,iw), 1_int32 )
               do i1 = n_el + 1_int32, n_fe
                  if (i1.ne. frm_cnf(iw)%i_fe) mat_GX(:,i1,iw) = mat_GX(:,i1,iw) + vec_GX_d(:,iw)
               enddo
            endif
         endif
      endif
   end subroutine updt_fermions_fermions_atomic_factors
   subroutine DD2_fermions_fermions_atomic_factors( iw, DD2_jst )
      integer(int32),              intent(in)    :: iw
      real(dp), dimension(4*n_fe), intent(inout) :: DD2_jst
      integer(int32) :: i1
      do i1 = 1_int32, n_fe
         call dgemv( 'T', n_jorbs, 4_int32, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*i1-3:4*i1), n_jorbs, &
         & mat_GX(:,i1,iw), 1_int32, 1.0_dp, DD2_jst(4*i1-3:4*i1), 1_int32)
      enddo
   end subroutine DD2_fermions_fermions_atomic_factors
   subroutine new_DD2_fermions_fermions_atomic_factors( iw, DD2_jst )
      integer(int32),         intent(in)    :: iw
      real(dp), dimension(4), intent(inout) :: DD2_jst
      call dgemv( 'T', n_jorbs, 4, 1.0_dp, jorbs_mat(iw)%DD2_o_new(:,1:4), n_jorbs,&
      &  mat_GX(:,frm_cnf(iw)%i_fe,iw), 1_int32, 1.0_dp, DD2_jst(:), 1_int32)
   end subroutine new_DD2_fermions_fermions_atomic_factors
   subroutine updt_DD2_fermions_fermions_atomic_factors( iw, DD2_jst  )
      integer(int32),                intent(in)    :: iw
      real(dp), dimension(4*n_fe), intent(inout) :: DD2_jst
      integer(int32) :: i1
      if ( frm_cnf(iw)%i_fe.ne.1_int32 ) then
         do i1 = 1_int32, frm_cnf(iw)%i_fe-1
            call dgemv( 'T', n_jorbs, 4_int32, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*i1-3:4*i1), n_jorbs, &
            & mat_GX(:,i1,iw), 1_int32, 1.0_dp, DD2_jst(4*i1-3:4*i1), 1_int32)
         enddo
      endif
      if ( frm_cnf(iw)%i_fe.ne.n_fe ) then
         do i1 = frm_cnf(iw)%i_fe+1, n_fe
            call dgemv( 'T', n_jorbs, 4_int32, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*i1-3:4*i1), n_jorbs, &
            & mat_GX(:,i1,iw), 1_int32, 1.0_dp, DD2_jst(4*i1-3:4*i1), 1_int32)
         enddo
      endif
   end subroutine updt_DD2_fermions_fermions_atomic_factors
   subroutine dg_fermions_fermions_atomic_factors( iw, vec_dg_2b )
      integer(int32),           intent(in)    :: iw
      real(dp), dimension(n_par_dyn2b),  intent(inout) :: vec_dg_2b
      integer(int32) :: ip
      ip = 1_int32
      if( ejd2_opt ) then
         if (g_ee_u_prs) then
            if (spin_j.eq.'U') then
               call cmp_dl_jffd( n_el_u, jorbs_mat(iw)%o(:,1:n_el_u), vec_oeu_s(:,iw), dg_2b(:,:,iw), vec_dg_2b(ip:ip+n_nz_G-1) )
            else
               call cmp_dl_jffd( n_el, jorbs_mat(iw)%o(:,1:n_el), vec_oeu_s(:,iw), dg_2b(:,:,iw), vec_dg_2b(ip:ip+n_nz_G-1) )
            endif
            ip = ip + n_nz_G
         endif
         if (g_ee_a_prs) then
            call cmp_dl_jped( vec_oed_s(:,iw), vec_oeu_s(:,iw), dg_2b(:,:,iw), vec_dg_2b(ip:ip+n_nz_G-1) )
            ip = ip + n_nz_G
         endif
         if (g_ee_d_prs) then
            call cmp_dl_jffd( n_el_d, jorbs_mat(iw)%o(:,n_el_u+1:n_el), vec_oed_s(:,iw), dg_2b(:,:,iw), vec_dg_2b(ip:ip+n_nz_G-1) )
            ip = ip + n_nz_G
         endif
      endif
      if ( pjd2_opt ) then
         if( g_pe_u_prs) then
            call cmp_dl_jped( vec_op_s(:,iw), vec_oeu_s(:,iw), dg_2b(:,:,iw), vec_dg_2b(ip:ip+n_nz_G-1) )
            ip = ip + n_nz_G
         endif
         if ( g_pe_d_prs) then
            call cmp_dl_jped( vec_op_s(:,iw), vec_oed_s(:,iw), dg_2b(:,:,iw), vec_dg_2b(ip:ip+n_nz_G-1) )
            ip = ip + n_nz_G
         endif
         if ( g_pp_prs ) then
            call cmp_dl_jffd( n_po, jorbs_mat(iw)%o(:,n_el+1:n_fe), vec_op_s(:,iw), dg_2b(:,:,iw), vec_dg_2b(ip:ip+n_nz_G-1) )
            ip = ip + n_nz_G
         endif
      endif
   contains
      subroutine cmp_dl_jffd( n_f, o, vec_of_s, mat_dg_2b, dg_2b )
         integer(int32),                            intent(in)    :: n_f
         real(dp), dimension(n_jorbs, 1:n_f),     intent(in)    :: o
         real(dp), dimension(n_jorbs),            intent(inout) :: vec_of_s
         real(dp), dimension(n_jorbs, 1:n_jorbs), intent(inout) :: mat_dg_2b
         real(dp), dimension(n_nz_G),             intent(inout) :: dg_2b
         integer(int32) :: i1
         mat_dg_2b = 0.0_dp
         do i1 = 1_int32, n_f
            vec_of_s = vec_of_s - o(:,i1)
            call dger( n_jorbs, n_jorbs, 1.0_dp, o(:,i1), 1_int32, &
            & vec_of_s, 1_int32, mat_dg_2b, n_jorbs)
            vec_of_s = vec_of_s + o(:,i1)
         enddo 
         do i1 = 1_int32, n_nz_G
            dg_2b(i1) = 0.5_dp * mat_dg_2b(nz_G(i1)%o1,nz_G(i1)%o2)
         enddo
      end subroutine cmp_dl_jffd
      subroutine cmp_dl_jped( vec_op_s, vec_oe_s, mat_dg_2b, dg_2b )
         real(dp), dimension(n_jorbs),            intent(inout) :: vec_op_s
         real(dp), dimension(n_jorbs),            intent(inout) :: vec_oe_s
         real(dp), dimension(n_jorbs, 1:n_jorbs), intent(inout) :: mat_dg_2b
         real(dp), dimension(n_nz_G),             intent(inout) :: dg_2b
         integer(int32) :: i1
         mat_dg_2b = 0.0_dp
         call dger( n_jorbs, n_jorbs, 1.0_dp, vec_op_s, 1_int32, &
         & vec_oe_s, 1_int32, mat_dg_2b, n_jorbs)
         do i1 = 1_int32, n_nz_G
            dg_2b(i1) = mat_dg_2b(nz_G(i1)%o1,nz_G(i1)%o2)
         enddo
      end subroutine cmp_dl_jped
   end subroutine dg_fermions_fermions_atomic_factors
   subroutine da_fermions_fermions_atomic_factors( iw, da_jst )
      integer(int32),           intent(in)    :: iw
      real(dp), dimension(n_par_jorbs),               intent(inout) :: da_jst
      integer(int32) :: i1, i2, i3, io, ip, fp
      do i1 = 1_int32, n_fe
         do i2 = 1_int32, n_at ; do i3 = 1_int32, jorbs(i2)%n_orbs
               io = jorbs(i2)%orb(i3)%i_orb
               if( jorbs(i2)%orb(i3)%n_par.gt.0_int32 ) then
                  ip = jorbs(i2)%orb(i3)%i_par
                  fp = ip + jorbs(i2)%orb(i3)%n_par - 1_int32
                  da_jst(ip:fp) = da_jst(ip:fp) + mat_GX(io,i1,iw) *  jorbs_mat(iw)%da_o(ip:fp,i1)
               endif
            enddo; enddo
      enddo 
   end subroutine da_fermions_fermions_atomic_factors
   subroutine read_fermions_fermions_atomic_factors_params( f_indx )
      integer(int32), intent(in) :: f_indx
      integer(int32) :: n_ele_2b
      integer(int32)             :: n_comp
      logical                    :: g_ee_u_p,  g_ee_a_p, g_ee_d_p, g_pe_u_p, g_pe_d_p, g_pp_p
      real(dp), allocatable, dimension(:) :: tmp_read
      integer(int32), allocatable, dimension(:), target :: indx_read
      real(dp) :: sum_par
      integer(int32), pointer :: i1dum, i2dum
      integer(int32) :: i1, i2
      if ( mpi_rank.eq.0_int32 ) then
         read(f_indx,*) ! Empty line
         read(f_indx,*) n_ele_2b, g_ee_u_p,  g_ee_a_p, g_ee_d_p, g_pe_u_p, g_pe_d_p, g_pp_p
         if ( n_ele_2b.gt.0_int32 ) then
            n_comp = 0_int32
            if ( g_ee_u_p ) n_comp = n_comp + 1_int32
            if ( g_ee_a_p ) n_comp = n_comp + 1_int32
            if ( g_ee_d_p ) n_comp = n_comp + 1_int32
            if ( g_pe_u_p ) n_comp = n_comp + 1_int32
            if ( g_pe_d_p ) n_comp = n_comp + 1_int32
            if (   g_pp_p ) n_comp = n_comp + 1_int32
            allocate( tmp_read (1:n_comp)) ; tmp_read = 0.0_dp
            allocate( indx_read(1:2)) ; indx_read = 0.0_dp
            do i1 = 1_int32, n_ele_2b
               read(f_indx,*) indx_read(1:2), tmp_read (1:n_comp)
               i1dum => indx_read(1);  i2dum => indx_read(2)
               i2 = 0_int32
               if ( g_ee_u_p ) then
                  i2 = i2 + 1_int32
                  g_ee_u(i1dum,i2dum) = tmp_read(i2)
                  if (spin_j.eq.'U') then
                     if ( .not.g_ee_a_p.and.g_ee_a_prs) then
                        g_ee_a(i1dum,i2dum) = tmp_read(i2)
                     endif
                     if ( .not.g_ee_d_p.and.g_ee_d_prs) then
                        g_ee_d(i1dum,i2dum) = tmp_read(i2)
                     endif
                  endif
               endif
               if ( g_ee_a_p ) then
                  i2 = i2 + 1_int32
                  g_ee_a(i1dum,i2dum) = tmp_read(i2)
               endif
               if ( g_ee_d_p ) then
                  i2 = i2 + 1_int32
                  g_ee_d(i1dum,i2dum) = tmp_read(i2)
               endif
               if ( g_pe_u_p ) then
                  i2 = i2 + 1_int32
                  g_pe_u(i1dum,i2dum) = tmp_read (i2)
                  if (spin_j.eq.'U') then
                     if ( .not.g_pe_d_p.and.g_pe_d_prs) then
                        g_pe_d(i1dum,i2dum) = tmp_read (i2)
                     endif
                  endif
               endif
               if ( g_pe_d_p ) then
                  i2 = i2 + 1_int32
                  g_pe_d(i1dum,i2dum) = tmp_read (i2)
               endif
               if ( g_pp_p ) then
                  i2 = i2 + 1_int32
                  g_pp(i1dum,i2dum) = tmp_read (i2)
               endif
            enddo
            deallocate( tmp_read, indx_read)
         endif
      endif !
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast( n_ele_2b, 1_int32, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr )
      if (n_ele_2b.gt.0_int32) then
         if ( g_ee_u_prs) call mpi_bcast( G_ee_u, n_jorbs**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,mpierr)
         if ( g_ee_a_prs) call mpi_bcast( G_ee_a, n_jorbs**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,mpierr)
         if ( g_ee_d_prs) call mpi_bcast( G_ee_d, n_jorbs**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,mpierr)
         if ( g_pe_u_prs) call mpi_bcast( G_pe_u, n_jorbs**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,mpierr)
         if ( g_pe_d_prs) call mpi_bcast( G_pe_d, n_jorbs**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,mpierr)
         if ( g_pp_prs)   call mpi_bcast( G_pp  , n_jorbs**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,mpierr)
      endif
#endif
      if (.not.jd2_opt) then
         sum_par = 0.0_dp
         if ( g_ee_u_prs) sum_par = sum_par + sum(abs(G_ee_u))
         if ( g_ee_a_prs) sum_par = sum_par + sum(abs(G_ee_a))
         if ( g_ee_d_prs) sum_par = sum_par + sum(abs(G_ee_d))
         if ( g_pe_u_prs) sum_par = sum_par + sum(abs(G_pe_u))
         if ( g_pe_d_prs) sum_par = sum_par + sum(abs(G_pe_d))
         if ( g_pp_prs)   sum_par = sum_par + sum(abs(G_pp))
         if ( sum_par.lt.10d-14 .and. jd2_prs ) then
            jd2_prs = .false.
            call write_variable_line(stdout,0,mpi_rank,2,"Two body coefficients smaller then 10d-14", jd2_prs,var_name="jd2_prs")
         endif
      endif
   end subroutine read_fermions_fermions_atomic_factors_params
   subroutine save_fermions_fermions_atomic_factors_params( f_indx )
      integer(int32), intent(in) :: f_indx
      integer(int32)             :: n_comp
      real(dp), allocatable, dimension(:) :: tmp_write
      integer(int32) :: i1, i2
      real(dp) :: sum_par
      n_comp = 0_int32
      if ( g_ee_u_prs) n_comp = n_comp + 1_int32
      if ( g_ee_a_prs) n_comp = n_comp + 1_int32
      if ( g_ee_d_prs) n_comp = n_comp + 1_int32
      if ( g_pe_u_prs) n_comp = n_comp + 1_int32
      if ( g_pe_d_prs) n_comp = n_comp + 1_int32
      if ( g_pp_prs)   n_comp = n_comp + 1_int32
      sum_par = 0.0_dp
      if ( g_ee_u_prs) sum_par = sum_par + sum(abs(G_ee_u))
      if ( g_ee_a_prs) sum_par = sum_par + sum(abs(G_ee_a))
      if ( g_ee_d_prs) sum_par = sum_par + sum(abs(G_ee_d))
      if ( g_pe_u_prs) sum_par = sum_par + sum(abs(G_pe_u))
      if ( g_pe_d_prs) sum_par = sum_par + sum(abs(G_pe_d))
      if ( g_pp_prs)   sum_par = sum_par + sum(abs(G_pp))
      if (sum_par.lt.10d-14) n_nz_G = 0_int32
      write(f_indx,'("# Coefficients of the 2Body linear Jastrow")')
      write(f_indx,'(I6,6L2)') n_nz_G, g_ee_u_prs, g_ee_a_prs, g_ee_d_prs, g_pe_u_prs, g_pe_d_prs, g_pp_prs
      if ( n_nz_G.gt.0_int32) then
         allocate(tmp_write(1:n_comp)) ; tmp_write = 0.0_dp
         do i1 = 1_int32, n_nz_G
            i2 = 0_int32
            tmp_write = 0.0_dp
            if ( g_ee_u_prs) then
               i2 = i2 + 1_int32
               tmp_write(i2) = G_ee_u(nz_G(i1)%o1, nz_G(i1)%o2)
            endif
            if ( g_ee_a_prs) then
               i2 = i2 + 1_int32
               tmp_write(i2)= G_ee_a(nz_G(i1)%o1, nz_G(i1)%o2)
            endif
            if ( g_ee_d_prs) then
               i2 = i2 + 1_int32
               tmp_write(i2) = G_ee_d(nz_G(i1)%o1, nz_G(i1)%o2)
            endif
            if ( g_pe_u_prs) then
               i2 = i2 + 1_int32
               tmp_write(i2) = G_pe_u(nz_G(i1)%o1, nz_G(i1)%o2)
            endif
            if ( g_pe_d_prs) then
               i2 = i2 + 1_int32
               tmp_write(i2)= G_pe_d(nz_G(i1)%o1, nz_G(i1)%o2)
            endif
            if ( g_pp_prs) then
               i2 = i2 + 1_int32
               tmp_write(i2)= G_pp(nz_G(i1)%o1, nz_G(i1)%o2)
            endif
            write(f_indx,'(2I6,6E18.9)') nz_G(i1)%o1, nz_G(i1)%o2, tmp_write(1:n_comp)
         enddo
         deallocate(tmp_write)
      endif
   end subroutine save_fermions_fermions_atomic_factors_params
   subroutine updt_fermions_fermions_atomic_factors_params( dp_s )
      real(dp), dimension(n_par_dyn2b_s), intent(in) :: dp_s
      integer(int32) :: i1, i2, ip
      integer(int32), pointer :: o1, o2
      ip = 0_int32
      if( ejd2_opt ) then
         if ( g_ee_u_prs ) then
            do i1 = 1_int32, n_sym_G_ee ; do i2 = 1_int32, sym_G_ee(i1)%n_c
                  o1 => nz_G(abs(sym_G_ee(i1)%c(i2)))%o1
                  o2 => nz_G(abs(sym_G_ee(i1)%c(i2)))%o2
                  G_ee_u(o1,o2) = G_ee_u(o1,o2) + dble(sign(1,sym_G_ee(i1)%c(i2))) * dp_s(ip+i1)
               enddo ; enddo
            ip = ip + n_sym_G_ee
         endif
         if( g_ee_a_prs ) then
            do i1 = 1_int32, n_sym_G_pe ; do i2 = 1_int32, sym_G_pe(i1)%n_c
                  o1 => nz_G(abs(sym_G_pe(i1)%c(i2)))%o1
                  o2 => nz_G(abs(sym_G_pe(i1)%c(i2)))%o2
                  G_ee_a(o1,o2) = G_ee_a(o1,o2) + dble(sign(1,sym_G_pe(i1)%c(i2))) * dp_s(ip+i1)
               enddo ;enddo
            ip = ip + n_sym_G_pe
         endif
         if ( g_ee_d_prs ) then
            do i1 = 1_int32, n_sym_G_ee ; do i2 = 1_int32, sym_G_ee(i1)%n_c
                  o1 => nz_G(abs(sym_G_ee(i1)%c(i2)))%o1
                  o2 => nz_G(abs(sym_G_ee(i1)%c(i2)))%o2
                  G_ee_d(o1,o2) = G_ee_d(o1,o2) + dble(sign(1,sym_G_ee(i1)%c(i2))) * dp_s(ip+i1)
               enddo ; enddo
            ip = ip + n_sym_G_ee
         endif
      endif 
      if( pjd2_opt ) then
         if ( g_pe_u_prs ) then
            do i1 = 1_int32, n_sym_G_pe ; do i2 = 1_int32, sym_G_pe(i1)%n_c
                  o1 => nz_G(abs(sym_G_pe(i1)%c(i2)))%o1
                  o2 => nz_G(abs(sym_G_pe(i1)%c(i2)))%o2
                  G_pe_u(o1,o2) = G_pe_u(o1,o2) + dble(sign(1,sym_G_pe(i1)%c(i2))) * dp_s(ip+i1)
               enddo ; enddo
            ip = ip + n_sym_G_pe
         endif
         if ( g_pe_d_prs ) then
            do i1 = 1_int32, n_sym_G_pe ; do i2 = 1_int32, sym_G_pe(i1)%n_c
                  o1 => nz_G(abs(sym_G_pe(i1)%c(i2)))%o1
                  o2 => nz_G(abs(sym_G_pe(i1)%c(i2)))%o2
                  G_pe_d(o1,o2) = G_pe_d(o1,o2) + dble(sign(1,sym_G_pe(i1)%c(i2))) * dp_s(ip+i1)
               enddo ; enddo
            ip = ip + n_sym_G_pe
         endif
         if ( g_pp_prs ) then
            do i1 = 1_int32, n_sym_G_ee ; do i2 = 1_int32, sym_G_ee(i1)%n_c
                  o1 => nz_G(abs(sym_G_ee(i1)%c(i2)))%o1
                  o2 => nz_G(abs(sym_G_ee(i1)%c(i2)))%o2
                  G_pp(o1,o2) = G_pp(o1,o2) + dble(sign(1,sym_G_ee(i1)%c(i2))) * dp_s(ip+i1)
               enddo ; enddo
            ip = ip + n_sym_G_ee
         endif 
      endif 
   end subroutine updt_fermions_fermions_atomic_factors_params
   subroutine applysym_fermions_fermions_atomic_factors_params( dp_ln_jst, dp_ln_jst_s )
      real(dp), dimension(n_par_dyn2b),   intent(in)    :: dp_ln_jst
      real(dp), dimension(n_par_dyn2b_s), intent(inout) :: dp_ln_jst_s
      integer(int32) :: i1, i2, i3, i4
      i3 = 0_int32 ; i4 = 0_int32
      if ( ejd2_opt ) then
         if ( g_ee_u_prs) then
            do i1 = 1_int32, n_sym_G_ee ; do i2 = 1_int32, sym_G_ee(i1)%n_c
                  dp_ln_jst_s(i3+i1) = dp_ln_jst_s(i3+i1) + dble(sign(1,sym_G_ee(i1)%c(i2))) * &
                  & dp_ln_jst(i4 + abs(sym_G_ee(i1)%c(i2)))
               enddo ; enddo
            i3 = i3 + n_sym_G_ee
            i4 = i4 + n_nz_G
         endif
         if ( g_ee_a_prs) then
            do i1 = 1_int32, n_sym_G_pe ; do i2 = 1_int32, sym_G_pe(i1)%n_c
                  dp_ln_jst_s(i3+i1) = dp_ln_jst_s(i3+i1) + dble(sign(1,sym_G_pe(i1)%c(i2))) * &
                  & dp_ln_jst(i4 + abs(sym_G_pe(i1)%c(i2)))
               enddo ; enddo
            i3 = i3 + n_sym_G_pe
            i4 = i4 + n_nz_G
         endif
         if ( g_ee_d_prs) then
            do i1 = 1_int32, n_sym_G_ee ; do i2 = 1_int32, sym_G_ee(i1)%n_c
                  dp_ln_jst_s(i3+i1) = dp_ln_jst_s(i3+i1) + dble(sign(1,sym_G_ee(i1)%c(i2))) * &
                  & dp_ln_jst(i4 + abs(sym_G_ee(i1)%c(i2)))
               enddo ; enddo
            i3 = i3 + n_sym_G_ee
            i4 = i4 + n_nz_G
         endif
      endif 
      if ( pjd2_opt ) then
         if ( g_pe_u_prs) then
            do i1 = 1_int32, n_sym_G_pe ; do i2 = 1_int32, sym_G_pe(i1)%n_c
                  dp_ln_jst_s(i3+i1) = dp_ln_jst_s(i3+i1) + dble(sign(1,sym_G_pe(i1)%c(i2))) * &
                  & dp_ln_jst(i4 + abs(sym_G_pe(i1)%c(i2)))
               enddo ; enddo
            i3 = i3 + n_sym_G_pe
            i4 = i4 + n_nz_G
         endif 
         if ( g_pe_d_prs) then
            do i1 = 1_int32, n_sym_G_pe ; do i2 = 1_int32, sym_G_pe(i1)%n_c
                  dp_ln_jst_s(i3+i1) = dp_ln_jst_s(i3+i1) + dble(sign(1,sym_G_pe(i1)%c(i2))) * &
                  & dp_ln_jst(i4 + abs(sym_G_pe(i1)%c(i2)))
               enddo ; enddo
            i3 = i3 + n_sym_G_pe
            i4 = i4 + n_nz_G
         endif 
         if ( g_pp_prs ) then
            do i1 = 1_int32, n_sym_G_ee ; do i2 = 1_int32, sym_G_ee(i1)%n_c
                  dp_ln_jst_s(i3+i1) = dp_ln_jst_s(i3+i1) + dble(sign(1,sym_G_ee(i1)%c(i2))) * &
                  & dp_ln_jst(i4 + abs(sym_G_ee(i1)%c(i2)))
               enddo ; enddo
            i3 = i3 + n_sym_G_ee
            i4 = i4 + n_nz_G
         endif
      endif 
   end subroutine applysym_fermions_fermions_atomic_factors_params
end module fermions_fermions_atomic_factors_m
