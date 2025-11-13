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
module fermionic_nuclei_atomic_factors_m
   use fortran_kinds_v,       only: dp, int32, stdout
   use openmp_mpi_m
   use write_lines_m
   use quantum_monte_carlo_v, only: n_wlk_max
   use orbitalsbased_jastrow_factors_params_m
   use fermionic_config_c,    only: frm_cnf
   use molecular_system_v,    only: n_at, n_el_u, n_el_d, n_po, n_el, n_fe
   use jastrow_factors_v,     only: jd1_prs, jd1_opt, ejd1_opt, pjd1_opt, n_par_jstd, n_par_jstd_s
   use orbbss_cls,            only: orbspa_t
   use jastrow_orbitals_m,    only: jorbs, jorbs_s, jorbs_mat, n_jorbs, n_par_jorbs, jc_opt, je_opt
   implicit none
   public :: init_fermionic_nuclei_atomic_factors, comp_fermionic_nuclei_atomic_factors, &
   & DD2_fermionic_nuclei_atomic_factors, dg_fermionic_nuclei_atomic_factors, da_fermionic_nuclei_atomic_factors, &
   & variation_fermionic_nuclei_atomic_factors, new_DD2_fermionic_nuclei_atomic_factors, &
   & comp_orbitals_sums, updt_orbitals_sums, da_orbitals_sums, save_fermionic_nuclei_atomic_factors_params, read_fermionic_nuclei_atomic_factors_params
contains
   subroutine init_fermionic_nuclei_atomic_factors( )
      if ( n_el.eq.0_int32 ) ejd1_opt = .false.
      if ( n_po.eq.0_int32 ) pjd1_opt = .false.
      if (n_jorbs.gt.0_int32) then
         allocate( g_en_u(1:n_jorbs) ) ; g_en_u = 0.0_dp
         if (spin_j.eq.'U') then
            allocate( g_en_d(1:n_jorbs) ) ; g_en_d = 0.0_dp
         endif
         if (n_po.gt.0_int32) then
            allocate( g_pn(1:n_jorbs) ) ; g_pn = 0.0_dp
         endif
         jd1_prs = .true.
      else
         jd1_prs = .false.
         ejd1_opt = .false.
         pjd1_opt = .false.
      endif
      if ( ejd1_opt ) then
         n_par_dyn1b   = n_par_dyn1b + n_jorbs
         n_par_dyn1b_s = n_par_dyn1b_s + jorbs_s%orbs_sym_tab%n_indp
         if (spin_j.eq.'U') then
            n_par_dyn1b   = n_par_dyn1b + n_jorbs
            n_par_dyn1b_s = n_par_dyn1b_s + jorbs_s%orbs_sym_tab%n_indp
         endif
      endif
      if (pjd1_opt ) then
         n_par_dyn1b   = n_par_dyn1b + n_jorbs
         n_par_dyn1b_s = n_par_dyn1b_s + jorbs_s%orbs_sym_tab%n_indp
      endif
      n_par_jstd   = n_par_jstd   + n_par_dyn1b
      n_par_jstd_s = n_par_jstd_s + n_par_dyn1b_s
#ifdef _PRNTALL
      if ( mpi_rank.eq.0_int32 .and. jd1_prs ) then
         call write_empty_line(stdout,0,mpi_rank)
         if ( allocated(g_en_u) ) &
         & write(*,'(3X,"Initialized Elec.-Nucl. dynamical Jastrow                      ")')
         if ( allocated(g_pn) ) &
         & write(*,'(3X,"Initialized Posi.-Nucl. dynamical Jastrow                      ")')
      endif
#endif
      allocate( vec_oeu_s(1:n_jorbs,1:n_wlk_max) ) ; vec_oeu_s = 0.0_dp
      if (spin_j.eq.'U') then
         allocate( vec_oed_s(1:n_jorbs,1:n_wlk_max) ) ; vec_oed_s = 0.0_dp
      endif
      if ( n_po.gt.0_int32 ) then
         allocate( vec_op_s(1:n_jorbs,1:n_wlk_max) )  ; vec_op_s = 0.0_dp
      endif
      if( jc_opt  ) then
         allocate( vec_da_oeu_s(1:n_par_jorbs,1:n_wlk_max) ) ; vec_da_oeu_s = 0.0_dp
         if (spin_j.eq.'U') then
            allocate( vec_da_oed_s(1:n_par_jorbs,1:n_wlk_max) ) ; vec_da_oed_s = 0.0_dp
         endif
         if(n_po.gt.0_int32) then
            allocate( vec_da_op_s(1:n_par_jorbs,1:n_wlk_max) ) ; vec_da_op_s = 0.0_dp
         endif
      endif
   end subroutine init_fermionic_nuclei_atomic_factors
   subroutine comp_orbitals_sums( iw )
      integer(int32), intent(in) :: iw
      integer(int32) :: i1
      if (spin_j.eq.'U') then
         vec_oeu_s(:,iw) = 0.0_dp
         do i1 = 1_int32, n_el_u
            vec_oeu_s(:,iw) = vec_oeu_s(:,iw) + jorbs_mat(iw)%o(:,i1)
         enddo
         vec_oed_s(:,iw) = 0.0_dp
         do i1 = n_el_u+1_int32, n_el
            vec_oed_s(:,iw) = vec_oed_s(:,iw) + jorbs_mat(iw)%o(:,i1)
         enddo
      else
         vec_oeu_s(:,iw) = 0.0_dp
         do i1 = 1_int32, n_el
            vec_oeu_s(:,iw) = vec_oeu_s(:,iw) + jorbs_mat(iw)%o(:,i1)
         enddo
      endif
      if (n_po.gt.0_int32) then
         vec_op_s(:,iw) = 0.0_dp
         do i1 = n_el+1_int32, n_fe
            vec_op_s(:,iw) = vec_op_s(:,iw) + jorbs_mat(iw)%o(:,i1)
         enddo
      endif
   end subroutine comp_orbitals_sums
   subroutine updt_orbitals_sums( iw )
      integer(int32), intent(in) :: iw
      if ( spin_j.eq.'U' ) then
         if (frm_cnf(iw)%i_fe.le.n_el) then
            if (frm_cnf(iw)%i_fe.le.n_el_u) then
               vec_oeu_s(:,iw) = vec_oeu_s(:,iw) + jorbs_mat(iw)%o_dif(:)
            else
               vec_oed_s(:,iw) = vec_oed_s(:,iw) + jorbs_mat(iw)%o_dif(:)
            endif
         else
            vec_op_s(:,iw) = vec_op_s(:,iw) + jorbs_mat(iw)%o_dif(:)
         endif
      else
         if (frm_cnf(iw)%i_fe.le.n_el) then
            vec_oeu_s(:,iw) = vec_oeu_s(:,iw) + jorbs_mat(iw)%o_dif(:)
         else
            vec_op_s(:,iw) = vec_op_s(:,iw) + jorbs_mat(iw)%o_dif(:)
         endif
      endif
   end subroutine updt_orbitals_sums
   subroutine da_orbitals_sums( iw )
      integer(int32), intent(in) :: iw
      integer(int32) :: i1
      if (spin_j.eq.'U') then
         vec_da_oeu_s(:,iw) = jorbs_mat(iw)%da_o(:,1)
         do i1 = 2_int32, n_el_u
            vec_da_oeu_s(:,iw) = vec_da_oeu_s(:,iw) + jorbs_mat(iw)%da_o(:,i1)
         enddo
         vec_da_oed_s(:,iw) = jorbs_mat(iw)%da_o(:,n_el_u+1)
         do i1 = n_el_u+2_int32, n_el
            vec_da_oed_s(:,iw) = vec_da_oed_s(:,iw) + jorbs_mat(iw)%da_o(:,i1)
         enddo
      else
         vec_da_oeu_s(:,iw) = jorbs_mat(iw)%da_o(:,1)
         do i1 = 2_int32, n_el
            vec_da_oeu_s(:,iw) = vec_da_oeu_s(:,iw) + jorbs_mat(iw)%da_o(:,i1)
         enddo
      endif
      if ( n_po.gt.0_int32 ) then
         vec_da_op_s(:,iw) = jorbs_mat(iw)%da_o(:,n_el+1)
         do i1 = n_el + 2_int32, n_fe
            vec_da_op_s(:,iw) = vec_da_op_s(:,iw) + jorbs_mat(iw)%da_o(:,i1)
         enddo
      endif
   end subroutine da_orbitals_sums
   subroutine comp_fermionic_nuclei_atomic_factors( iw, jst_1b )
      integer(int32), intent(in)    :: iw
      real(dp),       intent(inout) :: jst_1b
      if (spin_j.eq.'U') then
         jst_1b = jst_1b + dot_product( vec_oeu_s(:,iw) , g_en_u  ) + dot_product( vec_oed_s(:,iw) , g_en_d  )
      else
         jst_1b = jst_1b + dot_product( vec_oeu_s(:,iw) , g_en_u  )
      endif
      if (n_po.gt.0_int32) then
         jst_1b = jst_1b + dot_product( vec_op_s(:,iw) , g_pn )
      endif
   end subroutine comp_fermionic_nuclei_atomic_factors
   subroutine variation_fermionic_nuclei_atomic_factors( iw, var_1b )
      integer(int32), intent(in)    :: iw
      real(dp),       intent(inout) :: var_1b
      if (spin_j.eq.'U') then
         if (frm_cnf(iw)%i_fe.le.n_el) then
            if (frm_cnf(iw)%i_fe.le.n_el_u) then
               var_1b = var_1b + dot_product( jorbs_mat(iw)%o_dif(:) , g_en_u(:) )
            else
               var_1b = var_1b + dot_product( jorbs_mat(iw)%o_dif(:) , g_en_d(:) )
            endif
         else
            var_1b = var_1b + dot_product( jorbs_mat(iw)%o_dif(:) , g_pn(:) )
         endif
      else
         if (frm_cnf(iw)%i_fe.le.n_el) then
            var_1b = var_1b + dot_product( jorbs_mat(iw)%o_dif(:) , g_en_u(:) )
         else
            var_1b = var_1b + dot_product( jorbs_mat(iw)%o_dif(:) , g_pn(:) )
         endif
      endif
   end subroutine variation_fermionic_nuclei_atomic_factors
   subroutine DD2_fermionic_nuclei_atomic_factors( iw, DD2_1b )
      integer(int32),                intent(in)    :: iw
      real(dp), dimension(4*n_fe), intent(inout) :: DD2_1b
      if (spin_j.eq.'U') then
         call dgemv('T', n_jorbs, 4*n_el_u, 1.0_dp, jorbs_mat(iw)%DD2_o(:,1:4*n_el_u), n_jorbs, &
         & g_en_u, 1_int32, 1.0_dp, DD2_1b(1:4*n_el_u), 1_int32)
         call dgemv('T', n_jorbs, 4*n_el_d, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*n_el_u+1:4*n_el), n_jorbs, &
         & g_en_d, 1_int32, 1.0_dp, DD2_1b(4*n_el_u+1:4*n_el), 1_int32)
      else
         call dgemv('T', n_jorbs, 4*n_el, 1.0_dp, jorbs_mat(iw)%DD2_o(:,1:4*n_el), n_jorbs, &
         & g_en_u, 1_int32, 1.0_dp, DD2_1b(1:4*n_el), 1_int32)
      endif
      if ( n_po.gt.0_int32 ) then
         call dgemv('T', n_jorbs, 4*n_po, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*n_el+1:4*n_fe), n_jorbs, &
         & g_pn, 1_int32, 1.0_dp, DD2_1b(4*n_el+1:4*n_fe), 1_int32)
      endif
   end subroutine DD2_fermionic_nuclei_atomic_factors
   subroutine new_DD2_fermionic_nuclei_atomic_factors( iw, DD2_1b )
      integer(int32),           intent(in)    :: iw
      real(dp), dimension(4), intent(inout) :: DD2_1b
      if (spin_j.eq.'U') then
         if ( frm_cnf(iw)%i_fe.le.n_el ) then
            if (frm_cnf(iw)%i_fe.le.n_el_u) then
               call dgemv('T', n_jorbs, 4_int32, 1.0_dp, jorbs_mat(iw)%DD2_o_new(:,:), n_jorbs, &
               & g_en_u, 1_int32, 1.0_dp, DD2_1b, 1_int32)
            else
               call dgemv('T', n_jorbs, 4_int32, 1.0_dp, jorbs_mat(iw)%DD2_o_new(:,:), n_jorbs, &
               & g_en_d, 1_int32, 1.0_dp, DD2_1b, 1_int32)
            endif
         else
            call dgemv('T', n_jorbs, 4_int32, 1.0_dp, jorbs_mat(iw)%DD2_o_new(:,:), n_jorbs, &
            & g_pn, 1_int32, 1.0_dp, DD2_1b, 1_int32)
         endif
      else
         if ( frm_cnf(iw)%i_fe.le.n_el ) then
            call dgemv('T', n_jorbs, 4_int32, 1.0_dp, jorbs_mat(iw)%DD2_o_new(:,:), n_jorbs, &
            & g_en_u, 1_int32, 1.0_dp, DD2_1b, 1_int32)
         else
            call dgemv('T', n_jorbs, 4_int32, 1.0_dp, jorbs_mat(iw)%DD2_o_new(:,:), n_jorbs, &
            & g_pn, 1_int32, 1.0_dp, DD2_1b, 1_int32)
         endif
      endif
   end subroutine new_DD2_fermionic_nuclei_atomic_factors
   subroutine updt_DD2_fermionic_nuclei_atomic_factors( iw, DD2_1b )
      integer(int32),                intent(in)    :: iw
      real(dp), dimension(4*n_fe), intent(inout) :: DD2_1b
      integer(int32) :: i1
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         if (spin_j.eq.'U') then
            do i1 = 1_int32, n_el_u
               if (i1.ne.frm_cnf(iw)%i_fe) then
                  call dgemv('T', n_jorbs, 4, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*i1-3:4*i1), n_jorbs, &
                  & g_en_u, 1_int32, 1.0_dp, DD2_1b(4*i1-3:4*i1), 1_int32)
               endif
            enddo
            do i1 = n_el_u+1,n_el
               if (i1.ne.frm_cnf(iw)%i_fe) then
                  call dgemv('T', n_jorbs, 4, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*i1-3:4*i1), n_jorbs, &
                  & g_en_u, 1_int32, 1.0_dp, DD2_1b(4*i1-3:4*i1), 1_int32)
               endif
            enddo
         else
            do i1 = 1_int32, n_el
               if (i1.ne.frm_cnf(iw)%i_fe) then
                  call dgemv('T', n_jorbs, 4, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*i1-3:4*i1), n_jorbs, &
                  & g_en_u, 1_int32, 1.0_dp, DD2_1b(4*i1-3:4*i1), 1_int32)
               endif
            enddo
         endif
         if ( n_po.gt.0_int32 ) then
            call dgemv('T', n_jorbs, 4*n_po, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*n_el+1:4*n_fe), n_jorbs, &
            & g_pn, 1_int32, 1.0_dp, DD2_1b(4*n_el+1:4*n_fe), 1_int32)
         endif
      else
         if (spin_j.eq.'U') then
            call dgemv('T', n_jorbs, 4*n_el_u, 1.0_dp, jorbs_mat(iw)%DD2_o(:,1:4*n_el_u), n_jorbs, &
            & g_en_u, 1_int32, 1.0_dp, DD2_1b(1:4*n_el_u), 1_int32)
            call dgemv('T', n_jorbs, 4*n_el_d, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*n_el_u+1:4*n_el), n_jorbs, &
            & g_en_d, 1_int32, 1.0_dp, DD2_1b(4*n_el_u+1:4*n_el), 1_int32)
         else
            call dgemv('T', n_jorbs, 4*n_el, 1.0_dp, jorbs_mat(iw)%DD2_o(:,1:4*n_el), n_jorbs, &
            & g_en_u, 1_int32, 1.0_dp, DD2_1b(1:4*n_el), 1_int32)
         endif
         do i1 = n_el+1, n_fe
            if (i1.ne.frm_cnf(iw)%i_fe) then
               call dgemv('T', n_jorbs, 4, 1.0_dp, jorbs_mat(iw)%DD2_o(:,4*i1-3:4*i1), n_jorbs, &
               & g_pn, 1_int32, 1.0_dp, DD2_1b(4*i1-3:4*i1), 1_int32)
            endif
         enddo
      endif
   end subroutine updt_DD2_fermionic_nuclei_atomic_factors
   subroutine dg_fermionic_nuclei_atomic_factors( iw, dp_ln_jst )
      integer(int32),           intent(in)    :: iw
      real(dp), dimension(n_par_dyn1b),         intent(inout) :: dp_ln_jst
      integer(int32) :: ip
      ip = 1
      if ( ejd1_opt ) then
         dp_ln_jst(ip:ip+n_jorbs-1) = vec_oeu_s(:,iw)
         ip = ip + n_jorbs
         if (spin_j.eq.'U'.and.n_el_d.gt.0_int32) then
            dp_ln_jst(ip:ip+n_jorbs-1) = vec_oed_s(:,iw)
            ip = ip + n_jorbs
         endif
      endif
      if ( pjd1_opt ) then
         dp_ln_jst(ip:ip+n_jorbs-1) = vec_op_s(:,iw)
         ip = ip + n_jorbs
      endif
   end subroutine dg_fermionic_nuclei_atomic_factors
   subroutine da_fermionic_nuclei_atomic_factors( iw, da_jst )
      integer(int32),           intent(in)    :: iw
      real(dp), dimension(n_par_jorbs),           intent(inout) :: da_jst
      integer(int32) :: i1, i2, io, ip, fp
      do i1 = 1_int32, n_at ; do i2 = 1_int32, jorbs(i1)%n_orbs
            io = jorbs(i1)%orb(i2)%i_orb
            if( jorbs(i1)%orb(i2)%n_par.gt.0_int32 ) then
               ip = jorbs(i1)%orb(i2)%i_par
               fp = ip + jorbs(i1)%orb(i2)%n_par - 1_int32
               da_jst(ip:fp) = da_jst(ip:fp) + g_en_u(io) * vec_da_oeu_s(ip:fp,iw)
               if (spin_j.eq.'U') then
                  da_jst(ip:fp) = da_jst(ip:fp) + g_en_d(io) * vec_da_oed_s(ip:fp,iw)
               endif
               if ( n_po.gt.0_int32 ) then
                  da_jst(ip:fp) = da_jst(ip:fp) + g_pn(io) * vec_da_op_s(ip:fp,iw)
               endif
            endif
         enddo ; enddo 
   end subroutine da_fermionic_nuclei_atomic_factors
   subroutine read_fermionic_nuclei_atomic_factors_params( f_indx )
      integer(int32), intent(in) :: f_indx
      integer(int32)             :: n_ele_1b
      real(dp)                   :: sum_par
      integer(int32)             :: n_comp
      logical                    :: g_en_u_prs,  g_en_d_prs, g_pn_prs
      real(dp), allocatable, dimension(:,:) :: tmp_read
      integer(int32) :: i1
      if ( mpi_rank.eq.0_int32 ) then
         read(f_indx,*) ! Number of 1body elements
         read(f_indx,*) n_ele_1b, g_en_u_prs, g_en_d_prs, g_pn_prs
         n_comp = 0_int32
         if ( g_en_u_prs ) n_comp = n_comp + 1_int32
         if ( g_en_d_prs ) n_comp = n_comp + 1_int32
         if ( g_pn_prs ) n_comp = n_comp + 1_int32
         if ( n_ele_1b.gt.0_int32 ) then
            allocate( tmp_read(1:n_ele_1b,1:n_comp) ) ; tmp_read = 0.0_dp
            do i1 = 1_int32, n_ele_1b
               read(f_indx,*) tmp_read(i1,1:n_comp)
            enddo 
            n_comp = 0_int32
            if ( g_en_u_prs ) then
               n_comp = n_comp + 1_int32
               g_en_u = tmp_read(:,n_comp)
            endif
            if ( g_en_d_prs ) then
               n_comp = n_comp + 1_int32
               g_en_d = tmp_read(:,n_comp)
            endif
            if ( g_pn_prs ) then
               n_comp = n_comp + 1_int32
               g_pn = tmp_read(:,n_comp)
            endif
            deallocate( tmp_read )
         endif
      endif 
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast( n_ele_1b, 1_int32, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr )
      if ( n_ele_1b .ne. 0_int32 ) then
         if ( allocated(g_en_u) ) call mpi_bcast( g_en_u, n_jorbs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,mpierr)
         if ( allocated(g_en_d) ) call mpi_bcast( g_en_d, n_jorbs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,mpierr)
         if ( allocated(  g_pn) ) call mpi_bcast(   g_pn, n_jorbs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,mpierr)
      endif
#endif
      if (.not.jd1_opt ) then
         sum_par = 0.0_dp
         if ( allocated(g_en_u) ) sum_par = sum_par + sum(abs(g_en_u))
         if ( allocated(g_en_d) ) sum_par = sum_par + sum(abs(g_en_d))
         if ( allocated(g_pn) ) sum_par = sum_par + sum(abs(g_pn))
         if ( sum_par.lt.10d-14 .and. jd1_prs ) then
            jd1_prs = .false.
            call write_variable_line(stdout,0,mpi_rank,2,"One body coefficients smaller then 10d-14", jd1_prs,var_name="jd1_prs")
         endif
      endif
   end subroutine read_fermionic_nuclei_atomic_factors_params
   subroutine save_fermionic_nuclei_atomic_factors_params( f_indx )
      integer(int32), intent(in) :: f_indx
      integer(int32)             :: n_comp
      real(dp), allocatable, dimension(:,:) :: tmp_write
      integer(int32) :: i1
      n_comp = 0_int32
      if ( allocated(g_en_u) ) n_comp = n_comp + 1_int32
      if ( allocated(g_en_d) ) n_comp = n_comp + 1_int32
      if ( allocated(g_pn) ) n_comp = n_comp + 1_int32
      allocate( tmp_write(1:n_jorbs,1:n_comp) ) ; tmp_write = 0.0_dp
      n_comp = 0_int32
      if ( allocated(g_en_u) ) then
         n_comp = n_comp + 1_int32
         tmp_write(:,n_comp) = g_en_u
      endif
      if ( allocated(g_en_d) ) then
         n_comp = n_comp + 1_int32
         tmp_write(:,n_comp) = g_en_d
      endif
      if ( allocated(g_pn) ) then
         n_comp = n_comp + 1_int32
         tmp_write(:,n_comp) = g_pn
      endif
      write(f_indx,'("# Coefficients of the 1Body linear Jastrow ")')
      if (sum(abs(tmp_write)).gt.0.0_dp) then
         write(f_indx,'(I6,3L2)') n_jorbs, allocated(g_en_u), allocated(g_en_d), allocated(g_pn)
         do i1 = 1_int32, n_jorbs
            write(f_indx,'(3E18.9)') tmp_write(i1,1:n_comp)
         enddo
      else
         write(f_indx,'(I6,3L2)') 0_int32, allocated(g_en_u), allocated(g_en_d), allocated(g_pn)
      endif
      deallocate(tmp_write)
   end subroutine save_fermionic_nuclei_atomic_factors_params
end module fermionic_nuclei_atomic_factors_m
