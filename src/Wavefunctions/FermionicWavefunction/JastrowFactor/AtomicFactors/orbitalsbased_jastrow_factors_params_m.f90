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
module orbitalsbased_jastrow_factors_params_m
   use fortran_kinds_v,      only: int32, dp
   use geminal_symmetries_m, only: cmp_nzr_fee, cmp_sgm_sym, lmbd_ele_t, lmbd_sym_t
   use molecular_system_v,   only: n_el, n_at, n_fe, n_po, n_sys_sym
   use jastrow_orbitals_m,   only: jorbs, n_jorbs
   implicit none
   character(1), save, public     :: spin_j, cpl_j, cpl_m
   integer(int32),   save, public :: n_par_dyn1b, n_par_dyn1b_s
   logical, save, public :: g_ee_u_prs, g_ee_d_prs, g_ee_a_prs, g_pe_u_prs, g_pe_d_prs, g_pp_prs
   integer(int32),   save, public :: corr_type  
   integer(int32),   save, public :: order_coupling  
   integer(int32),   save, public :: n_par_dyn2b, n_par_dyn2b_s
   real(dp), save, public, allocatable, dimension(:,:) :: G_ee_u, G_ee_d, G_ee_a
   real(dp), save, public, allocatable, dimension(:,:) :: G_pp
   real(dp), save, public, allocatable, dimension(:,:) :: G_pe_u, G_pe_d
   integer(int32),   save, public :: n_nz_G
   type(lmbd_ele_t), save, public, allocatable, dimension(:), target :: nz_G
   integer(int32),   save, public :: n_sym_G_ee, n_sym_G_pe
   type(lmbd_sym_t), save, public, allocatable, dimension(:), target :: sym_G_ee, sym_G_pe
   real(dp), save, public, allocatable, dimension(:), target :: g_en_u, g_en_d
   real(dp), save, public, allocatable, dimension(:), target :: g_pn
   integer(int32), save, public :: n_nz_epa
   integer(int32), save, public :: n_sym_epa
   real(dp), save, public, allocatable, dimension(:,:)   :: e_T, p_T
   real(dp), save, public, allocatable, dimension(:,:)   :: e_Q, p_Q
   integer(int32), save, public :: n_nz_ep
   integer(int32), save, public :: n_sym_ep
   real(dp), save, public, allocatable, dimension(:,:) :: T_e
   real(dp), save, public, allocatable, dimension(:,:) :: T_p
   real(dp), save, public, allocatable, dimension(:,:) :: Q
   real(dp), save, public, allocatable, dimension(:,:) :: Z_ee_u
   real(dp), save, public, allocatable, dimension(:,:,:) :: Z_ee_u_b
   real(dp), save, public, allocatable, dimension(:,:)   :: vec_da_oeu_s
   real(dp), save, public, allocatable, dimension(:,:)   :: vec_da_oed_s
   real(dp), save, public, allocatable, dimension(:,:)   :: vec_da_op_s
   real(dp), save, public, allocatable, dimension(:,:)   :: vec_oeu_s
   real(dp), save, public, allocatable, dimension(:,:)   :: vec_oed_s
   real(dp), save, public, allocatable, dimension(:,:)   :: vec_op_s
   real(dp), save, public, allocatable, dimension(:,:)   :: epo_s_dif
   real(dp), save, public, allocatable, dimension(:,:,:) :: mat_epo_s
   public :: comp_fermions_fermions_atomic_factors_symtable
contains
   subroutine comp_fermions_fermions_atomic_factors_symtable
      type(lmbd_ele_t), allocatable, dimension(:) :: nz_G_tmp
      type(lmbd_sym_t), allocatable, dimension(:) :: sym_G_tmp
      integer(int32) :: i1
      if ( n_el.ge.2_int32 .or. n_po.ge.2_int32 ) then
         n_nz_G = n_jorbs**2
         allocate( nz_G_tmp(1:n_nz_G) )
         allocate(sym_G_tmp(1:n_nz_G))
         do i1 = 1_int32, n_nz_G
            sym_G_tmp(i1)%is_indp = .true.
            sym_G_tmp(i1)%n_c = 1_int32
            allocate( sym_G_tmp(i1)%c(1:2*n_sys_sym) )
            sym_G_tmp(i1)%c(1) = i1
         enddo
         n_sym_G_ee = n_nz_G
         call cmp_nzr_fee( cpl_j, n_jorbs, jorbs, n_nz_G, nz_G_tmp )
         call cmp_sgm_sym( .true., 'R', 's', n_jorbs, jorbs, n_nz_G, nz_G_tmp, n_sym_G_ee, sym_G_tmp  )
         allocate(nz_G(1:n_nz_G)) ; nz_G(1:n_nz_G) = nz_G_tmp(1:n_nz_G)
         deallocate(nz_G_tmp)
         allocate(sym_G_ee(1:n_sym_G_ee))
         do i1 = 1_int32, n_sym_G_ee
            sym_G_ee(i1)%n_c = sym_G_tmp(i1)%n_c
            allocate( sym_G_ee(i1)%c(1:sym_G_ee(i1)%n_c) )
            sym_G_ee(i1)%c(1:sym_G_ee(i1)%n_c) = sym_G_tmp(i1)%c(1:sym_G_ee(i1)%n_c)
         enddo
         deallocate(sym_G_tmp)
      endif 
      if (n_po.gt.0_int32.or.spin_j.eq.'U') then
         allocate(sym_G_tmp(1:n_jorbs**2))
         do i1 = 1_int32, n_jorbs**2
            sym_G_tmp(i1)%is_indp = .true.
            sym_G_tmp(i1)%n_c = 1_int32
            allocate( sym_G_tmp(i1)%c(1:2*n_sys_sym) )
            sym_G_tmp(i1)%c(1) = i1
         enddo
         if ( .not.allocated(nz_G) ) then
            n_nz_G = n_jorbs**2
            allocate( nz_G_tmp(1:n_nz_G) )
            n_sym_G_pe = n_nz_G
            call cmp_nzr_fee( cpl_j, n_jorbs, jorbs, n_nz_G, nz_G_tmp )
            call cmp_sgm_sym( .true., 'U', 's', n_jorbs, jorbs, n_nz_G, nz_G_tmp, n_sym_G_pe, sym_G_tmp )
            allocate(nz_G(1:n_nz_G)) ; nz_G(1:n_nz_G) = nz_G_tmp(1:n_nz_G)
            deallocate(nz_G_tmp)
         else
            n_sym_G_pe = n_nz_G
            call cmp_sgm_sym( .false., 'U', 's', n_jorbs, jorbs, n_nz_G, nz_G, n_sym_G_pe, sym_G_tmp )
         endif
         allocate(sym_G_pe(1:n_sym_G_pe))
         do i1 = 1_int32, n_sym_G_pe
            sym_G_pe(i1)%n_c = sym_G_tmp(i1)%n_c
            allocate( sym_G_pe(i1)%c(1:sym_G_pe(i1)%n_c) )
            sym_G_pe(i1)%c(1:sym_G_pe(i1)%n_c) = sym_G_tmp(i1)%c(1:sym_G_pe(i1)%n_c)
         enddo
         deallocate(sym_G_tmp)
      endif
   end subroutine comp_fermions_fermions_atomic_factors_symtable
end module orbitalsbased_jastrow_factors_params_m
