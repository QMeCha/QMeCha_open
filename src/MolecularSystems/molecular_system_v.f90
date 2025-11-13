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
module molecular_system_v
   use fortran_kinds_v, only: int32, dp, sp
   implicit none
   logical, public, save :: molsys_prs
   type, public :: atom_t
      character(4)   :: atm_name
      integer(int32) :: atm_z
      integer(int32) :: i_bs
      integer(int32) :: i_pp
   end type atom_t
   type, public :: dist_t
      real(dp), dimension(3) :: v
      real(dp)                 :: m
   end type dist_t
   integer(int32), save, public, target :: n_el
   integer(int32), save, public, target :: n_el_u
   integer(int32), save, public, target :: n_el_d
   integer(int32), save, public, target :: n_el_s
   integer(int32), save, public, target :: n_po
   integer(int32), save, public, target :: n_po_u
   integer(int32), save, public, target :: n_po_d
   integer(int32), save, public, target :: n_po_s
   integer(int32), save, public, target :: n_fe
   integer(int32), save, public, target :: n_at
   integer(int32), save, public, target :: chrg
   integer(int32), save, public, target :: nuc_chrg
   integer(int32), save, public, target :: mult
   real(dp),       save, public, target :: spin_e
   real(dp),       save, public, target :: spin_p
   real(dp),       save, public, target :: m_e
   real(dp),       save, public, target :: m_p
   real(dp),       save, public, target :: mu_ep
   real(dp),       save, public, allocatable, dimension(:,:), target :: r_at
   type(atom_t),   save, public, allocatable, dimension(:),   target :: atoms
   type(dist_t),   save, public, allocatable, dimension(:),   target :: d_nn
   integer(int32), save, public, allocatable, dimension(:) :: f_spin
   real(sp),       save, public, allocatable, dimension(:) :: f_chrg
   character(len=20), save, public :: fle_coord
   integer(int32), save, public                                      :: n_indp_at
   integer(int32), save, public                                      :: n_sys_sym
   integer(int32), save, public, target, allocatable, dimension(:,:) :: sys_sym_tab
   type, public :: atmsym_t
      integer(int32)                            :: n_cnct_at
      integer(int32), allocatable, dimension(:) :: cnct
   end type atmsym_t
   type(atmsym_t), save, public, allocatable, dimension(:), target :: atm_sym
end module molecular_system_v
