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
module basisset_v
   use fortran_kinds_v,  only: int32, dp
   use atomic_orbital_c, only: orb_t
   implicit none
   type, public :: bssset_t
      character(4)                             :: atm_name
      integer(int32)                           :: n_shlls
      integer(int32)                           :: n_orbs
      type(orb_t),   allocatable, dimension(:) :: shlls
   end type bssset_t
   integer(int32), public, save :: n_fbs
   integer(int32), public, save :: n_qbs
   type(bssset_t), public, save, target, allocatable, dimension(:) :: jst_bsst
   type(bssset_t), public, save, target, allocatable, dimension(:) :: frm_bsst
   type(bssset_t), public, save, target, allocatable, dimension(:) :: drd_bsst
   type(bssset_t), public, save, target, allocatable, dimension(:) :: jdr_bsst
   character(len=20), public ,save :: fle_basis
   character(len=20), public ,save :: fle_basis_qdo
end module basisset_v
