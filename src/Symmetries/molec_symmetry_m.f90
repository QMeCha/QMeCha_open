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
module molec_symmetry_m
   use fortran_kinds_v, only: int32, sp
   implicit none
   logical, save, public :: atm_rot
   logical, save, public :: atm_trs
   logical, save, public :: sys_rot
   logical, save, public :: atm_eqv
   integer(int32), save, public                                :: n_sym_op
   real(sp),       save, public, allocatable, dimension(:,:,:) :: sym_op
   public :: init_symmetry_operations, fnlz_symmetry_operations
contains
   subroutine init_symmetry_operations()
      integer(int32) :: is
      n_sym_op = 49_int32
      allocate( sym_op(1:3,1:3,1:n_sym_op) ) ; sym_op = 0.0_sp
      is = 1_int32
      sym_op(1,1,is) =  1.0_sp
      sym_op(2,2,is) =  1.0_sp
      sym_op(3,3,is) =  1.0_sp
      is = 2_int32
      sym_op(1,1,is) = -1.0_sp
      sym_op(2,2,is) = -1.0_sp
      sym_op(3,3,is) =  1.0_sp
      is = 3_int32
      sym_op(1,1,is) = -1.0_sp
      sym_op(2,2,is) =  1.0_sp
      sym_op(3,3,is) = -1.0_sp
      is = 4_int32
      sym_op(1,1,is) =  1.0_sp
      sym_op(2,2,is) = -1.0_sp
      sym_op(3,3,is) = -1.0_sp
      is = 5_int32
      sym_op(1,2,is) =  1.0_sp
      sym_op(2,1,is) =  1.0_sp
      sym_op(3,3,is) = -1.0_sp
      is = 6_int32
      sym_op(1,2,is) = -1.0_sp
      sym_op(2,1,is) = -1.0_sp
      sym_op(3,3,is) = -1.0_sp
      is = 7_int32
      sym_op(1,2,is) = -1.0_sp
      sym_op(2,1,is) =  1.0_sp
      sym_op(3,3,is) =  1.0_sp
      is = 8_int32
      sym_op(1,2,is) =  1.0_sp
      sym_op(2,1,is) = -1.0_sp
      sym_op(3,3,is) =  1.0_sp
      is = 9_int32
      sym_op(1,3,is) =  1.0_sp
      sym_op(2,2,is) = -1.0_sp
      sym_op(3,1,is) =  1.0_sp
      is = 10_int32
      sym_op(1,3,is) = -1.0_sp
      sym_op(2,2,is) = -1.0_sp
      sym_op(3,1,is) = -1.0_sp
      is = 11_int32
      sym_op(1,3,is) =  1.0_sp
      sym_op(2,2,is) =  1.0_sp
      sym_op(3,1,is) = -1.0_sp
      is = 12_int32
      sym_op(1,3,is) = -1.0_sp
      sym_op(2,2,is) =  1.0_sp
      sym_op(3,1,is) =  1.0_sp
      is = 13_int32
      sym_op(1,1,is) = -1.0_sp
      sym_op(2,3,is) =  1.0_sp
      sym_op(3,2,is) =  1.0_sp
      is = 14_int32
      sym_op(1,1,is) = -1.0_sp
      sym_op(2,3,is) = -1.0_sp
      sym_op(3,2,is) = -1.0_sp
      is = 15_int32
      sym_op(1,1,is) =  1.0_sp
      sym_op(2,3,is) =  1.0_sp
      sym_op(3,2,is) = -1.0_sp
      is = 16_int32
      sym_op(1,1,is) =  1.0_sp
      sym_op(2,3,is) = -1.0_sp
      sym_op(3,2,is) =  1.0_sp
      is = 17_int32
      sym_op(1,2,is) =  1.0_sp
      sym_op(2,3,is) =  1.0_sp
      sym_op(3,1,is) =  1.0_sp
      is = 18_int32
      sym_op(1,2,is) = -1.0_sp
      sym_op(2,3,is) =  1.0_sp
      sym_op(3,1,is) = -1.0_sp
      is = 19_int32
      sym_op(1,2,is) =  1.0_sp
      sym_op(2,3,is) = -1.0_sp
      sym_op(3,1,is) = -1.0_sp
      is = 20_int32
      sym_op(1,2,is) = -1.0_sp
      sym_op(2,3,is) = -1.0_sp
      sym_op(3,1,is) =  1.0_sp
      is = 21_int32
      sym_op(1,3,is) =  1.0_sp
      sym_op(2,1,is) =  1.0_sp
      sym_op(3,2,is) =  1.0_sp
      is = 22_int32
      sym_op(1,3,is) =  1.0_sp
      sym_op(2,1,is) = -1.0_sp
      sym_op(3,2,is) = -1.0_sp
      is = 23_int32
      sym_op(1,3,is) = -1.0_sp
      sym_op(2,1,is) = -1.0_sp
      sym_op(3,2,is) =  1.0_sp
      is = 24_int32
      sym_op(1,3,is) = -1.0_sp
      sym_op(2,1,is) =  1.0_sp
      sym_op(3,2,is) = -1.0_sp
      is = 25_int32
      sym_op(1,1,is) = -1.0_sp
      sym_op(2,2,is) = -1.0_sp
      sym_op(3,3,is) = -1.0_sp
      is = 26_int32
      sym_op(1,1,is) =  1.0_sp
      sym_op(2,2,is) =  1.0_sp
      sym_op(3,3,is) = -1.0_sp
      is = 27_int32
      sym_op(1,1,is) =  1.0_sp
      sym_op(2,2,is) = -1.0_sp
      sym_op(3,3,is) =  1.0_sp
      is = 28_int32
      sym_op(1,1,is) = -1.0_sp
      sym_op(2,2,is) =  1.0_sp
      sym_op(3,3,is) = -1.0_sp
      is = 29_int32
      sym_op(1,2,is) = -1.0_sp
      sym_op(2,1,is) = -1.0_sp
      sym_op(3,3,is) =  1.0_sp
      is = 30_int32
      sym_op(1,2,is) =  1.0_sp
      sym_op(2,1,is) =  1.0_sp
      sym_op(3,3,is) =  1.0_sp
      is = 31_int32
      sym_op(1,2,is) = -1.0_sp
      sym_op(2,1,is) =  1.0_sp
      sym_op(3,3,is) = -1.0_sp
      is = 32_int32
      sym_op(1,2,is) =  1.0_sp
      sym_op(2,1,is) = -1.0_sp
      sym_op(3,3,is) = -1.0_sp
      is = 33_int32
      sym_op(1,3,is) = -1.0_sp
      sym_op(2,2,is) =  1.0_sp
      sym_op(3,1,is) = -1.0_sp
      is = 34_int32
      sym_op(1,3,is) =  1.0_sp
      sym_op(2,2,is) =  1.0_sp
      sym_op(3,1,is) =  1.0_sp
      is = 35_int32
      sym_op(1,3,is) = -1.0_sp
      sym_op(2,2,is) = -1.0_sp
      sym_op(3,1,is) =  1.0_sp
      is = 36_int32
      sym_op(1,3,is) =  1.0_sp
      sym_op(2,2,is) = -1.0_sp
      sym_op(3,1,is) = -1.0_sp
      is = 37_int32
      sym_op(1,1,is) =  1.0_sp
      sym_op(2,3,is) = -1.0_sp
      sym_op(3,2,is) = -1.0_sp
      is = 38_int32
      sym_op(1,1,is) = 1.0_sp
      sym_op(2,3,is) = 1.0_sp
      sym_op(3,2,is) = 1.0_sp
      is = 39_int32
      sym_op(1,1,is) = -1.0_sp
      sym_op(2,3,is) = -1.0_sp
      sym_op(3,2,is) =  1.0_sp
      is = 40_int32
      sym_op(1,1,is) = -1.0_sp
      sym_op(2,3,is) =  1.0_sp
      sym_op(3,2,is) = -1.0_sp
      is = 41_int32
      sym_op(1,2,is) = -1.0_sp
      sym_op(2,3,is) = -1.0_sp
      sym_op(3,1,is) = -1.0_sp
      is = 42_int32
      sym_op(1,2,is) =  1.0_sp
      sym_op(2,3,is) = -1.0_sp
      sym_op(3,1,is) =  1.0_sp
      is = 43_int32
      sym_op(1,2,is) = -1.0_sp
      sym_op(2,3,is) =  1.0_sp
      sym_op(3,1,is) =  1.0_sp
      is = 44_int32
      sym_op(1,2,is) =  1.0_sp
      sym_op(2,3,is) =  1.0_sp
      sym_op(3,1,is) = -1.0_sp
      is = 45_int32
      sym_op(1,3,is) = -1.0_sp
      sym_op(2,1,is) = -1.0_sp
      sym_op(3,2,is) = -1.0_sp
      is = 46_int32
      sym_op(1,3,is) = -1.0_sp
      sym_op(2,1,is) =  1.0_sp
      sym_op(3,2,is) =  1.0_sp
      is = 47_int32
      sym_op(1,3,is) =  1.0_sp
      sym_op(2,1,is) =  1.0_sp
      sym_op(3,2,is) = -1.0_sp
      is = 48_int32
      sym_op(1,3,is) =  1.0_sp
      sym_op(2,1,is) = -1.0_sp
      sym_op(3,2,is) =  1.0_sp
      is = 49_int32
      sym_op(1,1,is) = -1.0_sp
      sym_op(2,2,is) =  1.0_sp
      sym_op(3,3,is) =  1.0_sp
   end subroutine init_symmetry_operations
   subroutine fnlz_symmetry_operations
      if ( allocated(sym_op) ) deallocate( sym_op )
   end subroutine fnlz_symmetry_operations
end module molec_symmetry_m
