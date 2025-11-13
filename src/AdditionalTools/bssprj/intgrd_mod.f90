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
module intgrd_mod
use fortran_kinds_v, only: dp, int32
use molecular_system_v, only: dist_t
implicit none
integer(int32), public, save :: n_int_points
real(dp)      , public, save :: orb_cut_off, nuc_cut_off
real(dp)      , public, save :: dorb 
real(dp)      , public, save, allocatable, dimension(:,:) :: trgt_grid
end module intgrd_mod