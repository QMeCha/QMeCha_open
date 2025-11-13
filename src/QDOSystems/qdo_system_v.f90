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
module qdo_system_v
  use fortran_kinds_v,    only: int8, int32, dp
  use molecular_system_v, only: dist_t
  implicit none
  logical, public ,save :: qdosys_prs
  type, public :: qdo_t
    character(4)   :: qdo_name
    real(dp)       :: qdo_q
    real(dp)       :: qdo_omega
    real(dp)       :: qdo_m 
    real(dp)       :: qdo_sigma_c
    real(dp)       :: qdo_sigma_d
    integer(int32) :: i_bsQ
  end type qdo_t
  type, public :: pchrg_t
    character(3)   :: pchrg_name 
    real(dp)       :: pchrg_q
    real(dp)       :: pchrg_sigma
  end type pchrg_t
  real(dp),       save, public, target :: el_sigma
  integer(int32), save, public, target :: n_qdo
  integer(int32), save, public, target :: n_pchrg
  character(3),   save, public, target :: qdoham
  real(dp),       save, public, allocatable, dimension(:,:), target :: r_qdo
  type(qdo_t),    save, public, allocatable, dimension(:),   target :: qdos
  type(dist_t),   save, public, allocatable, dimension(:),   target :: d_qq
  real(dp),       save, public, allocatable, dimension(:,:), target :: r_pchrg
  type(pchrg_t),  save, public, allocatable, dimension(:),   target :: pchrgs
  type(dist_t),   save, public, allocatable, dimension(:),   target :: d_pchpch
  type(dist_t),   save, public, allocatable, dimension(:,:), target :: d_qpch
  real(dp),  save, public, allocatable, dimension(:) :: d_chrg
  real(dp),  save, public, allocatable, dimension(:) :: d_mass
  character(len=20), save, public :: fle_coord_qdo
  ! ########################################################################
  type(dist_t),   save, public, allocatable, dimension(:,:), target :: d_qn
  type(dist_t),   save, public, allocatable, dimension(:,:), target :: d_npch
end module qdo_system_v
