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
module jastrow_orbitals_m
   use fortran_kinds_v, only: int32, dp, stdout
   use openmp_mpi_m, only: mpi_rank
   use write_lines_m
   use jastrow_factors_v, only: n_par_jst, n_par_jst_s
   use orbbss_cls, only: sysorb_t, orbspa_t
   use orbmat_cls, only: matorb_t
   use basisset_v, only: jst_bsst
   use quantum_monte_carlo_v, only: n_wlk_max
   implicit none
   integer(int32),  save, public,              pointer :: n_jorbs
   integer(int32),  save, public,              pointer :: n_par_jorbs,  n_par_jorbs_s
   logical,         save, public :: jc_opt, je_opt
   type(orbspa_t),  save, public, dimension(:), pointer :: jorbs
   type(sysorb_t),  save, public, target :: jorbs_s
   type(matorb_t),  save, public, allocatable, dimension(:), target :: jorbs_mat
   public :: ini_jorbs_mat
contains
   subroutine ini_jorbs_mat( )
      integer(int32) :: iw
      call jorbs_s%ini( jst_bsst, jc_opt, je_opt )
      jorbs         => jorbs_s%orbs
      n_jorbs       => jorbs_s%n_orbs
      n_par_jorbs   => jorbs_s%n_par_orbs
      n_par_jorbs_s => jorbs_s%n_par_orbs_s
      n_par_jst   = n_par_jorbs
      n_par_jst_s = n_par_jorbs_s
      if ( jorbs_s%n_orbs.gt.0_int32) then
         allocate( jorbs_mat(1:n_wlk_max) )
         do iw = 1_int32, n_wlk_max
            call jorbs_mat(iw)%ini( iw, jorbs_s )
         enddo  
      endif
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of atomic orbitals for Jastrow part", n_jorbs)
   end subroutine ini_jorbs_mat
end module jastrow_orbitals_m
