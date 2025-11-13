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
module drudonic_orbitals_m
   use fortran_kinds_v,  only: int32, dp
   use qdo_wavefunction_v, only: n_par_drd, n_par_drd_s
   use orbbssq_cls, only: qsysorb_t
   use orbmaq_cls,  only: matorbq_t
   use orbspq_cls,  only: orbspq_t
   implicit none
   logical,         save, public :: qc_opt, qe_opt
   integer(int32),  save, public,              pointer :: n_dorbs
   integer(int32),  save, public,              pointer :: n_par_dorbs
   integer(int32),  save, public,              pointer :: n_par_dorbs_s
   type(orbspq_t), save, public, dimension(:), pointer :: dorbs
   type(qsysorb_t),save, public,                            target :: dorbs_s
   type(matorbq_t),save, public, allocatable, dimension(:), target :: dorbs_mat
   public :: ini_dorbs_mat
contains
   subroutine ini_dorbs_mat( )
      use basisset_v, only: drd_bsst
      use quantum_monte_carlo_v, only: n_wlk_max
      integer(int32) :: i1
      call dorbs_s%ini( drd_bsst, qc_opt, qe_opt )
      dorbs         => dorbs_s%orbs
      n_dorbs       => dorbs_s%n_orbs
      n_par_dorbs   => dorbs_s%n_par_orbs
      n_par_dorbs_s => dorbs_s%n_par_orbs_s
      n_par_drd   = n_par_dorbs
      n_par_drd_s = n_par_dorbs_s
      if (n_dorbs.gt.0_int32) then
         allocate( dorbs_mat(1:n_wlk_max) )
         do i1 = 1_int32, n_wlk_max
            call dorbs_mat(i1)%ini( i1, dorbs_s )
         enddo
      endif
   end subroutine ini_dorbs_mat
end module drudonic_orbitals_m
