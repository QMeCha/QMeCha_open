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
module fermionic_orbitals_m
   use fortran_kinds_v, only: int32, dp, stdout
   use openmp_mpi_m, only: mpi_rank
   use write_lines_m
   use fermionic_wavefunction_v, only: n_par_frm, n_par_frm_s
   use orbbss_cls, only: sysorb_t, orbspa_t
   use orbmat_cls, only: matorb_t
   use quantum_monte_carlo_v, only: n_wlk_max
   use basisset_v, only: frm_bsst
   implicit none
   integer(int32), save, public,               pointer :: n_forbs
   integer(int32), save, public,               pointer :: n_par_forbs
   integer(int32), save, public,               pointer :: n_par_forbs_s
   logical,        save, public                        :: aoc_opt, aoe_opt
   type(orbspa_t), save, public, dimension(:), pointer :: forbs
   type(sysorb_t), save, public, target :: forbs_s
   type(matorb_t), save, public, allocatable, dimension(:), target :: forbs_mat
   public  :: ini_forbs_mat
contains
   subroutine ini_forbs_mat( )
      integer(int32) :: iw
      call forbs_s%ini( frm_bsst, aoc_opt, aoe_opt )
      forbs          => forbs_s%orbs
      n_forbs        => forbs_s%n_orbs
      n_par_forbs   => forbs_s%n_par_orbs
      n_par_forbs_s => forbs_s%n_par_orbs_s
      n_par_frm   = n_par_forbs
      n_par_frm_s = n_par_forbs_s
      if ( forbs_s%n_orbs.gt.0_int32) then
         allocate( forbs_mat(1:n_wlk_max) )
         do iw = 1_int32, n_wlk_max
            call forbs_mat(iw)%ini( iw, forbs_s )
         enddo
      endif
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of atomic orbitals for fermionic part", n_forbs)
   end subroutine ini_forbs_mat
end module fermionic_orbitals_m
