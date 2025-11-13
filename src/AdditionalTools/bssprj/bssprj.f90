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
program bssprj
use fortran_kinds_v, only: int32
use openmp_mpi_m, only: init_ompmpi_env, fnlz_ompmpi_env
use argpio_mod, only: read_comm_args
use datpio_mod, only: trg_basis, new_basis
use basisset_v, only: fle_basis 
use io_fermionic_basisset_m, only: read_bssst_fle, dall_bssst_var, prnt_bssst_fle
use mersenne_twister19937_m,  only: init_mt19937
use molec_symmetry_m, only: init_symmetry_operations
use molecular_system_m, only: read_coord_fle, init_molec_sys
use trgtgo_mod, only: ini_trgt_orbs
use optorb_mod, only: ini_opt_orbs, run_opt_orbs, prnt_opt_orbs
use quantum_monte_carlo_v, only: n_wlk, n_wlk_max
implicit none
integer(int32) :: ierr 
call init_ompmpi_env( )
call read_comm_args( ierr )
if ( ierr .ne. 0_int32 ) then
  call fnlz_ompmpi_env( )
  stop
endif
call init_mt19937(n_wlk, n_wlk_max  )
call init_symmetry_operations()
call read_coord_fle(  )
call init_molec_sys()
fle_basis = trg_basis
call read_bssst_fle( )
call ini_trgt_orbs()
call dall_bssst_var()
fle_basis = new_basis
call read_bssst_fle( )
call ini_opt_orbs()
call run_opt_orbs()
call prnt_opt_orbs()
call prnt_bssst_fle()
call fnlz_ompmpi_env()
end program bssprj