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
program QMeCha
   use openmp_mpi_m,             only: init_ompmpi_env, fnlz_ompmpi_env
   use timings_m,                only: init_simulation_timings, fnlz_simulation_timings
   use io_datasheet_m,           only: read_data_fle
   use io_arguments_m,           only: read_command_args
   use mersenne_twister19937_m,  only: init_mt19937 
   use molecular_system_v,       only: molsys_prs
   use molecular_system_m,       only: read_coord_fle, init_molec_sys
   use qdo_system_v,             only: qdosys_prs
   use qdo_system_m,             only: read_qcoord_fle, ini_qdo_sys
   use io_fermionic_basisset_m,  only: read_bssst_fle, dall_bssst_var
   use io_drudonic_basisset_m,   only: read_qbssst_fle
   use pseudopotentials_v,       only: read_psdpot_fle
   use pseudopotentials_m,       only: init_pseudopotentials
   use molec_symmetry_m,         only: init_symmetry_operations, fnlz_symmetry_operations
   use fermionic_wavefunction_m, only: init_fermionic_wavefunction_v
   use qdo_wavefunction_m,       only: ini_qdowvf_par
   use electronic_properties_m,  only: init_mlt_mmn
   use cubic_harmonics_m,        only: init_cubic_hrm
   use quantum_monte_carlo_v,    only: n_wlk, n_wlk_max
   use quantum_monte_carlo_m,    only: init_qmc_run, exec_qmc_run
   use local_energy_m,           only: ini_locene
   implicit none
   call init_ompmpi_env()
   call qmecha_version()
   call init_simulation_timings()
   call read_command_args()
   call read_coord_fle()
   if ( molsys_prs ) call read_psdpot_fle()
   if ( molsys_prs ) call read_bssst_fle()
   call read_qcoord_fle()
   if(qdosys_prs) call read_qbssst_fle()
   call init_mt19937( n_wlk, n_wlk_max )
   call init_symmetry_operations()
   if(molsys_prs) then
      call init_pseudopotentials()
      call init_molec_sys()
      call init_mlt_mmn()
   endif
   if(qdosys_prs) call ini_qdo_sys()
   call init_cubic_hrm()
   if(molsys_prs) then
      call init_fermionic_wavefunction_v()
      call dall_bssst_var()
   endif
   if (qdosys_prs) then
      call ini_qdowvf_par()
   endif
   call ini_locene()
   call fnlz_symmetry_operations()
   call init_qmc_run()
   call exec_qmc_run()
   call fnlz_simulation_timings()
   call fnlz_ompmpi_env()
end program QMeCha
