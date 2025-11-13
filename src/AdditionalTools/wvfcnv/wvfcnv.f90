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
program wvfcnv
use fortran_kinds_v, only: int32, stdout
use openmp_mpi_m, only: init_ompmpi_env, fnlz_ompmpi_env, mpi_rank
use write_lines_m
use basisset_v, only: fle_basis 
use io_fermionic_basisset_m, only: read_bssst_fle, dall_bssst_var, prnt_bssst_fle
use mersenne_twister19937_m,  only: init_mt19937
use cubic_harmonics_m, only: init_cubic_hrm
use molec_symmetry_m, only: init_symmetry_operations
use molecular_system_m, only: read_coord_fle, init_molec_sys
use fermionic_wavefunction_m, only: init_fermionic_wavefunction_v, ini_frmwvf_red, save_frmwvf_red
use io_datasheet_m, only: dflt_data_var, bcst_data_var
use wvfcnv_mod, only: eval_wvfcnv, bcst_input_wvfcnv, wf_type_e_new, wf_type_p_new, wf_name_cnv
use molecular_system_v, only: fle_coord
use pseudopotentials_m, only: init_pseudopotentials
use pseudopotentials_v, only: fle_pspot, read_psdpot_fle
use quantum_monte_carlo_v, only: restart, n_wlk, n_wlk_max
implicit none
integer(int32) :: i1
integer(int32) :: n_args
character(len=20) :: args
integer(int32), external :: wvfn_cod
call init_ompmpi_env( )
call init_symmetry_operations()
call dflt_data_var() 
call bcst_data_var()
if( mpi_rank.eq.0_int32 ) then
  call write_empty_line(stdout,0,mpi_rank)
  write(*,'(3X,"================================================================")')
  write(*,'(3X,"                    READING INPUT ARGUMENTS                     ")')
  call write_empty_line(stdout,0,mpi_rank)
  n_args = command_argument_count() 
  do i1 = 1_int32, n_args-1
    call get_command_argument(i1, args)
    select case ( trim(args) )
    case('-m') 
      call get_command_argument(i1+1, args)
      fle_coord = trim(args)
    case('-b') 
      call get_command_argument(i1+1, args)
      fle_basis = trim(args)
    case('-p') 
      call get_command_argument(i1+1, args)
      fle_pspot = trim(args)
    case('-w') 
      call get_command_argument(i1+1, args)
      wf_name_cnv = trim(args)
    end select
  enddo ! i1 (n_args)
endif ! mpi_rank.eq.0_int32
#if defined _MPI || defined _MPIh || defined _MPI08
  call bcst_input_wvfcnv()
#endif
call init_mt19937(n_wlk, n_wlk_max)
call read_coord_fle(  )
call read_psdpot_fle( )
call init_pseudopotentials()
call init_molec_sys()
call read_bssst_fle( )
call init_cubic_hrm()
restart = .true.
call init_fermionic_wavefunction_v()
wf_type_e_new = wvfn_cod( wf_name_cnv )
wf_type_p_new = wvfn_cod( wf_name_cnv )
call ini_frmwvf_red( wf_type_e_new, wf_type_p_new )
call eval_wvfcnv( )
call save_frmwvf_red(wf_type_e_new, wf_type_p_new )
call fnlz_ompmpi_env()
end program wvfcnv
