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
module quantum_monte_carlo_v
   use fortran_kinds_v, only: dp, int32, stdout
   use openmp_mpi_m
   use write_lines_m
   implicit none
   character(100), save, public :: sysname
   logical, save, public :: restart
   integer(int32), save, public :: save_status
   integer(int32), save, public :: save_step
   integer(int32), save, public :: n_wlk_max
   integer(int32), save, public :: n_tot_wlk
   real(dp),       save, public :: G_w_trgt
   integer(int32), save, public :: n_wlk
   integer(int32), save, public :: bin_l
   integer(int32), save, public :: n_bin
   integer(int32), save, public :: bin_l_old
   integer(int32), save, public :: n_bin_old
   integer(int32), public, save, target :: n_tot_smpl, n_tot_smpl_old
   integer(int32), save, public :: n_trm
   integer(int32), save, public :: n_bra, n_bra_vmc
   integer(int32), save, public :: n_scr
   character(3),   save, public :: qmc_mthd, qmc_mthd_old
   character(3),   save, public :: smp_mthd
   logical,        save, public :: spc_map
   logical,        save, public :: two_step_sampling
   integer(int32), save, public :: smp_type
   integer(int32), save, public :: mapping_type
   logical,        save, public :: nuc_corr
   real(dp)      , save, public, target :: dt_e, dt_p, dt_d, dt
   real(dp)      , save, public :: arate, a_drft
   logical       , save, public :: var_dt, grad_updat
   public :: save_qmc_var, load_qmc_var
contains
   subroutine save_qmc_var( step )
      integer(int32), intent(in) :: step
      character(100) :: svf_fle
      if (mpi_rank.eq.0_int32) then
         svf_fle = 'qmc.save/qmc_variables.sav'
         open(unit=1,file=svf_fle,action='write',form='unformatted',status='unknown',access='sequential')
         write(1) save_status, qmc_mthd
         write(1) n_wlk, n_wlk_max, n_tot_wlk
         write(1) bin_l, step+n_bin_old, n_trm, n_bra , n_bra_vmc, n_scr
         write(1) smp_mthd, spc_map, two_step_sampling, smp_type, mapping_type, nuc_corr
         write(1) dt_e, dt_p, dt_d, dt
         write(1) arate, a_drft
         write(1) var_dt, grad_updat
         close(1)
      endif
   end subroutine save_qmc_var
   subroutine load_qmc_var( )
      character(100) :: svf_fle
      logical        :: file_present
      integer(int32) :: n_wlk_old, n_wlk_max_old, n_tot_wlk_old
      file_present = .false.
      if (mpi_rank.eq.0_int32) then
         svf_fle = 'qmc.save/qmc_variables.sav'
         inquire(file=svf_fle, exist=file_present)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(file_present, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif
      if (file_present) then
         call write_empty_line(stdout,0,mpi_rank)
         call write_line_check(stdout,0,mpi_rank,2,"Loading quantum Monte Carlo variables")
         if (mpi_rank.eq.0_int32) then
            open(unit=1,file=svf_fle,action='read',form='unformatted',status='old',access='sequential')
            read(1) save_status, qmc_mthd_old
            read(1) n_wlk_old, n_wlk_max_old, n_tot_wlk_old
            read(1) bin_l_old, n_bin_old, n_trm, n_bra, n_bra_vmc, n_scr
            read(1) smp_mthd, spc_map, two_step_sampling, smp_type, mapping_type, nuc_corr
            read(1) dt_e, dt_p, dt_d, dt
            read(1) arate, a_drft
            read(1) var_dt, grad_updat
            close(1)
         endif
         call write_done(stdout,0,mpi_rank)
         call write_simple_line(stdout,0,mpi_rank,2,"c","__________________________")
         call write_simple_line(stdout,0,mpi_rank,2,"c","Found previous calculation")
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_line(stdout,0,mpi_rank,2,"Previous calculation", qmc_mthd_old,var_name="qmc_mthd")
         select case ( save_status )
          case(0)
            call write_variable_line(stdout,0,mpi_rank,2,"Status, Need restart from scratch", save_status)
          case(1)
            call write_variable_line(stdout,0,mpi_rank,2,"Status, Restart after VMC thermalization", save_status)
          case(2)
            call write_variable_line(stdout,0,mpi_rank,2,"Status, Restart after completed (or partial) VMC calculation", save_status)
          case(3)
            call write_variable_line(stdout,0,mpi_rank,2,"Status, Restart after DMC thermalization", save_status)
          case(4)
            call write_variable_line(stdout,0,mpi_rank,2,"Status, Restart after completed (or partial) DMC calculation", save_status)
          case default
         end select
         call write_variable_line(stdout,0,mpi_rank,2,"Total number of walkers", n_tot_wlk_old,var_name="n_tot_wlk")
         call write_variable_line(stdout,0,mpi_rank,2,"Number of independent walkers per MPI task", n_wlk_old,var_name="n_wlk")
         call write_variable_line(stdout,0,mpi_rank,2,"Maximum number of independent walkers per MPI task", n_wlk_max_old,var_name="n_wlk_max")
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_bcast(save_status,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(  qmc_mthd_old, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(n_tot_smpl_old, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(     bin_l_old, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(     n_bin_old, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
         if ( bin_l_old .ne. bin_l) then
            bin_l = bin_l_old
            n_bin = int(n_tot_smpl / bin_l_old)
            call write_empty_line(stdout,0,mpi_rank)
            call write_variable_line(stdout,0,mpi_rank,2,"Resetting bin length to previous calculation", bin_l,var_name="bin_l")
            call write_variable_line(stdout,0,mpi_rank,2,"Setting number of bins", n_bin,var_name="n_bin")
         endif
         n_tot_smpl = n_bin * bin_l
         n_tot_smpl_old = bin_l_old * n_bin_old
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_bcast(n_wlk,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(n_tot_wlk,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(n_wlk_max,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(n_trm,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(n_bra,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(n_bra_vmc,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(n_scr,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(qmc_mthd,     3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(spc_map ,     1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(smp_mthd,     3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(two_step_sampling,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(smp_type,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(mapping_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(nuc_corr,     1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(dt,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(dt_e,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(dt_p,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(dt_d,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(arate,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(a_drft,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(var_dt ,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(grad_updat,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
#endif
      else
         call write_empty_line(stdout,0,mpi_rank)
         call write_simple_line(stdout,0,mpi_rank,2,"l","No previous calculation found")
         save_status = 0_int32
      endif
   end subroutine load_qmc_var
end module quantum_monte_carlo_v
