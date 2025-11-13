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
module openmp_mpi_m
   use fortran_kinds_v, only: int32, stdout
   use write_lines_m,   only: write_separator_line, write_simple_line, write_empty_line, &
   & write_variable_line
#ifdef _MPI
   use mpi
#endif
#ifdef _MPI08
   use mpi_f08
#endif
   implicit none
#ifdef _MPIh
   include 'mpif.h'
#endif
   integer(int32), public, save :: mpi_rank
   integer(int32), public, save :: n_mpi_tasks
   integer(int32), public, save :: n_omp_tasks
   integer       , public :: mpierr
#ifdef _OMP
   integer(int32), external :: omp_get_num_threads, omp_get_thread_num
#endif
   public :: init_ompmpi_env, fnlz_ompmpi_env
contains
   subroutine init_ompmpi_env()
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_init(mpierr)
      call mpi_comm_size(MPI_COMM_WORLD,n_mpi_tasks,mpierr)
      call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,mpierr)
      call mpi_barrier  (MPI_COMM_WORLD,mpierr)
#else
      n_mpi_tasks = 1_int32
      mpi_rank    = 0_int32
#endif
#ifdef _OMP
!$omp parallel default(shared) 
      n_omp_tasks = omp_get_num_threads()
!$omp end parallel 
#else
      n_omp_tasks = 1_int32
#endif
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","INITIALIZING MPI/OMP ENVIRONMENT")
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Number of MPI tasks", n_mpi_tasks,var_name="n_mpi_tasks")
      call write_variable_line(stdout,0,mpi_rank,2,"Number of OPENMP threads per MPI task", n_omp_tasks,var_name="n_omp_tasks")
   end subroutine init_ompmpi_env
   subroutine fnlz_ompmpi_env()
#if defined _MPI || defined _MPIh || defined _MPI08
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","FINALIZING MPI TASKS")
      call mpi_finalize(mpierr)
      call write_separator_line(stdout,0,mpi_rank,2,"=")
#endif
   end subroutine fnlz_ompmpi_env
end module openmp_mpi_m
