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
module io_chgoptpar_m
   use fortran_kinds_v,             only: dp, int32
   use openmp_mpi_m
   use wavefunction_optimization_v, only: da, eps, n_opt_step, eps_auto_tune
   use quantum_monte_carlo_v,       only: bin_l
   use write_lines_m
   implicit none
   public :: read_chng_fle, write_chng_fle, delete_chng_fle
contains
   subroutine write_chng_fle
      if ( mpi_rank .eq. 0_int32 ) then
         open(unit=1,file='opt.chck',action='write',form='formatted',status='unknown')
         write(1,'("&optpar")')
         write(1,'("       da = ",E15.6)') da
         if(.not.eps_auto_tune) write(1,'("      eps = ",E15.6)') eps
         write(1,'("    bin_l = ",I8)') bin_l
         write(1,'(" stop_opt = .false.")')
         write(1,'("/")')
         close(1)
      endif
   end subroutine write_chng_fle
   subroutine delete_chng_fle
      if ( mpi_rank .eq. 0_int32 ) call system("rm -rf opt.chck")
   end subroutine delete_chng_fle
   subroutine read_chng_fle( stop_opt )
      logical, intent(inout) :: stop_opt
      logical                :: file_exists
      real(dp)               :: da_old, eps_old
      integer(int32)         :: bin_l_old
      integer(int32)         :: ierr
      namelist /optpar/ da, eps, stop_opt, bin_l
      stop_opt = .false.
      ierr        = 0
      file_exists = .false.
      if ( mpi_rank .eq. 0_int32 ) then
         da_old     = da
         eps_old    = eps
         bin_l_old  = bin_l
         inquire(file='opt.chck', exist=file_exists)
         if( file_exists  ) then
            open(unit=1,file='opt.chck',action='read',form='formatted',status='unknown')
            read(1,nml=optpar,IOSTAT=ierr)
            if ( ierr.ne.0 ) then
               write(*,'(3X,"ERROR while reading opt.chck file. (Ignoring changes)            ")')
               da     = da_old
               eps    = eps_old
               bin_l  = bin_l_old
            else
               !if ( eps_auto_tune ) eps = eps_old
               if ( da.ne.da_old ) call write_variable_line(stdout,0,mpi_rank,2,"Changing parameter da to", da,var_name="da")
               !if ( (.not.eps_auto_tune).and.(eps.ne.eps_old) ) call write_variable_line(stdout,0,mpi_rank,2,"Changing regularization eps to", eps,var_name="eps")
               if ( bin_l.ne.bin_l_old ) call write_variable_line(stdout,0,mpi_rank,2,"Changing bin length to", bin_l,var_name="bin_l")
            endif
            close(1)
         endif
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(file_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(ierr,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      if ( file_exists.and.(ierr.eq.0) ) then
         call mpi_bcast(bin_l,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(da   ,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(eps  ,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(stop_opt , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      endif
#endif
   end subroutine read_chng_fle
end module io_chgoptpar_m
