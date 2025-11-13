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
module timings_m
   use fortran_kinds_v, only: stdout, sp, dp, int32
   use openmp_mpi_m
   use write_lines_m
   implicit none
   character(8), target, public, save :: t_date
   character(10),target, public, save :: t_time
   real(dp),     target, public, save :: t_cpu_tmp
   real(dp),             public, save :: t_tot
   real(dp),             public, save :: t_therm_vmc, t_therm_dmc
   real(dp),             public, save :: t_dmc_step, t_dmc_tot
   real(dp),             public, save :: t_lnalg_step, t_lnalg_tot
   real(dp),             public, save :: t_opt_step, t_opt_tot
   real(dp),             public, save :: t_crs_step, t_crs_tot
   real(dp),             public, save :: t_avg_step, t_avg_tot
   real(dp),             public, save :: t_sampl_step, t_sampl_tot, t_resampl_tot
   real(dp),             public, save :: t_brnch_step, t_brnch_tot
   integer(int32), dimension(8), target, save, public :: t_intg_vec_tmp
   public  :: time_compute, time_compute_dateandtime, time_compute_cputime,&
   & init_simulation_timings, fnlz_simulation_timings
   private :: time_convert_intg_vec_to_sec
contains
   subroutine time_compute( time_sec, time_intg_vec )
      real(dp),                                       intent(inout) :: time_sec
      integer(int32), dimension(8), optional, target, intent(inout) :: time_intg_vec
      integer(int32), dimension(:), pointer :: time_intg_vec_p
      if (present(time_intg_vec)) then
         time_intg_vec_p => time_intg_vec
      else
         time_intg_vec_p => t_intg_vec_tmp
      endif
      call date_and_time(VALUES=time_intg_vec_p)
      call time_convert_intg_vec_to_sec(time_intg_vec_p,time_sec)
   end subroutine time_compute
   subroutine time_compute_difference( time_cpu_dif, time_cpu_ini, time_cpu_fin, time_cpu_tot )
      real(dp),                   intent(inout) :: time_cpu_dif
      real(dp), optional, target, intent(inout) :: time_cpu_ini
      real(dp), optional, target, intent(inout) :: time_cpu_fin
      real(dp), optional, target, intent(inout) :: time_cpu_tot
      real(dp) :: time_cpu_tmp
      if (.not.present(time_cpu_fin)) then
         call date_and_time(VALUES=t_intg_vec_tmp)
         call time_convert_intg_vec_to_sec(t_intg_vec_tmp,time_cpu_tmp)
      endif
      if (present(time_cpu_ini)) then
         if (present(time_cpu_fin)) then
            time_cpu_dif = time_cpu_fin - time_cpu_ini
         else
            time_cpu_dif = time_cpu_tmp - time_cpu_ini
         endif
      else
         if (present(time_cpu_fin)) then
            time_cpu_dif = time_cpu_fin - time_cpu_dif
         else
            time_cpu_dif = time_cpu_tmp - time_cpu_dif
         endif
      endif
      if (present(time_cpu_tot)) time_cpu_tot = time_cpu_tot + time_cpu_dif
   end subroutine time_compute_difference
   subroutine time_convert_intg_vec_to_sec( time_intg_vec, time_sec )
      integer(int32), dimension(8), intent(in)    :: time_intg_vec
      real(dp),                     intent(inout) :: time_sec
      time_sec = (dble(time_intg_vec(3))*24.0d0+dble(time_intg_vec(5)))*3600.0d0 + &
      & dble(time_intg_vec(6))*60.0d0 + dble(time_intg_vec(7)) + 0.001d0*dble(time_intg_vec(8))
   end subroutine time_convert_intg_vec_to_sec
   subroutine time_compute_dateandtime( time_date, time_time )
      character(8),  optional, target, intent(inout) :: time_date
      character(10), optional, target, intent(inout) :: time_time
      character(8),  pointer :: time_date_p
      character(10), pointer :: time_time_p
      if (present(time_date)) then
         time_date_p => time_date
      else
         time_date_p => t_date
      endif
      if (present(time_time)) then
         time_time_p => time_time
      else
         time_time_p => t_time
      endif
      call date_and_time(DATE=time_date_p,TIME=time_time_p)
   end subroutine time_compute_dateandtime
   subroutine time_compute_cputime( time_cpu )
      real(dp), optional, target, intent(inout) :: time_cpu
      real(dp),  pointer :: time_cpu_p
      if (present(time_cpu)) then
         time_cpu_p => time_cpu
      else
         time_cpu_p => t_cpu_tmp
      endif
      call cpu_time( time_cpu_p )
   end subroutine time_compute_cputime
   subroutine time_compute_cputime_difference( time_cpu_dif, time_cpu_ini, time_cpu_fin, time_cpu_tot )
      real(dp),                   intent(inout) :: time_cpu_dif
      real(dp), optional, target, intent(inout) :: time_cpu_ini
      real(dp), optional, target, intent(inout) :: time_cpu_fin
      real(dp), optional, target, intent(inout) :: time_cpu_tot
      real(dp) :: time_cpu_tmp
      if (.not.present(time_cpu_fin)) call cpu_time( time_cpu_tmp )
      if (present(time_cpu_ini)) then
         if (present(time_cpu_fin)) then
            time_cpu_dif = time_cpu_fin - time_cpu_ini
         else
            time_cpu_dif = time_cpu_tmp - time_cpu_ini
         endif
      else
         if (present(time_cpu_fin)) then
            time_cpu_dif = time_cpu_fin - time_cpu_dif
         else
            time_cpu_dif = time_cpu_tmp - time_cpu_dif
         endif
      endif
      if (present(time_cpu_tot)) time_cpu_tot = time_cpu_tot + time_cpu_dif
   end subroutine time_compute_cputime_difference
   subroutine init_simulation_timings( )
      character(51) :: full_string = ''
      t_date = ''
      t_time = ''
      t_tot = 0.0_dp
      t_cpu_tmp = 0.0_dp
      t_intg_vec_tmp = 0_int32
      t_lnalg_step = 0.0_dp
      t_lnalg_tot = 0.0_dp
      t_crs_tot = 0.0_dp
      t_opt_tot = 0.0_dp
      t_therm_vmc = 0.0_dp
      t_sampl_tot = 0.0_dp
      t_resampl_tot = 0.0_dp
      t_avg_tot = 0.0_dp
      t_therm_dmc = 0.0_dp
      t_dmc_tot = 0.0_dp
      t_brnch_tot = 0.0_dp
      call time_compute_dateandtime()
      call write_empty_line(stdout,0,mpi_rank)
      write(full_string,'("Starting the simulation on the ",A8, " at ",A8)') t_date(7:8)//'.'//t_date(5:6)//'.'//t_date(1:4), t_time(1:2)//':'//t_time(3:4)//':'//t_time(5:6)
      call write_simple_line(stdout,0,mpi_rank,2,"l",full_string)
      call time_compute(t_tot)
   end subroutine init_simulation_timings
   subroutine fnlz_simulation_timings( )
      character(51) :: full_string = ''
      call time_compute_difference(t_tot)
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","SIMULATION TIMINGS")
      call time_compute_dateandtime()
      call write_empty_line(stdout,0,mpi_rank)
      write(full_string,'("Ending the simulation on the ",A8, " at ",A8)') t_date(7:8)//'.'//t_date(5:6)//'.'//t_date(1:4), t_time(1:2)//':'//t_time(3:4)//':'//t_time(5:6)
      call write_simple_line(stdout,0,mpi_rank,2,"l",full_string)
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: VMC thermalization",t_therm_vmc,units='sec.' )
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: DMC thermalization",t_therm_dmc,units='sec.' )
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: Parameter solver",t_lnalg_tot,units='sec.' )
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: Optimization steps",t_opt_tot,units='sec.' )
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: Optimization steps with corr. samp.",t_crs_tot,units='sec.' )
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: Monte Carlo sampling",t_sampl_tot,units='sec.' )
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: Monte Carlo resampling",t_resampl_tot,units='sec.' )
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: Observable averages",t_avg_tot,units='sec.' )
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: DMC sampling",t_dmc_tot,units='sec.' )
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: DMC branching",t_brnch_tot,units='sec.' )
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"MPI AVG: Total simulation time",t_tot,units='sec.' )
   end subroutine fnlz_simulation_timings
end module timings_m
