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
module monte_carlo_checkpointfiles_m
   use fortran_kinds_v,             only: dp, int32, stdout
   use openmp_mpi_m,                only: mpi_rank 
   use write_lines_m
   use mersenne_twister19937_m,     only: load_mt19937, save_mt19937
   use quantum_monte_carlo_v,       only: bin_l, n_wlk_max, save_status, load_qmc_var, save_qmc_var
   use molecular_system_v,          only: molsys_prs
   use molecular_system_m,          only: load_fermionic_confs, save_fermionic_confs
   use qdo_system_v,                only: qdosys_prs 
   use qdo_system_m,                only: load_drudonic_confs, save_drudonic_confs
   use monte_carlo_averages_m,      only: load_qmcsmp, save_qmcsmp
   use diffusion_monte_carlo_v,     only: load_dmcmthd_var, save_dmcmthd_var
   implicit none
   logical, save, public :: checkpointdir_present
   public :: load_checkpointfiles, save_checkpointfiles
contains
   subroutine load_checkpointfiles()
      checkpointdir_present = .false. 
      inquire (file='./.qmc.save/qmc_variables.sav', EXIST=checkpointdir_present )
      if ( checkpointdir_present ) then
         call write_separator_line(stdout,0,mpi_rank,2,"_")
         call write_simple_line(stdout,0,mpi_rank,2,"c","LOADING CHECKPOINT FILES")
         call load_qmc_var( )
         call write_empty_line(stdout,0,mpi_rank)
         call load_mt19937( n_wlk_max )
         if (molsys_prs) call load_fermionic_confs( )
         if (qdosys_prs) call load_drudonic_confs( )
         call load_qmcsmp( )
         if ( save_status.ge.3_int32 ) call load_dmcmthd_var( n_wlk_max, bin_l )
      else 
         call write_empty_line(stdout,0,mpi_rank)
         call write_simple_line(stdout,0,mpi_rank,2,"l","No previous checkpoint directory found.")
         save_status = 0_int32
      endif
   end subroutine load_checkpointfiles
   subroutine save_checkpointfiles( save_status_in, step_in, create_backup )
      integer(int32), intent(in) :: save_status_in, step_in 
      logical,        intent(in) :: create_backup
      logical :: temp_directory
      temp_directory = .false.
      save_status = save_status_in
      checkpointdir_present = .false. 
      inquire (file='./qmc.save/qmc_variables.sav', EXIST=checkpointdir_present )
      if ( checkpointdir_present ) then
         if (create_backup) then
            if (step_in.eq.0_int32) call write_empty_line(stdout,0,mpi_rank)
            call write_line_check(stdout,0,mpi_rank,2,"Creating backup of Checkpoint directory ./qmc.save")
            inquire (file='./.qmc.save/qmc_variables.sav', EXIST=temp_directory) 
            if (temp_directory) then 
               if (mpi_rank.eq.0_int32) call system('rm -rf .qmc.save')
            endif 
            if (mpi_rank.eq.0_int32) call system('cp -rf qmc.save .qmc.save')
            if (step_in.eq.0_int32) call write_done(stdout,0,mpi_rank)
         endif
      else 
         if (step_in.eq.0_int32) call write_empty_line(stdout,0,mpi_rank)
         if (step_in.eq.0_int32) call write_line_check(stdout,0,mpi_rank,2,"Creating Checkpoint directory ./qmc.save")
         if (mpi_rank.eq.0_int32) call system('mkdir -p qmc.save')
         if (step_in.eq.0_int32) call write_done(stdout,0,mpi_rank)
      endif 
      if (step_in.eq.0_int32) then
         call write_line_check(stdout,0,mpi_rank,2,"Saving Checkpoint files") 
      endif
      call save_qmc_var(step_in )
      call save_mt19937( n_wlk_max )
      if (molsys_prs) call save_fermionic_confs( )
      if (qdosys_prs) call save_drudonic_confs( )
      call save_qmcsmp( )
      if ( save_status.ge.3_int32 ) call save_dmcmthd_var( n_wlk_max, bin_l )
      if (step_in.eq.0_int32) call write_done(stdout,0,mpi_rank)
   end subroutine save_checkpointfiles
   subroutine clean_checkpointfiles( )
      logical :: temp_directory
      inquire (file='./.qmc.save/qmc_variables.sav', EXIST=temp_directory) 
      if (temp_directory) then 
         call write_empty_line(stdout,0,mpi_rank)
         call write_line_check(stdout,0,mpi_rank,2,"Deleating backup Checkpoint directory ./.qmc.save")
         if (mpi_rank.eq.0_int32) call system('rm -rf .qmc.save')
         call write_done(stdout,0,mpi_rank)
      endif 
   end subroutine clean_checkpointfiles
end module monte_carlo_checkpointfiles_m
