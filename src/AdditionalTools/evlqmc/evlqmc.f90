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
program evlqmc
  use fortran_kinds_v, only: int32
  use quantum_monte_carlo_v, only: qmc_mthd
  use molecular_system_m, only: read_coord_fle, init_molec_sys
  use molecular_system_v, only: fle_coord
  use openmp_mpi_m, only: init_ompmpi_env, fnlz_ompmpi_env
  use monte_carlo_averages_m
  implicit none
  integer(int32)    :: n_args
  character(len=32) :: arg
  integer(int32)    :: i1
  call init_ompmpi_env()
  call write_empty_line(stdout,0,mpi_rank)
  write(*,'(3X,"================================================================")')
  write(*,'(3X,"                    READING INPUT ARGUMENTS                     ")')
  call write_empty_line(stdout,0,mpi_rank)
  fle_coord = 'none'
  n_args = command_argument_count() 
  if (n_args.gt.0_int32) then
    do i1 = 1_int32, n_args, 2_int32
      call get_command_argument(i1, arg)
      select case ( trim(arg) )
      case('-m')
        call get_command_argument(i1+1, arg)
        fle_coord = trim(arg) 
        write(*,'(" Reading coordinate file : ",A20)') fle_coord
      case default
        write(*,'(3X,"ERROR input takes flag : -s name_of_smpl_file ")')
        write(*,'(3X,"      name_of_smpl_file : file containing VMC sampling")')
      end select ! args
    enddo
  else
    write(*,'(3X,"WARNING input is empty, will try to read the default file : vmc_smpl.dat")')
  endif
  !if (use_sym) then
  !  call ini_rndnmb()
  !  call ini_symop()
  !  call read_coord_fle( ierr )
  !  if ( ierr .ne. 0_int32 ) then
  !    call fnlz_ompmpi_env()
  !    stop
  !  endif
  !  call init_molec_sys()
  !endif
  call read_qmcsmp(qmc_mthd)
  call init_qmcsmp( )
  do i1 = 1_int32, n_tot_smpl
    read(10) obs_b(:,1)
    if(qmc_mthd.eq.'vmc') then
      call accu_vmcavg( obs_b(:,1) )
    else 
      call accu_dmcavg( obs_b(:,1) )
    endif
  enddo
  call write_qmcavg( )
  call fnlz_ompmpi_env()
end program evlqmc
