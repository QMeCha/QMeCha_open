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
module argpio_mod
use fortran_kinds_v, only: dp, int32, stdout
use openmp_mpi_m
use write_lines_m
implicit none
integer(int32) :: n_args
character(len=20) :: args
contains 
subroutine read_comm_args( ierr ) 
#if defined _MPI || defined _MPIh || defined _MPI08
  use datpio_mod, only: read_input, bcast_input
#else
  use datpio_mod, only: read_input
#endif
  integer(int32), intent(inout) :: ierr
  if( mpi_rank.eq.0_int32 ) then
    call write_empty_line(stdout,0,mpi_rank)
    write(*,'(3X,"================================================================")')
    write(*,'(3X,"                    READING INPUT ARGUMENTS                     ")')
    call write_empty_line(stdout,0,mpi_rank)
    n_args = command_argument_count() 
    call get_command_argument(1_int32, args)
    select case ( trim(args) )
    case('-i')
      ierr = 0
      call get_command_argument(2, args)
      write(*,'(3X,"Reading input file : ",A20)') trim(args) 
      call read_input( trim(args), ierr )
    case default
      write(*,'(3X,"ERROR input takes flags : -i opt")')
      write(*,'(3X,"      -i opt run opt input file")')
      ierr = 1
    end select ! args
  endif ! mpi_rank.eq.0_int32
#if defined _MPI || defined _MPIh || defined _MPI08
  call bcast_input()
#endif
end subroutine read_comm_args
end module argpio_mod