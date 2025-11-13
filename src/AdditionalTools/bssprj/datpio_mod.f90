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
module datpio_mod
use fortran_kinds_v, only: int32, dp
use openmp_mpi_m
use molecular_system_v, only: fle_coord
use intgrd_mod, only: n_int_points, orb_cut_off, nuc_cut_off
implicit none
character(len=20), save, public :: trg_basis
character(len=20), save, public :: new_basis
real(dp),       save, public :: b_par, l_par
logical,        save, public :: oe_opt, jc1_opt
namelist /optmet/ n_int_points, orb_cut_off, oe_opt, jc1_opt, b_par, l_par, nuc_cut_off
namelist /basset/ trg_basis, new_basis
namelist /molsys/ fle_coord
contains
subroutine read_input( fle_name, ierr )
  character(*),   intent(in)    :: fle_name
  integer(int32), intent(inout) :: ierr
  oe_opt = .false.
  jc1_opt = .false.
  b_par = 0.0_dp
  l_par = 1.0_dp
  orb_cut_off = 0.1d-7
  nuc_cut_off = 0.1d-1
  if( mpi_rank.eq.0_int32 ) then
    call write_empty_line(stdout,0,mpi_rank)
    write(*,'(3X,"================================================================")')
    write(*,'(3X,"                        READING DATAFILE                        ")')
    call write_empty_line(stdout,0,mpi_rank)
    open(unit=1, file=fle_name, action='read', form='formatted', status='old', iostat=ierr )
    if(ierr.ne.0_int32) then
      write(*,'(3X,"ERROR cannot open input file ")')
    endif
    if( ierr.eq.0_int32 ) then
      read(1,nml=optmet, iostat=ierr )
      if ( ierr .lt. 0_int32 ) then
        write(*,'(3X,"ERROR cannot find section %optmet ")')
      elseif ( ierr .gt. 0_int32 ) then
        write(*,'(3X,"ERROR while reading section %optmet ")')
      endif
      rewind(1)
    endif
    if( ierr.eq.0_int32 ) then
      read(1,nml=molsys, iostat=ierr )
      if ( ierr .lt. 0_int32 ) then
        write(*,'(3X,"ERROR cannot find section %molsys ")')
      elseif ( ierr .gt. 0_int32 ) then
        write(*,'(3X,"ERROR while reading section %molsys ")')
      endif
      rewind(1)
    endif
    if( ierr.eq.0_int32 ) then
      read(1,nml=basset, iostat=ierr )
      if ( ierr .lt. 0_int32 ) then
        write(*,'(3X,"ERROR cannot find section %basset ")')
      elseif ( ierr .gt. 0_int32 ) then
        write(*,'(3X,"ERROR while reading section %basset ")')
      endif
    endif
    close(1)
  endif ! ( mpi_rank.eq.0_int32 )
end subroutine read_input
#if defined _MPI || defined _MPIh || defined _MPI08
subroutine bcast_input(  ) 
  call mpi_bcast(orb_cut_off , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
  call mpi_bcast(nuc_cut_off , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
  call mpi_bcast(n_int_points, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr) 
  call mpi_bcast(oe_opt, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr) 
  call mpi_bcast(jc1_opt, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr) 
  call mpi_bcast(b_par, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr) 
  call mpi_bcast(l_par, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr) 
end subroutine bcast_input
#endif
end module datpio_mod