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
subroutine check_namelist_optional( ierr, namelist_name )
   use fortran_kinds_v, only: int32, stdout
   use openmp_mpi_m,    only: mpi_rank
   use write_lines_m,   only: write_variable_line
   implicit none
   integer(int32), intent(inout) :: ierr
   character(*),   intent(in)    :: namelist_name
   if ( ierr .lt. 0_int32 ) then
      call write_variable_line(stdout,0,mpi_rank,2,"WARNING, Cannot find optional section", namelist_name)
      ierr = 0_int32
   elseif ( ierr .gt. 0_int32 ) then
      call write_variable_line(stdout,0,mpi_rank,2,"ERROR!!! While reading section", namelist_name)
   else
      call write_variable_line(stdout,0,mpi_rank,2,"Reading section", namelist_name)
   endif
end subroutine check_namelist_optional
subroutine check_namelist_mandatory( ierr, namelist_name )
   use fortran_kinds_v, only: int32, stdout
   use openmp_mpi_m,    only: mpi_rank
   use write_lines_m,   only: write_variable_line
   implicit none
   integer(int32), intent(in)    :: ierr
   character(*),   intent(in)    :: namelist_name
   if ( ierr .lt. 0_int32 ) then
      call write_variable_line(stdout,0,mpi_rank,2,"ERROR!!! Cannot find mandatory section", namelist_name)
   elseif ( ierr .gt. 0_int32 ) then
      call write_variable_line(stdout,0,mpi_rank,2,"ERROR!!! While reading section", namelist_name)
   else
      call write_variable_line(stdout,0,mpi_rank,2,"Reading section", namelist_name)
   endif
end subroutine check_namelist_mandatory
