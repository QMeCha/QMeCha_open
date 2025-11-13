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
module write_lines_m
   use fortran_kinds_v, only: int32,sp,dp
   implicit none
   integer(int32), private, parameter :: line_length = 70_int32
   interface write_variable_line
      procedure write_variable_intg_line
      procedure write_variable_dble_line
      procedure write_variable_strg_line
      procedure write_variable_logc_line
   end interface write_variable_line
   public  :: write_separator_line, write_simple_line, write_empty_line, write_variable_line,&
   & write_variable_error_dble_line
   private :: write_variable_intg_line, write_variable_dble_line, write_variable_strg_line, &
   & write_variable_logc_line
contains
   subroutine write_separator_line(writing_unit,writing_mpi_rank,mpi_rank,init_space,repeated_char)
      integer(int32), intent(in) :: writing_unit
      integer(int32), intent(in) :: writing_mpi_rank
      integer(int32), intent(in) :: mpi_rank
      integer(int32), intent(in) :: init_space
      character(*),   intent(in) :: repeated_char
      integer(int32) :: char_length, n_repetitions, n_remaining_char, i1
      char_length = len(repeated_char)
      n_repetitions = int(line_length / char_length)
      n_remaining_char = mod( line_length, char_length)
      if ( writing_mpi_rank >= 0_int32) then
         if ( mpi_rank.eq.writing_mpi_rank) &
         & write(writing_unit,'(100A)') (" ", i1 = 1_int32, init_space), (trim(repeated_char), i1 = 1_int32, n_repetitions),&
         & trim(repeated_char(1:n_remaining_char))
      else
         write(writing_unit,'(100A)') (" ", i1 = 1_int32, init_space), (trim(repeated_char), i1 = 1_int32, n_repetitions),&
         & trim(repeated_char(1:n_remaining_char))
      endif
   end subroutine write_separator_line
   subroutine write_simple_line(writing_unit,writing_mpi_rank,mpi_rank,init_space,title_alignment,title)
      integer(int32),   intent(in) :: writing_unit
      integer(int32),   intent(in) :: writing_mpi_rank
      integer(int32),   intent(in) :: mpi_rank
      integer(int32),   intent(in) :: init_space
      character(len=1), intent(in) :: title_alignment
      character(*),     intent(in) :: title
      integer(int32) :: n_spaces, i1
      select case(title_alignment)
       case('c')
         n_spaces = int((line_length - len(title))/2_int32)
       case('l')
         n_spaces = 0_int32
       case('r')
         n_spaces = line_length - len(title)
      end select
      if ( writing_mpi_rank >= 0_int32) then
         if ( mpi_rank.eq.writing_mpi_rank) &
         & write(writing_unit,'(100A)') (" ", i1 = 1_int32, init_space), (" ", i1 = 1_int32, n_spaces), trim(title)
      else
         write(writing_unit,'(100A)') (" ", i1 = 1_int32, init_space), (" ", i1 = 1_int32, n_spaces), trim(title)
      endif
   end subroutine write_simple_line
   subroutine write_line_check(writing_unit,writing_mpi_rank,mpi_rank,init_space,text)
      integer(int32),   intent(in) :: writing_unit
      integer(int32),   intent(in) :: writing_mpi_rank
      integer(int32),   intent(in) :: mpi_rank
      integer(int32),   intent(in) :: init_space
      character(*),     intent(in) :: text
      integer(int32) :: n_spaces, i1
      integer(int32) :: string_length 
      string_length = line_length - 6_int32
      if (len(text) < string_length ) then
         n_spaces = string_length - len(text)
         if ( writing_mpi_rank >= 0_int32) then
            if ( mpi_rank.eq.writing_mpi_rank) &
            & write(writing_unit,'(100A)',advance='no') (" ", i1 = 1_int32, init_space), trim(text), (" ", i1 = 1_int32, n_spaces)
         else
            write(writing_unit,'(100A)',advance='no') (" ", i1 = 1_int32, init_space), trim(text), (" ", i1 = 1_int32, n_spaces)
         endif
      else  
         if ( writing_mpi_rank >= 0_int32) then
            if ( mpi_rank.eq.writing_mpi_rank) &
            & write(writing_unit,'(100A)',advance='no') (" ", i1 = 1_int32, init_space), trim(text(1:string_length))
         else
            write(writing_unit,'(100A)',advance='no') (" ", i1 = 1_int32, init_space), trim(text(1:string_length))
         endif
      endif
   end subroutine write_line_check
   subroutine write_done(writing_unit,writing_mpi_rank,mpi_rank)
      integer(int32),   intent(in) :: writing_unit
      integer(int32),   intent(in) :: writing_mpi_rank
      integer(int32),   intent(in) :: mpi_rank
         if ( writing_mpi_rank >= 0_int32) then
            if ( mpi_rank.eq.writing_mpi_rank) &
            & write(writing_unit,'(6A)') "(done)"
         else
            write(writing_unit,'(6A)') "(done)"
         endif
   end subroutine write_done
   subroutine write_empty_line(writing_unit,writing_mpi_rank,mpi_rank)
      integer(int32), intent(in) :: writing_unit
      integer(int32), intent(in) :: writing_mpi_rank
      integer(int32), intent(in) :: mpi_rank
      if ( writing_mpi_rank >= 0_int32) then
         if ( mpi_rank.eq.writing_mpi_rank) write(writing_unit,*)
      else
         write(writing_unit,*)
      endif
   end subroutine write_empty_line
   subroutine write_variable_intg_line(writing_unit,writing_mpi_rank,mpi_rank,init_space,text,variable,units,var_name)
      integer(int32), intent(in) :: writing_unit
      integer(int32), intent(in) :: writing_mpi_rank
      integer(int32), intent(in) :: mpi_rank
      integer(int32), intent(in) :: init_space
      character(*),   intent(in) :: text
      integer(int32), intent(in) :: variable
      character(*),   optional, intent(in) :: units
      character(*),   optional, intent(in) :: var_name
      integer(int32) :: text_length, units_length, max_text_length, n_spaces, var_name_length, i1
      character(len=13)  :: num_string
      character(len=100)  :: full_string
      text_length = len(text)
      max_text_length = line_length - 16
      if ( present(units) ) then
         units_length = len(units)
         max_text_length = max_text_length - units_length - 2
      else
         units_length = 0_int32
      endif
      if ( present(var_name) ) then
         var_name_length = len(var_name)
         max_text_length = max_text_length - var_name_length - 2
      else
         var_name_length = 0_int32
      endif
      if ( max_text_length .gt. text_length) then
         n_spaces = max_text_length - text_length
      else
         text_length = max_text_length
         n_spaces = 0_int32
      endif
      write(num_string,'(I13)') variable
      if ( units_length > 0_int32) then
         if ( var_name_length > 0_int32) then
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"(",var_name,")[",units,"] : ",num_string
         else
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"[",units,"] : ",num_string
         endif
      else
         if ( var_name_length > 0_int32) then
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"(",var_name,") : ",num_string
         else
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces)," : ",num_string
         endif
      endif
      if ( writing_mpi_rank >= 0_int32) then
         if ( mpi_rank.eq.writing_mpi_rank) write(writing_unit,'(100A)') trim(full_string)
      else
         write(writing_unit,'(100A)') trim(full_string)
      endif
   end subroutine write_variable_intg_line
   subroutine write_variable_dble_line(writing_unit,writing_mpi_rank,mpi_rank,init_space,text,variable,units,var_name)
      integer(int32), intent(in) :: writing_unit
      integer(int32), intent(in) :: writing_mpi_rank
      integer(int32), intent(in) :: mpi_rank
      integer(int32), intent(in) :: init_space
      character(*),   intent(in) :: text
      real(dp),       intent(in) :: variable
      character(*),   optional, intent(in) :: units
      character(*),   optional, intent(in) :: var_name
      integer(int32) :: text_length, units_length, max_text_length, n_spaces, var_name_length, i1
      character(len=13)  :: num_string
      character(len=100)  :: full_string
      text_length = len(text)
      max_text_length = line_length - 16
      if ( present(units) ) then
         units_length = len(units)
         max_text_length = max_text_length - units_length - 2
      else
         units_length = 0_int32
      endif
      if ( present(var_name) ) then
         var_name_length = len(var_name)
         max_text_length = max_text_length - var_name_length - 2
      else
         var_name_length = 0_int32
      endif
      if ( max_text_length .gt. text_length) then
         n_spaces = max_text_length - text_length
      else
         text_length = max_text_length
         n_spaces = 0_int32
      endif
      write(num_string,'(E13.6)') variable
      if ( units_length > 0_int32) then
         if ( var_name_length > 0_int32) then
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"(",var_name,")[",units,"] : ",num_string
         else
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"[",units,"] : ",num_string
         endif
      else
         if ( var_name_length > 0_int32) then
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"(",var_name,") : ",num_string
         else
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces)," : ",num_string
         endif
      endif
      if ( writing_mpi_rank >= 0_int32) then
         if ( mpi_rank.eq.writing_mpi_rank) write(writing_unit,'(100A)') trim(full_string)
      else
         write(writing_unit,'(100A)') trim(full_string)
      endif
   end subroutine write_variable_dble_line
   subroutine write_variable_error_dble_line(writing_unit,writing_mpi_rank,mpi_rank,init_space,text,variable,variable_err,units)
      integer(int32), intent(in) :: writing_unit
      integer(int32), intent(in) :: writing_mpi_rank
      integer(int32), intent(in) :: mpi_rank
      integer(int32), intent(in) :: init_space
      character(*),   intent(in) :: text
      real(dp),       intent(in) :: variable
      real(dp),       intent(in) :: variable_err
      character(*),   optional, intent(in) :: units
      integer(int32) :: text_length, units_length, max_text_length, n_spaces, i1
      character(len=46)  :: num_string
      character(len=100)  :: full_string
      text_length = len(text)
      max_text_length = line_length - 47 - init_space
      if ( present(units) ) then
         units_length = len(units)
         max_text_length = max_text_length - units_length - 2
      else
         units_length = 0_int32
      endif
      if ( max_text_length .gt. text_length) then
         n_spaces = max_text_length - text_length
      else
         text_length = max_text_length
         n_spaces = 0_int32
      endif
      write(num_string,'(E21.14," +/-",E21.14)') variable, variable_err
      if ( units_length > 0_int32) then
         write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"[",units,"] : ",num_string
      else
         write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces)," : ",num_string
      endif
      if ( writing_mpi_rank >= 0_int32) then
         if ( mpi_rank.eq.writing_mpi_rank) write(writing_unit,'(100A)') trim(full_string)
      else
         write(writing_unit,'(100A)') trim(full_string)
      endif
   end subroutine write_variable_error_dble_line
   subroutine write_variable_strg_line(writing_unit,writing_mpi_rank,mpi_rank,init_space,text,variable,units,var_name)
      integer(int32), intent(in) :: writing_unit
      integer(int32), intent(in) :: writing_mpi_rank
      integer(int32), intent(in) :: mpi_rank
      integer(int32), intent(in) :: init_space
      character(*),   intent(in) :: text
      character(*),   intent(in) :: variable
      character(*),   optional, intent(in) :: units
      character(*),   optional, intent(in) :: var_name
      integer(int32) :: text_length, units_length, max_text_length, n_spaces, var_name_length, i1
      character(len=100)  :: full_string
      text_length = len(text)
      max_text_length = line_length - 16
      if ( present(units) ) then
         units_length = len(units)
         max_text_length = max_text_length - units_length - 2
      else
         units_length = 0_int32
      endif
      if ( present(var_name) ) then
         var_name_length = len(var_name)
         max_text_length = max_text_length - var_name_length - 2
      else
         var_name_length = 0_int32
      endif
      if ( max_text_length .gt. text_length) then
         n_spaces = max_text_length - text_length
      else
         text_length = max_text_length
         n_spaces = 0_int32
      endif
      if ( units_length > 0_int32) then
         if ( var_name_length > 0_int32) then
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"(",var_name,")[",units,"] : ",trim(variable)
         else
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"[",units,"] : ",trim(variable)
         endif
      else
         if ( var_name_length > 0_int32) then
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"(",var_name,") : ",trim(variable)
         else
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces)," : ",trim(variable)
         endif
      endif
      if ( writing_mpi_rank >= 0_int32) then
         if ( mpi_rank.eq.writing_mpi_rank) write(writing_unit,'(100A)') trim(full_string)
      else
         write(writing_unit,'(100A)') trim(full_string)
      endif
   end subroutine write_variable_strg_line
   subroutine write_variable_logc_line(writing_unit,writing_mpi_rank,mpi_rank,init_space,text,variable,units,var_name)
      integer(int32), intent(in) :: writing_unit
      integer(int32), intent(in) :: writing_mpi_rank
      integer(int32), intent(in) :: mpi_rank
      integer(int32), intent(in) :: init_space
      character(*),   intent(in) :: text
      logical,        intent(in) :: variable
      character(*),   optional, intent(in) :: units
      character(*),   optional, intent(in) :: var_name
      integer(int32) :: text_length, units_length, max_text_length, n_spaces, var_name_length, i1
      character(len=13)  :: num_string = ""
      character(len=100)  :: full_string = ""
      text_length = len(text)
      max_text_length = line_length - 16
      if ( present(units) ) then
         units_length = len(units)
         max_text_length = max_text_length - units_length - 2
      else
         units_length = 0_int32
      endif
      if ( present(var_name) ) then
         var_name_length = len(var_name)
         max_text_length = max_text_length - var_name_length - 2
      else
         var_name_length = 0_int32
      endif
      if ( max_text_length .gt. text_length) then
         n_spaces = max_text_length - text_length
      else
         text_length = max_text_length
         n_spaces = 0_int32
      endif
      write(num_string,'(12X,L1)') variable
      if ( units_length > 0_int32) then
         if ( var_name_length > 0_int32) then
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"(",var_name,")[",units,"] : ",num_string
         else
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"[",units,"] : ",num_string
         endif
      else
         if ( var_name_length > 0_int32) then
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces),"(",var_name,") : ",num_string
         else
            write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces)," : ",num_string
         endif
      endif
      if ( writing_mpi_rank >= 0_int32) then
         if ( mpi_rank.eq.writing_mpi_rank) write(writing_unit,'(100A)') trim(full_string)
      else
         write(writing_unit,'(100A)') trim(full_string)
      endif
   end subroutine write_variable_logc_line
   subroutine write_simple_strn_line(writing_unit,writing_mpi_rank,mpi_rank,init_space,text,variable,max_text_length)
      integer(int32), intent(in) :: writing_unit
      integer(int32), intent(in) :: writing_mpi_rank
      integer(int32), intent(in) :: mpi_rank
      integer(int32), intent(in) :: init_space
      character(*),   intent(in) :: text
      character(*),   intent(in) :: variable
      integer(int32), intent(in) :: max_text_length
      integer(int32) :: text_length, n_spaces, i1
      character(len=100)  :: full_string
      text_length = len(text)
      if ( text_length >= max_text_length ) then
         text_length = max_text_length
         n_spaces = 0
      else
         n_spaces = max_text_length - text_length
      endif
      write(full_string,'(100A)') (" ", i1 = 1_int32, init_space),text(1:text_length),(" ",i1=1_int32, n_spaces)," : ", trim(variable)
      if ( writing_mpi_rank >= 0_int32) then
         if ( mpi_rank.eq.writing_mpi_rank) write(writing_unit,'(100A)') trim(full_string)
      else
         write(writing_unit,'(100A)') trim(full_string)
      endif
   end subroutine write_simple_strn_line
end module write_lines_m
