!===============================================================================
!> SUBROUTINE qmecha_version (QMeCha BETA VERSION 2023)
!>  
!> @Author Matteo Barborini (matteo.barborini@gmail.com)
!>  
!> @Private repository https://github.com/mbarborini/QMeCha_beta
!>
subroutine qmecha_version()

   use fortran_kinds_v, only: stdout
   use openmp_mpi_m,    only: mpi_rank
   use write_lines_m,   only: write_simple_strn_line, write_empty_line, write_simple_line
   implicit none 
   call write_empty_line(stdout,0,mpi_rank)
   call write_simple_line(stdout,0,mpi_rank,8,"l"," _______   _______         _______  _                    ")
   call write_simple_line(stdout,0,mpi_rank,8,"l","(_______) (_______)       (_______)| |                   ")
   call write_simple_line(stdout,0,mpi_rank,8,"l"," _    _    _  _  _  _____  _       | |__   _____         ")
   call write_simple_line(stdout,0,mpi_rank,8,"l","| |  | |  | ||_|| || ___ || |      |  _ \ (____ |        ")
   call write_simple_line(stdout,0,mpi_rank,8,"l","| |__| |_ | |   | || ____|| |_____ | | | |/ ___ |        ")
   call write_simple_line(stdout,0,mpi_rank,8,"l","\________)|_|   |_||_____)\_______)|_| |_|\_____|        ")
   call write_simple_line(stdout,0,mpi_rank,8,"l","                                                  (Beta) ")
   call write_empty_line(stdout,0,mpi_rank)
   call write_simple_line(stdout,0,mpi_rank,2,"l","@Author             : Matteo Barborini (matteo.barborini@gmail.com)")
   call write_empty_line(stdout,0,mpi_rank)
   call write_simple_line(stdout,0,mpi_rank,2,"l","@Private repository : https://github.com/mbarborini/QMeCha_beta    ")
   call write_empty_line(stdout,0,mpi_rank)
   call write_simple_strn_line(stdout,0,mpi_rank,2,"Git version","",25)
   call write_simple_strn_line(stdout,0,mpi_rank,2,"Compilation date","04-nov-2025",25)

end subroutine qmecha_version