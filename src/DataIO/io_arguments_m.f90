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
module io_arguments_m
   use fortran_kinds_v,             only: dp, int32, stdout
   use openmp_mpi_m
   use write_lines_m,               only: write_separator_line, write_simple_line, &
   & write_empty_line, write_simple_strn_line, &
   & write_variable_line
   use quantum_monte_carlo_v,       only: qmc_mthd, sysname
   use wavefunction_optimization_v, only: opt_mthd
   use orbitalsbased_jastrow_factors_params_m, only: cpl_j, cpl_m, corr_type
   use molecular_system_v,          only: fle_coord
   use qdo_system_v,                only: fle_coord_qdo
   use basisset_v,                  only: fle_basis, fle_basis_qdo
   use pseudopotentials_v,          only: fle_pspot
   use io_datasheet_m,              only: fle_dataio, dflt_data_var, bcst_data_var,&
   & write_data_fle, read_data_fle
   implicit none
   public  :: read_command_args
   private :: read_input_variables, comm_h_arg, comm_i_arg, comm_pi_arg, comm_s_arg
contains
   subroutine read_command_args
      integer(int32)     :: n_args, ierr
      character(len=100) :: args
      call dflt_data_var()
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","READING INPUT ARGUMENTS")
      call write_empty_line(stdout,0,mpi_rank)
      ierr = 0_int32
      if (mpi_rank.eq.0) then
         n_args = command_argument_count()
         if (n_args.eq.0) then
            call comm_d_arg()
            ierr = -1
         else
            call get_command_argument(1, args)
            select case ( trim(args) )
             case('-h')
               call comm_h_arg( )
               ierr = -1_int32
             case('-pi')
               call get_command_argument(2, args)
               call comm_pi_arg( args )
               ierr = -1_int32
             case('-s')
               call comm_s_arg( )
               call read_input_variables( n_args, 2_int32, args, ierr )
             case('-i')
               call comm_i_arg( args )
               call read_input_variables( n_args, 3_int32, args, ierr )
               if (ierr.ge.0) call read_data_fle( fle_dataio, ierr )
             case default
               call write_simple_line(stdout,0,mpi_rank,5,"l","Error!!! Input flags are wrongly defined.")
               call comm_h_arg( )
               ierr = -1
            end select 
         endif
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
      if (ierr.lt.0) then
         call fnlz_ompmpi_env()
         stop
      else
#if defined _MPI || defined _MPIh || defined _MPI08
         call bcst_data_var()
#endif
      endif
   end subroutine read_command_args
   subroutine comm_pi_arg( fle_to_prnt )
      character(*), intent(in) :: fle_to_prnt
      call write_simple_line(stdout,0,mpi_rank,2,"c","____________________________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","QMeCha evoked printing of input file")
      call write_empty_line(stdout,0,mpi_rank)
      select case( trim(fle_to_prnt) )
       case('opt')
         if ( mpi_rank.eq.0_int32 ) call write_data_fle(dflt=.true., inp_type='opt', fle_name='opt.input')
       case('vmc')
         if ( mpi_rank.eq.0_int32 ) call write_data_fle(dflt=.true., inp_type='vmc', fle_name='vmc.input')
       case('dmc')
         if ( mpi_rank.eq.0_int32 ) call write_data_fle(dflt=.true., inp_type='dmc', fle_name='dmc.input')
       case default
         call write_simple_line(stdout,0,mpi_rank,2,"l","ERROR!!! -pi flag takes options : vmc, dmc, opt ")
         call write_empty_line(stdout,0,mpi_rank)
         call write_simple_line(stdout,0,mpi_rank,8,"l","vmc : Variational Monte Carlo run    ")
         call write_simple_line(stdout,0,mpi_rank,8,"l","dmc : Diffusion Monte Carlo run    ")
         call write_simple_line(stdout,0,mpi_rank,8,"l","opt : Wave function optimization run ")
      end select
   end subroutine comm_pi_arg
   subroutine comm_s_arg
      call write_simple_line(stdout,0,mpi_rank,2,"c","__________________________________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","QMeCha evoked wave function initialization")
      call write_empty_line(stdout,0,mpi_rank)
      qmc_mthd   = 'non' ; opt_mthd   = 'res'
   end subroutine comm_s_arg
   subroutine comm_i_arg( args )
      character(*),    intent(inout) :: args
      call write_simple_line(stdout,0,mpi_rank,2,"c","___________________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","QMeCha evoked job execution")
      call write_empty_line(stdout,0,mpi_rank)
      call get_command_argument(2, args) ; fle_dataio = trim(args)
      call write_variable_line(stdout,0,mpi_rank,2,"Input file", fle_dataio)
   end subroutine comm_i_arg
   subroutine comm_h_arg( )
      call write_simple_line(stdout,0,mpi_rank,2,"c","___________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","Help command evoked")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"l","Specify the following variables :")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,5,"l","NUM_THREADS=1 (Number of threads per mpi task)")
      call write_simple_line(stdout,0,mpi_rank,5,"l","NUM_SOCKETS=4 (Number of mpi task)")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,5,"l","QMECHA_BIN = $DIRECTORY (your bin directory)")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,5,"l","export OMP_NUM_THREADS = $NUM_THREADS (number of cores/socket)")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"l","Run the code as follows :")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,5,"l","mpiexec -np $NUM_SOCKETS ${QMECHA_BIN}/QMeCha -i input > out &")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"l","QMeCha accepts the optional additional flags in the input line that are")
      call write_simple_line(stdout,0,mpi_rank,2,"l","eventually overwritten by the input variables:")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,5,"l","-b   ${basis_set_file}")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-m   ${molecular_input_file}")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-p   ${pseudopotential_file}")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-j   ${dynamical_jastrow_coupling_type} (ma,mf,md,ca,cf,cd)")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-bq  ${basis_set_file_for_qdos}")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-mq  ${molecular_input_file_for_qdos}")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-sn  (optional) ${molecular system name}")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"l","Please read QMeCha manual for further details.")
   end subroutine comm_h_arg
   subroutine comm_d_arg()
      call write_simple_line(stdout,0,mpi_rank,2,"c","___________________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","QMeCha evoked without flags")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"l","The QMeCha code expects the following initial flags:")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,3,"l","-i   ${qmecha_input_file}")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,3,"l","-s   the wave function initialization is evoked. ")
      call write_simple_line(stdout,0,mpi_rank,3,"l","     This command must be followed by the flags: ")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,5,"l","-b   ${basis_set_file}")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-m   ${molecular_input_file}")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-bq  (optional) ${basis_set_file_for_qdos}")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-mq  (optional) ${molecular_input_file_for_qdos}")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-p   (optional) ${pseudopotential_file}")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-j   (optional) ${dynamical_jastrow_coupling_type} (ma,mf,md,ca,cf,cd)")
      call write_simple_line(stdout,0,mpi_rank,5,"l","-sn  (optional) ${molecular system name}")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,3,"l","-pi  ${input_type} prints example of input file (vmc,dmc,opt)")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,3,"l","-h   call help")
      call write_empty_line(stdout,0,mpi_rank)
      call write_simple_line(stdout,0,mpi_rank,2,"l","Please read QMeCha manual for further details.")
   end subroutine comm_d_arg
   subroutine read_input_variables( n_args, first_arg, args, ierr )
      integer(int32), intent(in)    :: n_args
      integer(int32), intent(in)    :: first_arg
      character(*),   intent(inout) :: args
      integer(int32), intent(inout) :: ierr
      character(3)   :: j_tmp
      integer(int32) :: i1
      do i1 = first_arg, n_args, 2_int32
         call get_command_argument(i1, args)
         select case ( trim(args) )
          case('-m')
            call get_command_argument(i1+1, args)
            fle_coord = trim(args)
            call write_variable_line(stdout,0,mpi_rank,2,"Molecular system file", fle_coord)
          case('-b')
            call get_command_argument(i1+1, args)
            fle_basis = trim(args)
            call write_variable_line(stdout,0,mpi_rank,2,"Fermionic basis set file", fle_basis)
          case('-p')
            call get_command_argument(i1+1, args)
            fle_pspot = trim(args)
            call write_variable_line(stdout,0,mpi_rank,2,"Pseudopotential file", fle_pspot)
          case('-mq')
            call get_command_argument(i1+1, args)
            fle_coord_qdo = trim(args)
            call write_variable_line(stdout,0,mpi_rank,2,"Drudonic molecular file", fle_coord_qdo)
          case('-bq')
            call get_command_argument(i1+1, args)
            fle_basis_qdo = trim(args)
            call write_variable_line(stdout,0,mpi_rank,2,"Drudonic basis set file", fle_basis_qdo)
          case('-j')
            call get_command_argument(i1+1, args)
            j_tmp = trim(args)
            cpl_m = 'm'
            cpl_j = j_tmp(1:1)
            corr_type = 0_int32
            call write_variable_line(stdout,0,mpi_rank,2,"Jastrow coupling type", cpl_j)
         case('-sn')
           call get_command_argument(i1+1, args)
           sysname = trim(args)
           call write_variable_line(stdout,0,mpi_rank,2,"System name used for save files", sysname)
          case default
            call write_variable_line(stdout,0,mpi_rank,2,"Error!!! Flag not recongnized", trim(args))
            ierr = -1
            exit
         end select
      enddo
   end subroutine read_input_variables
end module io_arguments_m
