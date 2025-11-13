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
module branching_m
   use fortran_kinds_v,         only: dp, int32
   use openmp_mpi_m
   use timings_m
   use write_lines_m,           only: write_simple_line
   use quantum_monte_carlo_v,   only: n_wlk, n_tot_wlk, n_wlk_max
   use diffusion_monte_carlo_v, only: w, restart_walker, w_chk, w_max, &
   & log_w_max
   use mersenne_twister19937_m, only: rndn
   use molecular_system_v,      only: n_fe, molsys_prs
   use fermionic_config_c,      only: frm_cnf
   use qdo_system_v,            only: n_qdo, qdosys_prs
   use drudonic_config_c,       only: drd_cnf
   implicit none
   integer(int32), save, private :: n_empty, n_repli
   integer(int32), save, private :: n_tot_empty, n_tot_repli, n_max_repli
   integer(int32), save, private, allocatable, dimension(:,:) :: empty_lst
   integer(int32), save, private, allocatable, dimension(:,:) :: repli_lst
#if defined _MPI || defined _MPIh || defined _MPI08
   real(dp),       save, private, allocatable, dimension(:)   :: w_shift
   integer(int32), save, private, allocatable, dimension(:,:) :: er_lst, er_indx
#endif
   public  :: Init_Branching, Exec_Branching, Fnlz_Branching
#if defined _MPI || defined _MPIh || defined _MPI08
   private :: Init_ReplicaTable, copy_walker, redistribute_InterNode, redistribute_IntraNode, VariableWalkersBranching,&
   & StochasticReconfigurationBranching
#else
   private :: Init_ReplicaTable, copy_walker, redistribute_IntraNode, VariableWalkersBranching, StochasticReconfigurationBranching
#endif
contains
   subroutine Init_Branching()
      n_max_repli = n_tot_wlk
      allocate( empty_lst(1:2,1:n_max_repli) )
      allocate( repli_lst(1:2,1:n_max_repli) )
#if defined _MPI || defined _MPIh || defined _MPI08
      allocate( er_lst(1:n_mpi_tasks,1:2) ) ; er_lst = 0_int32
      allocate( er_indx(1:n_mpi_tasks,1:2) ) ; er_lst = 0_int32
      allocate( w_shift(1:n_mpi_tasks) ) ; w_shift = 0.0_dp
#endif
      log_w_max = log(w_max)
   end subroutine Init_Branching
   subroutine Exec_Branching( brnch_type )
      integer(int32), intent(in)    :: brnch_type
      call time_compute( t_brnch_step )
      select case ( brnch_type )
       case(1)
         call VariableWalkersBranching( )
       case(2)
         call StochasticReconfigurationBranching( )
       case default
      end select
      call time_compute_difference(t_brnch_step,time_cpu_tot=t_brnch_tot)
   end subroutine Exec_Branching
   subroutine chk_brnchn( do_branching )
      logical, intent(inout) :: do_branching
      real(dp) :: w_single_max, w_single_min
      do_branching = .false.
      w_single_max = maxval(w(1:n_wlk))
      w_single_min = minval(w(1:n_wlk))
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE, w_single_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(MPI_IN_PLACE, w_single_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr )
#endif
      if (w_single_max.ge.(1.0_dp+w_chk) .or. w_single_min.le.(1.0_dp-w_chk)) do_branching = .true.
   end subroutine chk_brnchn
   subroutine VariableWalkersBranching
      integer(int32)  :: G_w_tot
      integer(int32)  :: n_wlk_a
      integer(int32)  :: i1
      do i1 = 1_int32, n_wlk
         w(i1) = w(i1) + rndn(i1)%rndc()
         if (w(i1).gt.w_max) w(i1) = w_max
         w(i1) = dble(int(w(i1)))
      enddo
      G_w_tot = int(sum(w(1:n_wlk)))
#if defined _MPI || defined _MPIh || defined _MPI08
      if (n_mpi_tasks.gt.1_int32) then
         call mpi_allreduce(MPI_IN_PLACE, G_w_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr )
         n_wlk_a = G_w_tot / n_mpi_tasks
         if ( mpi_rank.le. mod(G_w_tot, n_mpi_tasks)-1 )  n_wlk_a = n_wlk_a + 1_int32
      else
#endif
         n_wlk_a = G_w_tot
#if defined _MPI || defined _MPIh || defined _MPI08
      endif
#endif
      if (dble(G_w_tot)/ dble(n_tot_wlk) .lt.10.0d-3) then
         call write_simple_line(stdout,0,mpi_rank,2,"l","ERROR!!! Weight of the configurations is too low")
         call fnlz_simulation_timings()
         call fnlz_ompmpi_env()
         stop
      endif
      call Init_ReplicaTable( n_tot_wlk )
      if ( n_wlk.gt.n_wlk_a ) then
         i1 = n_wlk
         do while ( i1 .gt. n_wlk_a )
            if ( w(i1).ge.1.0_dp ) then
               n_repli = n_repli + 1_int32
               repli_lst(1,n_repli) = i1
               repli_lst(2,n_repli) = mpi_rank
            endif
            i1 = i1 - 1_int32
         enddo
      elseif ( n_wlk.lt.n_wlk_a ) then
         i1 = n_wlk + 1_int32
         do while ( i1 .le. n_wlk_a )
            n_empty = n_empty + 1_int32
            empty_lst(1,n_empty) = i1
            empty_lst(2,n_empty) = mpi_rank
            i1 = i1 + 1_int32
         enddo
      endif
      if ( n_repli.gt.0_int32 .and. n_empty.gt.0_int32 ) call redistribute_IntraNode( )
#if defined _MPI || defined _MPIh || defined _MPI08
      if ( n_mpi_tasks.gt.1_int32 ) call redistribute_InterNode(  )
#endif
      n_wlk = n_wlk_a
      w(1:n_wlk) = 1.0_dp ; w(n_wlk+1:n_wlk_max) = 0.0_dp
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(n_wlk, n_tot_wlk, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
   end subroutine VariableWalkersBranching
   subroutine StochasticReconfigurationBranching
      real(dp)        :: G_w_scale, G_w_tot
      real(dp)        :: z
      real(dp)        :: z_max
      integer(int32)  :: i1, i2
!$omp parallel
!$omp do
      do i1 = 1_int32, n_wlk
         if ( w(i1).gt.w_max ) w(i1) = w_max
      enddo
!$omp end do
!$omp end parallel
      G_w_tot = sum(w(1:n_wlk))
      G_w_scale = G_w_tot
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE, G_w_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      if (G_w_tot/ dble(n_tot_wlk) .lt.10.0d-3) then
         call write_simple_line(stdout,0,mpi_rank,2,"l","ERROR!!! Weight of the configurations is too low")
         call fnlz_simulation_timings()
         call fnlz_ompmpi_env()
         stop
      endif
      w = w / G_w_tot
      G_w_scale = G_w_scale / G_w_tot
#if defined _MPI || defined _MPIh || defined _MPI08
      w_shift = 0.0_dp
      call mpi_allgather( G_w_scale, 1, MPI_DOUBLE_PRECISION, w_shift, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpierr )
      do i1 = 2_int32, n_mpi_tasks
         w_shift(i1) = w_shift(i1) + w_shift(i1-1)
      enddo
      if (mpi_rank.gt.0_int32) w(1) = w(1) + w_shift(mpi_rank)
#endif
      do i1 = 2_int32, n_wlk
         w(i1) = w(i1) + w(i1-1)
      enddo
      if (mpi_rank.eq.0_int32) z = (1.0 - rndn(1)%rndo()) / dble(n_tot_wlk)
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast( z, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr )
#endif
      z_max = ( dble(n_tot_wlk) - 1.0_dp ) / dble(n_tot_wlk) + z
#if defined _MPI || defined _MPIh || defined _MPI08
      if (mpi_rank.gt.0_int32) then
         do while (z.lt.w_shift(mpi_rank))
            z = z + 1.0_dp / dble(n_tot_wlk)
         enddo
      endif
#endif
      i1 = 1_int32 ; i2 = 0_int32
      do while ( i1.le.n_wlk )
         if ( z .le. w(i1) .and. z .le. z_max ) then
            i2 = i2 + 1_int32
            z = z + 1.0_dp / dble(n_tot_wlk)
         else
            w(i1) = dble(i2)
            i2 = 0_int32
            i1 = i1 + 1_int32
         endif
      enddo 
      call Init_ReplicaTable( n_tot_wlk )
      restart_walker = .false.
      if ( n_repli.gt.0_int32 .and. n_empty.gt.0_int32 ) call redistribute_IntraNode( )
#if defined _MPI || defined _MPIh || defined _MPI08
      if ( n_mpi_tasks.gt.1_int32 ) call redistribute_InterNode( )
#endif
      w = G_w_tot / dble(n_tot_wlk)
   end subroutine StochasticReconfigurationBranching
   subroutine Init_ReplicaTable ( n_max_weight )
      integer(int32), intent(in) :: n_max_weight
      integer(int32)  :: i1, i2
      if (n_max_weight .gt. n_max_repli) then
         n_max_repli = n_max_weight
         deallocate( repli_lst, empty_lst )
         allocate( repli_lst(1:2,1:n_max_repli), empty_lst(1:2,1:n_max_repli) )
      endif
      n_empty = 0_int32 ; n_repli = 0_int32
      repli_lst = 0_int32 ; empty_lst = 0_int32
      do i1 = 1_int32, n_wlk
         if ( w(i1).lt.1.0_dp) then
            n_empty = n_empty + 1_int32
            empty_lst(1,n_empty) = i1
            empty_lst(2,n_empty) = mpi_rank
         elseif ( w(i1).ge.2.0_dp ) then
            do i2 = 1_int32, int(w(i1)) - 1_int32
               n_repli = n_repli + 1_int32
               repli_lst(1,n_repli) = i1
               repli_lst(2,n_repli) = mpi_rank
            enddo
         endif
      enddo
   end subroutine Init_ReplicaTable
   subroutine redistribute_IntraNode
      integer(int32)  :: i1
      if ( n_repli.ge.n_empty ) then
         do i1 = 1_int32, n_empty
            call copy_walker( 0, empty_lst(1,i1), mpi_rank, repli_lst(1,i1), mpi_rank )
            empty_lst(:,i1) = 0_int32 ; repli_lst(:,i1) = 0_int32
         enddo
         repli_lst(1:2,1:n_repli-n_empty) = repli_lst(1:2,n_empty+1:n_repli)
         n_repli = n_repli-n_empty
         n_empty = 0_int32
      else
         do i1 = 1_int32, n_repli
            call copy_walker( 0, empty_lst(1,i1), mpi_rank, repli_lst(1,i1), mpi_rank )
            empty_lst(:,i1) = 0_int32 ; repli_lst(:,i1) = 0_int32
         enddo
         empty_lst(1:2,1:n_empty-n_repli) = empty_lst(1:2,n_repli+1:n_empty)
         n_empty = n_empty-n_repli
         n_repli = 0_int32
      endif
   end subroutine redistribute_IntraNode
#if defined _MPI || defined _MPIh || defined _MPI08
   subroutine redistribute_InterNode
      integer(int32)  :: i1
      call mpi_allreduce(n_empty, n_tot_empty, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr )
      call mpi_allreduce(n_repli, n_tot_repli, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr )
      if ( n_tot_repli.gt.0_int32 ) then
         er_indx = 0_int32 ; er_lst = 0_int32
         call mpi_allgather( 2*n_empty, 1, MPI_INTEGER, er_lst(:,1), 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr )
         call mpi_allgather( 2*n_repli, 1, MPI_INTEGER, er_lst(:,2), 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr )
         do i1 = 2_int32, n_mpi_tasks
            er_indx(i1,1) = er_indx(i1-1,1) + er_lst(i1-1,1)
            er_indx(i1,2) = er_indx(i1-1,2) + er_lst(i1-1,2)
         enddo
         if (mpi_rank.eq.0_int32) then
            call mpi_gatherv( MPI_IN_PLACE, 2*n_empty, MPI_INTEGER, empty_lst, &
            & er_lst(:,1), er_indx(:,1), MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr )
            call mpi_gatherv( MPI_IN_PLACE, 2*n_repli, MPI_INTEGER, repli_lst, &
            & er_lst(:,2), er_indx(:,2), MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr )
         else
            call mpi_gatherv( empty_lst(1:2,1:n_empty), 2*n_empty, MPI_INTEGER, empty_lst, &
            & er_lst(:,1), er_indx(:,1), MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr )
            call mpi_gatherv( repli_lst(1:2,1:n_repli), 2*n_repli, MPI_INTEGER, repli_lst, &
            & er_lst(:,2), er_indx(:,2), MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr )
         endif
         call mpi_bcast (repli_lst(1:2,1:n_tot_repli), 2*n_tot_repli, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr )
         call mpi_bcast (empty_lst(1:2,1:n_tot_empty), 2*n_tot_empty, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr )
         do i1 = 1_int32, n_tot_repli
            call copy_walker( i1, empty_lst(1,i1), empty_lst(2,i1), repli_lst(1,i1), repli_lst(2,i1) )
         enddo 
      endif
   end subroutine redistribute_InterNode
#endif
   subroutine copy_walker( icom, wlk1, node1, wlk2, node2 )
      integer(int32), intent(in) :: icom
      integer(int32), intent(in) :: wlk1
      integer(int32), intent(in) :: node1
      integer(int32), intent(in) :: wlk2
      integer(int32), intent(in) :: node2
#if defined _MPI || defined _MPIh || defined _MPI08
      if ( (node1.eq.mpi_rank) .and. (node2.eq.mpi_rank) ) then
#endif
         if ( wlk1.ne.wlk2 ) then
            if (molsys_prs) frm_cnf(wlk1)%r_fe = frm_cnf(wlk2)%r_fe
            if (qdosys_prs) drd_cnf(wlk1)%r_drd = drd_cnf(wlk2)%r_drd
            restart_walker(wlk1) = .true.
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
      else
         if ( mpi_rank.eq.node2 ) then
            if (molsys_prs) then
               call mpi_send( frm_cnf(wlk2)%r_fe, 3*n_fe, MPI_DOUBLE_PRECISION, node1, &
               & icom, MPI_COMM_WORLD, mpierr )
            endif
            if (qdosys_prs) then
               call mpi_send( drd_cnf(wlk2)%r_drd, 3*n_qdo, MPI_DOUBLE_PRECISION, node1, &
               & icom+n_wlk_max, MPI_COMM_WORLD, mpierr )
            endif
         else if ( mpi_rank.eq.node1 ) then
            if (molsys_prs) then
               call mpi_recv( frm_cnf(wlk1)%r_fe, 3*n_fe, MPI_DOUBLE_PRECISION, node2, &
               & icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr )
            endif
            if (qdosys_prs) then
               call mpi_recv( drd_cnf(wlk1)%r_drd, 3*n_qdo, MPI_DOUBLE_PRECISION, node2, &
               & icom+n_wlk_max, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr )
            endif
            restart_walker(wlk1) = .true.
         endif
      endif
#endif
   end subroutine copy_walker
   subroutine Fnlz_Branching
      deallocate( empty_lst, repli_lst )
#if defined _MPI || defined _MPIh || defined _MPI08
      deallocate( w_shift, er_lst, er_indx )
#endif
   end subroutine Fnlz_Branching
end module branching_m
