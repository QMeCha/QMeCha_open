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
module diffusion_monte_carlo_v
   use fortran_kinds_v,         only: dp, int32, stdout
   use openmp_mpi_m
   use write_lines_m
   implicit none
   integer(int32), save, public :: brnch_type
   integer(int32), save, public :: encrr_type
   integer(int32), save, public :: n_trm_dmc
   real(dp),       save, public :: e_r, e_r_avg, we_r
   real(dp),       save, public :: g_r
   real(dp),       save, public, allocatable, dimension(:) :: e_r_vec
   real(dp),       save, public :: g_e
   real(dp),       save, public, allocatable, dimension(:)   :: e_l_old
   real(dp),       save, public, allocatable, dimension(:,:) :: e_l_pp_old
   real(dp),       save, public, allocatable, dimension(:)   :: w
   integer(int32), save, public, allocatable, dimension(:)   :: wlk_life
   integer(int32), save, public                              :: max_wlk_life
   integer(int32), save, public, allocatable, dimension(:,:) :: fe_i_at_old
   logical,        save, public, allocatable, dimension(:)   :: restart_walker
   real(dp),       save, public                              :: dt_eff_avg, wdt_eff
   real(dp),       save, public, allocatable, dimension(:)   :: dt_eff
   real(dp),       save, public, allocatable, dimension(:)   :: pdr2, dr2
   real(dp),       save, public, allocatable, dimension(:)   :: V2
   real(dp),       save, public, allocatable, dimension(:,:) :: V2_pp
   real(dp),       save, public, pointer                   :: e_best, e_best_var
   real(dp),       save, public                            :: e_cut
   real(dp),       save, public                            :: w_blck, w2_blck
   real(dp),       save, public, target                    :: e_blck, e_var_blck
   real(dp),       save, public, allocatable, dimension(:,:), target :: spe_blck
   real(dp),       save, public, allocatable, dimension(:),   target :: spe_var_blck
   real(dp),       save, public, pointer, dimension(:)     :: spe_best, spe_best_var
   real(dp),       save, public                            :: w_chk, w_max, log_w_max
   integer(int32), save, public                            :: t_bra, trotter_order
   public :: allocate_dmcmthd_var, deallocate_dmcmthd_var
contains
   subroutine allocate_dmcmthd_var( n_particles, n_wlk, n_wlk_max, n_pp_obs, bin_l )
      integer(int32), intent(in) :: n_particles
      integer(int32), intent(in) :: n_wlk
      integer(int32), intent(in) :: n_wlk_max
      integer(int32), intent(in) :: n_pp_obs
      integer(int32), intent(in) :: bin_l
      allocate(w(1:n_wlk_max)) ; w(1:n_wlk) = 1.0_dp
      if (n_wlk.ne.n_wlk_max) w(n_wlk+1:n_wlk_max) = 0.0_dp
      allocate(e_l_old(1:n_wlk_max))      ; e_l_old = 0.0_dp
      allocate(e_l_pp_old(1:n_particles,1:n_wlk_max)) ; e_l_pp_old = 0.0_dp
      allocate( e_r_vec(1:bin_l) ) ; e_r_vec = 0.0_dp
      allocate( wlk_life(1:n_wlk_max) ) ; wlk_life = 0_int32
      allocate(restart_walker(1:n_wlk_max)) ; restart_walker=.true.
      allocate( dr2(1:n_wlk_max) ) ; dr2 = 0.0_dp
      allocate( pdr2(1:n_wlk_max) ) ; pdr2 = 0.0_dp
      select case(encrr_type)
       case(3,4,5)
         allocate(V2(1:n_wlk_max)) ; V2(1:n_wlk_max) = 0.0_dp
       case(13,14,15,24,25)
         allocate(V2_pp(1:n_particles,1:n_wlk_max)) ; V2_pp = 0.0_dp
       case default
      end select
      allocate(dt_eff(1:n_wlk_max)) ; dt_eff = 0.0_dp
      if (encrr_type.ge.20_int32) then
         allocate(spe_blck(1:n_pp_obs,1:4)) ; spe_blck = 0.0_dp
         allocate(spe_var_blck(1:n_pp_obs) ) ; spe_var_blck = 0.0_dp
         allocate(fe_i_at_old(1:n_particles,1:n_wlk_max) ) ; fe_i_at_old = 0_int32
      else if (encrr_type.ge.10_int32) then
         allocate(spe_blck(1:n_pp_obs,1:2)) ; spe_blck = 0.0_dp
         allocate(spe_var_blck(1:n_pp_obs) ) ; spe_var_blck = 0.0_dp
      endif
   end subroutine allocate_dmcmthd_var
   subroutine save_dmcmthd_var( n_wlk_max, bin_l )
      integer(int32), intent(in) :: n_wlk_max
      integer(int32), intent(in) :: bin_l
      real(dp),       allocatable, dimension(:) :: w_buf
      integer(int32), allocatable, dimension(:) :: i_buf
      character(100) :: svf_fle
      integer(int32) :: i_mpi_task
      svf_fle = 'qmc.save/dmc_variables.sav'
      if (mpi_rank.eq.0_int32) then
         open(unit=1,file=svf_fle,action='write',form='unformatted',status='unknown',access='sequential')
         !open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
         write(1) encrr_type, brnch_type, n_trm_dmc, trotter_order, g_r, g_e
      endif 
      if (mpi_rank.eq.0_int32) then
         write(1) w(1:n_wlk_max)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      if ( n_mpi_tasks.gt.1_int32 ) then
         if ( mpi_rank.eq.0_int32 ) allocate( w_buf(1:n_wlk_max) )
         do i_mpi_task = 1_int32, n_mpi_tasks-1
            if (mpi_rank.eq.i_mpi_task) then
               call mpi_send (w, n_wlk_max, MPI_DOUBLE_PRECISION, 0, i_mpi_task, MPI_COMM_WORLD, mpierr)
            else if (mpi_rank.eq.0) then
               call mpi_recv (w_buf, n_wlk_max, MPI_DOUBLE_PRECISION, i_mpi_task, i_mpi_task, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
            endif
            call mpi_barrier(MPI_COMM_WORLD, mpierr)
            if ( mpi_rank.eq.0_int32 ) write(1) w_buf(1:n_wlk_max)
         enddo
         if (allocated(w_buf)) deallocate( w_buf )
      endif
#endif
      if (mpi_rank.eq.0_int32) then
         write(1) dt_eff_avg, e_r, e_r_avg, e_r_vec(1:bin_l)
      endif
      if (mpi_rank.eq.0_int32) then
         write(1) wlk_life(1:n_wlk_max)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      if ( n_mpi_tasks.gt.1_int32 ) then
         if ( mpi_rank.eq.0_int32 ) allocate( i_buf(1:n_wlk_max) )
         do i_mpi_task = 1_int32, n_mpi_tasks-1
            if (mpi_rank.eq.i_mpi_task) then
               call mpi_send (wlk_life, n_wlk_max, MPI_INTEGER, 0, i_mpi_task, MPI_COMM_WORLD, mpierr)
            else if (mpi_rank.eq.0) then
               call mpi_recv (i_buf, n_wlk_max, MPI_INTEGER, i_mpi_task, i_mpi_task, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
            endif
            call mpi_barrier(MPI_COMM_WORLD, mpierr)
            if ( mpi_rank.eq.0_int32 ) write(1) i_buf(1:n_wlk_max)
         enddo
         if (allocated(i_buf)) deallocate( i_buf )
      endif
#endif
      if (mpi_rank.eq.0_int32) close(1)
   end subroutine save_dmcmthd_var
   subroutine load_dmcmthd_var( n_wlk_max, bin_l )
      integer(int32), intent(in) :: n_wlk_max
      integer(int32), intent(in) :: bin_l
      real(dp),       allocatable, dimension(:) :: w_buf
      integer(int32), allocatable, dimension(:) :: i_buf
      character(100) :: svf_fle
      logical        :: file_present
      integer(int32) :: i_mpi_task
      file_present = .false.
      if (mpi_rank.eq.0_int32) then
         svf_fle = 'qmc.save/dmc_variables.sav'
         inquire(file=svf_fle, exist=file_present)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(file_present, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif
      if (file_present) then
         call write_line_check(stdout,0,mpi_rank,2,"Loading Diffusion Monte Carlo variables")
         if (mpi_rank.eq.0_int32) then
            open(unit=1,file=svf_fle,action='read',form='unformatted',status='old',access='sequential')
            read(1) encrr_type, brnch_type, n_trm_dmc, trotter_order, g_r, g_e
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_bcast(encrr_type,    1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(brnch_type,    1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(n_trm_dmc,     1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(trotter_order, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(g_r,           1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(g_e,           1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
         if (mpi_rank.eq.0_int32) then
            read(1) w(1:n_wlk_max)
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
      if ( n_mpi_tasks.gt.1_int32 ) then
         if ( mpi_rank.eq.0_int32 ) allocate( w_buf(1:n_wlk_max) )
         do i_mpi_task = 1_int32, n_mpi_tasks-1
            if (mpi_rank.eq.0) then
               read(1) w_buf(1:n_wlk_max)
               call mpi_send (w_buf, n_wlk_max, MPI_DOUBLE_PRECISION, i_mpi_task, i_mpi_task, MPI_COMM_WORLD, mpierr)
            else if (mpi_rank.eq.i_mpi_task) then
               call mpi_recv (w, n_wlk_max, MPI_DOUBLE_PRECISION, 0, i_mpi_task, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
            endif
            call mpi_barrier(MPI_COMM_WORLD, mpierr)
         enddo
         if (allocated(w_buf)) deallocate( w_buf)
      endif
#endif
      if ( mpi_rank.eq.0_int32 ) then
         read(1) dt_eff_avg, e_r, e_r_avg, e_r_vec(1:bin_l)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(e_r_vec, bin_l, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(e_r,         1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(e_r_avg,     1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(dt_eff_avg,  1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
      if (mpi_rank.eq.0_int32) then
         read(1) wlk_life(1:n_wlk_max)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      if ( n_mpi_tasks.gt.1_int32 ) then
         if ( mpi_rank.eq.0_int32 ) allocate( i_buf(1:n_wlk_max) )
         do i_mpi_task = 1_int32, n_mpi_tasks-1
            if (mpi_rank.eq.0) then
               read(1) i_buf(1:n_wlk_max)
               call mpi_send (i_buf, n_wlk_max, MPI_INTEGER, i_mpi_task, i_mpi_task, MPI_COMM_WORLD, mpierr)
            else if (mpi_rank.eq.i_mpi_task) then
               call mpi_recv (wlk_life, n_wlk_max, MPI_INTEGER, 0, i_mpi_task, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
            endif
            call mpi_barrier(MPI_COMM_WORLD, mpierr)
         enddo
         if (allocated(i_buf)) deallocate( i_buf )
      endif
#endif
         if (mpi_rank.eq.0_int32) close(1)
         call write_done(stdout,0,mpi_rank)
      endif
   end subroutine load_dmcmthd_var
   subroutine deallocate_dmcmthd_var()
      if (allocated(w)) deallocate(w)
      if (allocated(e_l_old)) deallocate(e_l_old)
      if (allocated(e_l_pp_old)) deallocate(e_l_pp_old)
      if (allocated(e_r_vec)) deallocate( e_r_vec )
      if (allocated(wlk_life)) deallocate( wlk_life )
      if (allocated(restart_walker)) deallocate(restart_walker)
      if (allocated(dr2)) deallocate( dr2 )
      if (allocated(pdr2)) deallocate( pdr2 )
      if (allocated(V2)) deallocate( V2 )
      if (allocated(V2_pp) ) deallocate( V2_pp )
      if (allocated(dt_eff) ) deallocate( dt_eff )
      if (allocated(spe_blck) ) deallocate( spe_blck )
      if (allocated(spe_var_blck) ) deallocate( spe_var_blck )
      if (allocated(fe_i_at_old) ) deallocate( fe_i_at_old )
   end subroutine deallocate_dmcmthd_var
end module diffusion_monte_carlo_v
