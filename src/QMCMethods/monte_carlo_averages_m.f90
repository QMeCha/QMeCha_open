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
module monte_carlo_averages_m
   use fortran_kinds_v,         only: int32, dp
   use openmp_mpi_m
   use write_lines_m
   use merge_sort_m,            only: mrg_srt
   use physical_constants_v,    only: au_to_debye, au_to_ang
   use monte_carlo_averages_v
   use timings_m
   use molecular_system_v,      only: n_fe, n_el, n_po, n_at
   use qdo_system_m,            only: n_qdo
   use quantum_monte_carlo_v,   only: G_w_trgt, n_tot_wlk, n_wlk, sysname, &
   & dt_d, dt_e, dt_p, dt, n_bin, bin_l, n_bin_old, bin_l_old, qmc_mthd, n_wlk_max, &
   & save_status, save_step, n_tot_smpl, n_tot_smpl_old
   use diffusion_monte_carlo_v, only: w, encrr_type, dt_eff, e_r_vec, &
   & dt_eff_avg, e_r, e_r_avg
   use local_energy_m,          only: v_ext, e_l, v_l, k_l, e_l_pp, v_fn, v_ee, v_dd, v_dq
   use electronic_properties_m, only: dipole, quadrupole, dip_mom, quad_mom, expec_dist, expec_pote
   use fermionic_config_c,      only: frm_cnf
   use drudonic_config_c,       only: drd_cnf
   implicit none
   public  :: init_qmcsmp, empty_qmcsmp, accu_qmcavg, write_qmcavg, &
   & bld_obsvec, save_qmcsmp, load_qmcsmp
contains
   subroutine init_qmcsmp( )
      integer(int32) :: i1, i2, i3, i4, i5, ierr
      integer(int32) :: order_of_primes(4), n_max_block_size
      integer(int32), allocatable, dimension(:) :: tmp_block_sizes
      call decomposition_primes (n_tot_smpl, order_of_primes(1:4) )
      n_max_block_size = 1_int32
      do i1 = 1_int32, 4_int32
         if (order_of_primes(i1).gt.0_int32) n_max_block_size = n_max_block_size * (order_of_primes(i1)+1)
      enddo
      allocate(tmp_block_sizes(1:n_max_block_size)) ; tmp_block_sizes = 0_int32
      tmp_block_sizes(1) = 1_int32
      i5 = 1_int32
      do i1 = 0_int32,order_of_primes(1)
         do i2 = 0_int32,order_of_primes(2)
            do i3 = 0_int32,order_of_primes(3)
               do i4 = 0_int32,order_of_primes(4)
                  tmp_block_sizes(i5) = 1_int32
                  tmp_block_sizes(i5) = tmp_block_sizes(i5) * 2**(i1)
                  tmp_block_sizes(i5) = tmp_block_sizes(i5) * 3**(i2)
                  tmp_block_sizes(i5) = tmp_block_sizes(i5) * 5**(i3)
                  tmp_block_sizes(i5) = tmp_block_sizes(i5) * 7**(i4)
                  i5 = i5 + 1_int32
               enddo
            enddo
         enddo
      enddo
      call mrg_srt(n_max_block_size, tmp_block_sizes(1:n_max_block_size) )
      i2 = 1_int32
      do i1 = 1_int32, n_max_block_size-1
         if ( tmp_block_sizes(i2).gt.10000) exit
         i2 = i2 + 1_int32
      enddo
      n_max_block_size = i2 - 1_int32
      min_blck_size = tmp_block_sizes(1)
      max_blck_size = tmp_block_sizes(n_max_block_size)
      n_blck_sizes = n_max_block_size
      allocate( blck_sizes(1:n_blck_sizes) ) ; blck_sizes(1:n_blck_sizes) = tmp_block_sizes(1:n_blck_sizes)
      deallocate( tmp_block_sizes )
      n_obs = 6_int32
      if ( expec_pote ) then
         if ( n_fe.gt.0_int32 )  n_obs = n_obs + n_fe + 1
         if ( n_qdo.gt.0_int32 ) n_obs = n_obs + n_qdo + 1
      endif
      if (dip_mom) n_obs = n_obs + 3_int32
      if (quad_mom) n_obs = n_obs + 6_int32
      if ( expec_dist ) then
         n_obs = n_obs + 3_int32 * ( n_qdo + n_fe )
         if ( n_el.ge.2_int32 )  n_obs = n_obs + 1_int32
         if ( n_po.ge.2_int32 )  n_obs = n_obs + 1_int32
         if ( n_qdo.ge.2_int32 ) n_obs = n_obs + 1_int32
         if (n_el.gt.0_int32.and.n_po.gt.0_int32)  n_obs = n_obs + 1_int32
      endif
      if (qmc_mthd.eq.'dmc') then
         n_obs = n_obs + 1_int32 
      endif
      allocate( obs_avg(1:n_obs) ) ; obs_avg = 0.0_dp
      allocate( obs_var(1:2) )  ; obs_var = 0.0_dp
      n_pp_obs = 0_int32
      if (qmc_mthd.eq.'dmc') then
         if (encrr_type.ge.20_int32) then
            if ( n_el.gt.0_int32 )  n_pp_obs = n_pp_obs + n_at
            if ( n_po.gt.0_int32 )  n_pp_obs = n_pp_obs + n_at
            if ( n_qdo.gt.0_int32 ) n_pp_obs = n_pp_obs + n_qdo
            allocate( obs_pp_avg(1:n_pp_obs,1:4) ) ; obs_pp_avg = 0.0_dp
            allocate( obs_pp_b(1:n_pp_obs,1:4,1:bin_l) ) ; obs_pp_b = 0.0_dp
         else if (encrr_type.ge.10_int32) then
            n_pp_obs = n_qdo+n_fe
            allocate( obs_pp_avg(1:n_pp_obs,1:2) )  ; obs_pp_avg = 0.0_dp
            allocate( obs_pp_b(1:n_pp_obs,1:2,1:bin_l) ) ; obs_pp_b = 0.0_dp
         endif
         if (encrr_type.ge.10_int32) then
            allocate( obs_pp_var(1:n_pp_obs) ) ; obs_pp_var = 0.0_dp
         endif
      endif
      allocate( obs_b(1:n_obs,1:bin_l) ) ; obs_b = 0.0_dp
      allocate( blck_l(1:n_blck_sizes) ) ; blck_l = 0_int32
      allocate( blck_avg(1:n_obs-2,1:n_blck_sizes) )  ; blck_avg = 0.0_dp
      allocate( blck_var(1:n_obs-2,1:n_blck_sizes) )  ; blck_var = 0.0_dp
      allocate( blck_err(1:n_obs-2,1:n_blck_sizes) )  ; blck_err = 0.0_dp
      i1 = 1_int32
      w_avg    => obs_avg(i1) ; i1 = i1 + 1_int32
      e_avg    => obs_avg(i1) ; i1 = i1 + 1_int32
      k_avg    => obs_avg(i1) ; i1 = i1 + 1_int32
      v_avg    => obs_avg(i1) ; i1 = i1 + 1_int32
      if ( expec_pote ) then
         if ( n_fe.gt.0_int32 )  then
            v_ee_avg => obs_avg(i1) ; i1 = i1 + 1_int32
            v_fn_avg => obs_avg(i1:i1+n_fe-1) ; i1 = i1 + n_fe
         endif
         if ( n_qdo.gt.0_int32 ) then
            v_dd_avg => obs_avg(i1) ; i1 = i1 + 1_int32
            v_dq_avg => obs_avg(i1:i1+n_qdo-1) ; i1 = i1 + n_qdo
         endif
      endif
      if (dip_mom) then
         dip_avg => obs_avg(i1:i1+2) ; i1 = i1 + 3_int32
      endif
      if (quad_mom) then
         qud_avg => obs_avg(i1:i1+5) ; i1 = i1 + 6_int32
      endif
      if (expec_dist) then
         if (n_fe.gt.0_int32) then
            r_fe_avg => obs_avg(i1:i1+3*n_fe-1) ; i1 = i1 + 3*n_fe
         endif
         if (n_qdo.gt.0_int32) then
            r_drd_avg => obs_avg(i1:i1+3*n_qdo-1) ; i1 = i1 + 3*n_qdo
         endif
         if (n_el.ge.2_int32) then
            d_ee_avg  => obs_avg(i1) ; i1 = i1 + 1_int32
         endif
         if (n_po.ge.2_int32) then
            d_pp_avg  => obs_avg(i1) ; i1 = i1 + 1_int32
         endif
         if (n_el.gt.0_int32.and.n_po.gt.0_int32) then
            d_ep_avg  => obs_avg(i1) ; i1 = i1 + 1_int32
         endif
         if (n_qdo.gt.1_int32) then
            d_dd_avg  => obs_avg(i1) ; i1 = i1 + 1_int32
         endif
      endif
      w2_avg   => obs_avg(n_obs-1)
      e_var    => obs_var(2)
      if (qmc_mthd.eq.'dmc') then
         if (encrr_type.ge.20_int32) then
            e_pp_avg => obs_pp_avg(1:n_pp_obs,2)
         else if (encrr_type.ge.10_int32) then
            e_pp_avg => obs_pp_avg(1:n_pp_obs,1)
         endif
         if (encrr_type.ge.10_int32) then
            e_pp_var => obs_pp_var(1:n_pp_obs)
         endif
      endif
   end subroutine init_qmcsmp
   subroutine save_qmcsmp( )
      character(50) :: svf_fle
      svf_fle = 'qmc.save/eval_stats.sav'
      if (mpi_rank.eq.0_int32) then
         open(unit=1,file=svf_fle,action='write',form='unformatted',status='unknown',access='sequential')
         if ( save_status .ge. 2_int32 ) then
            write(1) dip_mom, quad_mom, expec_dist, expec_pote
            write(1) n_blck_sizes, n_obs, n_pp_obs
            write(1) blck_sizes(1:n_blck_sizes)
            write(1) blck_var(1:n_obs-2,1:n_blck_sizes)
            write(1) obs_avg(1:n_obs)
            write(1) obs_var(1:2)
            if (qmc_mthd.eq.'dmc') then
               if (encrr_type.ge.20_int32) then
                  write(1) obs_pp_avg(1:n_pp_obs,1:4)
               else if (encrr_type.ge.10_int32) then
                  write(1) obs_pp_avg(1:n_pp_obs,1:2)
               endif
               if (encrr_type.ge.10_int32) then
                  write(1) obs_pp_var(1:n_pp_obs)
               endif
               write(1) obs_b(1:n_obs,1:bin_l)
               if (encrr_type.ge.20_int32) then
                  write(1) obs_pp_b(1:n_pp_obs,1:4,1:bin_l)
               else if (encrr_type.ge.10_int32) then
                  write(1) obs_pp_b(1:n_pp_obs,1:2,1:bin_l)
               endif
            endif
         endif
         close(1)
      endif
   end subroutine save_qmcsmp
   subroutine load_qmcsmp(  )
      logical        :: file_present
      integer(int32) :: n_blck_sizes_old, n_obs_old, n_pp_obs_old
      integer(int32) :: iw, i_mpi_task, i_block, n_compatible_blocks, i_block_old
      character(3)   :: qmc_mthd_old
      integer(int32) :: encrr_type_old, n_tot_smpl_tot
      character(50)  :: svf_fle
      integer(int32) :: order_of_primes(4)
      integer(int32), allocatable, dimension(:)   :: blck_sizes_old
      logical,        allocatable, dimension(:)   :: block_compatibility
      real(dp),       allocatable, dimension(:,:) :: blck_var_old
      real(dp),       allocatable, dimension(:)   :: w_buf
      file_present = .false.
      if ( ( save_status.eq.2_int32 .and. qmc_mthd.eq.'vmc').or.&
      & ( save_status.eq.4_int32 .and. qmc_mthd.eq.'dmc') ) then
         if (mpi_rank.eq.0_int32) then
            svf_fle = 'qmc.save/eval_stats.sav'
            inquire(file=svf_fle, exist=file_present)
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_bcast(file_present, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif
         if (allocated(obs_b)) then
            deallocate(obs_b)
            allocate( obs_b(1:n_obs,1:bin_l) ) ; obs_b = 0.0_dp
         endif
         if (qmc_mthd.eq.'dmc') then
            deallocate(e_r_vec)
            allocate(e_r_vec(1:bin_l)) ; e_r_vec = 0.0_dp
            if (allocated(obs_pp_b)) then
               deallocate(obs_pp_b)
               if (encrr_type.ge.20_int32) then
                  allocate( obs_pp_b(1:n_pp_obs,1:4,1:bin_l) )
               else if (encrr_type.ge.10_int32) then
                  allocate( obs_pp_b(1:n_pp_obs,1:2,1:bin_l) )
               endif
            endif
         endif
         if (file_present) then
            call write_line_check(stdout,0,mpi_rank,2,"Loading Monte Carlo sampling")
            if (mpi_rank.eq.0_int32) then
               open(unit=1,file=svf_fle,action='read',form='unformatted',status='old',access='sequential')
               read(1) dip_mom, quad_mom, expec_dist, expec_pote
               read(1) n_blck_sizes_old, n_obs_old, n_pp_obs_old
            endif
#if defined _MPI || defined _MPIh || defined _MPI08
            call mpi_bcast(    dip_mom, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(   quad_mom, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast( expec_dist, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast( expec_pote, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast( n_blck_sizes_old, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(        n_obs_old, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(     n_pp_obs_old, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
            allocate( blck_sizes_old(1:n_blck_sizes_old) )
            allocate( block_compatibility(1:n_blck_sizes_old) ) ; block_compatibility = .false.
            if (mpi_rank.eq.0_int32) then
               read(1) blck_sizes_old(1:n_blck_sizes_old)
            endif
#if defined _MPI || defined _MPIh || defined _MPI08
            call mpi_bcast(blck_sizes_old, n_blck_sizes_old, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
            n_compatible_blocks = 0_int32
            n_tot_smpl_tot = n_tot_smpl + n_tot_smpl_old
            do i_block = 1_int32, n_blck_sizes_old
               if ( mod( n_tot_smpl_tot, blck_sizes_old(i_block)).eq.0_int32 ) then
                  block_compatibility(i_block) = .true.
                  n_compatible_blocks = n_compatible_blocks + 1_int32
               endif
            enddo
            if ( n_compatible_blocks.ne.n_blck_sizes_old ) then
               deallocate( blck_l, blck_avg, blck_var, blck_err, blck_sizes )
               n_blck_sizes = n_compatible_blocks
               allocate( blck_l(1:n_blck_sizes) ) ; blck_l = 0_int32
               allocate( blck_avg(1:n_obs-2,1:n_blck_sizes) )  ; blck_avg = 0.0_dp
               allocate( blck_var(1:n_obs-2,1:n_blck_sizes) )  ; blck_var = 0.0_dp
               allocate( blck_err(1:n_obs-2,1:n_blck_sizes) )  ; blck_err = 0.0_dp
               allocate( blck_sizes(1:n_blck_sizes) )
               allocate( blck_var_old(1:n_obs-2,1:n_blck_sizes_old))
               if (mpi_rank.eq.0_int32) then
                  read(1) blck_var_old(1:n_obs-2,1:n_blck_sizes_old)
               endif
#if defined _MPI || defined _MPIh || defined _MPI08
               call mpi_bcast( blck_var_old, (n_obs-2)*n_blck_sizes_old, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
               i_block = 1_int32
               do i_block_old = 1_int32, n_blck_sizes_old
                  if ( block_compatibility(i_block_old) ) then
                     blck_sizes(i_block) = blck_sizes_old(i_block_old)
                     blck_var(1:n_obs-2,i_block) = blck_var_old(1:n_obs-2,i_block_old)
                     i_block = i_block + 1_int32
                  endif
               enddo
               deallocate (blck_var_old)
            else
               if (n_blck_sizes_old.ne.n_blck_sizes) then 
                  n_blck_sizes = n_blck_sizes_old
                  deallocate( blck_l, blck_avg, blck_var, blck_err, blck_sizes )
                  allocate( blck_l(1:n_blck_sizes) ) ; blck_l = 0_int32
                  allocate( blck_avg(1:n_obs-2,1:n_blck_sizes) )  ; blck_avg = 0.0_dp
                  allocate( blck_var(1:n_obs-2,1:n_blck_sizes) )  ; blck_var = 0.0_dp
                  allocate( blck_err(1:n_obs-2,1:n_blck_sizes) )  ; blck_err = 0.0_dp
                  allocate( blck_sizes(1:n_blck_sizes) )
               endif
               blck_sizes(1:n_blck_sizes) = blck_sizes_old(1:n_blck_sizes)
               if (mpi_rank.eq.0_int32) then
                  read(1) blck_var(1:n_obs-2,1:n_blck_sizes)
               endif
#if defined _MPI || defined _MPIh || defined _MPI08
               call mpi_bcast(blck_var, (n_obs-2)*n_blck_sizes, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
            endif
            if (mpi_rank.eq.0_int32) then
               read(1) obs_avg(1:n_obs)
               read(1) obs_var(1:2)
               if (qmc_mthd.eq.'dmc') then
                  if (encrr_type.ge.20_int32) then
                     read(1) obs_pp_avg(1:n_pp_obs,1:4)
                  else if (encrr_type.ge.10_int32) then
                     read(1) obs_pp_avg(1:n_pp_obs,1:2)
                  endif
                  if (encrr_type.ge.10_int32) then
                     read(1) obs_pp_var(1:n_pp_obs)
                  endif
                  read(1) obs_b(1:n_obs,1:bin_l)
                  if (encrr_type.ge.20_int32) then
                     read(1) obs_pp_b(1:n_pp_obs,1:4,1:bin_l)
                  else if (encrr_type.ge.10_int32) then
                     read(1) obs_pp_b(1:n_pp_obs,1:2,1:bin_l)
                  endif
               endif
            endif
#if defined _MPI || defined _MPIh || defined _MPI08
            call mpi_bcast(obs_avg,     n_obs, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(obs_var,         2, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(obs_b, n_obs*bin_l, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
            if (qmc_mthd.eq.'dmc') then
               if (encrr_type.ge.20_int32) then
                  call mpi_bcast(obs_pp_avg,       n_pp_obs*4, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
                  call mpi_bcast(  obs_pp_b, n_pp_obs*4*bin_l, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
               else if (encrr_type.ge.10_int32) then
                  call mpi_bcast(obs_pp_avg,       n_pp_obs*2, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
                  call mpi_bcast(  obs_pp_b, n_pp_obs*2*bin_l, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
               endif
               if (encrr_type.ge.10_int32) then
                  call mpi_bcast(obs_pp_var, n_pp_obs, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
               endif
            endif
#endif
            if (mpi_rank.eq.0_int32) close(1)
            call write_done(stdout,0,mpi_rank)
         endif
      else
         call write_empty_line(stdout,0,mpi_rank)
         call write_simple_line(stdout,0,mpi_rank,2,"l","No requirement to load previous sampling.")
      endif
   end subroutine load_qmcsmp
   subroutine empty_qmcsmp
      if(allocated(obs_pp_avg)) obs_pp_avg = 0.0_dp
      if(allocated(obs_pp_var)) obs_pp_var = 0.0_dp
      obs_avg = 0.0_dp ; obs_var = 0.0_dp
      blck_l = 0_int32
      blck_avg = 0.0_dp
      blck_var = 0.0_dp
      blck_err = 0.0_dp
   end subroutine empty_qmcsmp
   subroutine bld_obsvec( obs_vec, obs_pp_vec )
      real(dp),           dimension(n_obs),       intent(out) :: obs_vec
      real(dp), optional, dimension(n_pp_obs,*),  intent(out) :: obs_pp_vec
      integer(int32)                               :: io, iw, i1, i_min
      call time_compute(t_avg_step)
      obs_vec = 0.0_dp
      io = 1_int32
      obs_vec(io) = sum(w(1:n_wlk)) ; io = io + 1_int32
      obs_vec(io) = dot_product ( w(1:n_wlk) , e_l(1:n_wlk) ) ; io = io + 1_int32
      obs_vec(io) = dot_product ( w(1:n_wlk) , k_l(1:n_wlk) ) ; io = io + 1_int32
      obs_vec(io) = dot_product ( w(1:n_wlk) , v_l(1:n_wlk) ) ; io = io + 1_int32
      if ( expec_pote ) then
         if ( n_fe.gt.0_int32 )  then
            do iw = 1_int32, n_wlk
               obs_vec(io) = obs_vec(io) + w(iw) *sum(v_ee(:,iw))
               obs_vec(io+1:io+n_fe) = obs_vec(io+1:io+n_fe) + w(iw) *v_fn(:,iw)
            enddo
            io = io + n_fe + 1
         endif
         if ( n_qdo.gt.0_int32 ) then
            do iw = 1_int32, n_wlk
               obs_vec(io) = obs_vec(io) + w(iw) * sum(v_dd(:,iw))
               obs_vec(io+1:io+n_qdo) = obs_vec(io+1:io+n_qdo) + w(iw) * v_dq(:,iw)
            enddo
            io = io + n_qdo + 1
         endif
      endif
      if (dip_mom) then
         do iw = 1_int32, n_wlk
            obs_vec(io:io+2_int32) = obs_vec(io:io+2_int32) + dipole(1:3,iw) * w(iw)
         enddo
         io = io + 3_int32
      endif
      if (quad_mom) then
         do iw = 1_int32, n_wlk
            obs_vec(io:io+5_int32) = obs_vec(io:io+5_int32) + quadrupole(1:6,iw) * w(iw)
         enddo
         io = io + 6_int32
      endif
      if (expec_dist) then
         if (n_fe.gt.0_int32) then
            do i1 = 1_int32, n_fe ; do iw = 1_int32, n_wlk
                  obs_vec(io:io+2_int32) = obs_vec(io:io+2_int32) + frm_cnf(iw)%r_fe(1:3,i1) * w(iw)
               enddo
               io = io + 3_int32
            enddo
         endif
         if (n_qdo.gt.0_int32) then
            do i1 = 1_int32, n_qdo ; do iw = 1_int32, n_wlk
                  obs_vec(io:io+2_int32) = obs_vec(io:io+2_int32) + drd_cnf(iw)%r_drd(1:3,i1) * w(iw)
               enddo
               io = io + 3_int32
            enddo
         endif
         if (n_el.gt.1_int32) then
            do iw = 1_int32, n_wlk
               obs_vec(io) = obs_vec(io) + sum( frm_cnf(iw)%d_ee(:)%m)/ dble(n_el*(n_el-1)/2) * w(iw)
            enddo
            io = io + 1_int32
         endif
         if (n_po.gt.1_int32) then
            do iw = 1_int32, n_wlk
               obs_vec(io) = obs_vec(io) + sum( frm_cnf(iw)%d_pp(:)%m) / dble(n_po*(n_po-1)/2) * w(iw)
            enddo
            io = io + 1_int32
         endif
         if (n_el.gt.0_int32.and.n_po.gt.0_int32) then
            do iw = 1_int32, n_wlk
               obs_vec(io) = obs_vec(io) + sum( frm_cnf(iw)%d_ep(:,:)%m) * w(iw)/ dble(n_po*n_el)
            enddo
            io = io + 1_int32
         endif
         if (n_qdo.gt.1_int32) then
            do iw = 1_int32, n_wlk
               obs_vec(io) = obs_vec(io) + sum( drd_cnf(iw)%d_dd(1:n_qdo*(n_qdo-1)/2)%m)/ dble(n_qdo*(n_qdo-1)/2) * w(iw)
            enddo
            io = io + 1_int32
         endif
      endif
      if (qmc_mthd.eq.'dmc') then
         obs_vec(io) = dot_product ( dt_eff(1:n_wlk) , w(1:n_wlk) )
         io = io + 1_int32
      endif
      obs_vec(io) = dot_product ( w(1:n_wlk) , w(1:n_wlk) ) / G_w_trgt ; io = io + 1_int32
      obs_vec(io) = dot_product ( w(1:n_wlk) , e_l(1:n_wlk)**2 ) ; io = io + 1_int32
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_allreduce(MPI_IN_PLACE, obs_vec, n_obs, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
      obs_vec = obs_vec / G_w_trgt
      if (present(obs_pp_vec)) then
         io = 0_int32
         if (encrr_type.ge.20_int32) then
            obs_pp_vec(1:n_pp_obs,1:4) = 0.0_dp
            if (n_el.gt.0_int32) then
               do iw = 1_int32, n_wlk
                  do i1 = 1_int32, n_el
                     i_min = frm_cnf(iw)%fe_i_at(i1)
                     obs_pp_vec(io+i_min,1) = obs_pp_vec(io+i_min,1) + w(iw)
                     obs_pp_vec(io+i_min,2) = obs_pp_vec(io+i_min,2) + w(iw) * e_l_pp(i1,iw)
                     obs_pp_vec(io+i_min,3) = obs_pp_vec(io+i_min,3) + w(iw)**2  / G_w_trgt
                     obs_pp_vec(io+i_min,4) = obs_pp_vec(io+i_min,4) + w(iw) * e_l_pp(i1,iw) **2
                  enddo
               enddo
               io = io + n_at
            endif
            if (n_po.gt.0_int32) then
               do iw = 1_int32, n_wlk
                  do i1 = n_el+1_int32, n_fe
                     i_min = frm_cnf(iw)%fe_i_at(i1)
                     obs_pp_vec(io+i_min,1) = obs_pp_vec(io+i_min,1) + w(iw)
                     obs_pp_vec(io+i_min,2) = obs_pp_vec(io+i_min,2) + w(iw) * e_l_pp(i1,iw)
                     obs_pp_vec(io+i_min,3) = obs_pp_vec(io+i_min,3) + w(iw)**2  / G_w_trgt
                     obs_pp_vec(io+i_min,4) = obs_pp_vec(io+i_min,4) + w(iw) * e_l_pp(i1,iw) **2
                  enddo
               enddo
               io = io + n_at
            endif
#if defined _MPI || defined _MPIh || defined _MPI08
            call mpi_allreduce(MPI_IN_PLACE, obs_pp_vec(1:n_pp_obs,1:4), n_pp_obs*4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
            obs_pp_vec(1:n_pp_obs,1:4) = obs_pp_vec(1:n_pp_obs,1:4) / G_w_trgt
         else if (encrr_type.ge.10_int32) then
            obs_pp_vec(1:n_pp_obs,1:2) = 0.0_dp
            do iw = 1_int32, n_wlk
               obs_pp_vec(1:n_pp_obs,1) = obs_pp_vec(1:n_pp_obs,1) + w(iw) * e_l_pp(1:n_pp_obs,iw)
               obs_pp_vec(1:n_pp_obs,2) = obs_pp_vec(1:n_pp_obs,2) + w(iw) * e_l_pp(1:n_pp_obs,iw)**2
            enddo
#if defined _MPI || defined _MPIh || defined _MPI08
            call mpi_allreduce(MPI_IN_PLACE, obs_pp_vec(1:n_pp_obs,1:2), n_pp_obs*2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr )
#endif
            obs_pp_vec(1:n_pp_obs,1:2) = obs_pp_vec(1:n_pp_obs,1:2) / G_w_trgt
         endif
      endif 
      call time_compute_difference(t_avg_step,time_cpu_tot=t_avg_tot)
   end subroutine bld_obsvec
   subroutine accu_qmcavg( obs_vec, obs_pp_vec )
      real(dp),           dimension(n_obs),      intent(in) :: obs_vec
      real(dp), optional, dimension(n_pp_obs,*), intent(in) :: obs_pp_vec
      real(dp)       :: R
      integer(int32) :: i1
      obs_avg(1) = obs_avg(1) + obs_vec(1)
      obs_avg(2:n_obs-2) = obs_avg(2:n_obs-2) + obs_vec(1) / obs_avg(1) * ( obs_vec(2:n_obs-2) / obs_vec(1)-obs_avg(2:n_obs-2))
      obs_avg(n_obs-1) = obs_avg(n_obs-1) + obs_vec(n_obs-1)
      obs_avg(n_obs) = obs_avg(n_obs) + obs_vec(1) / obs_avg(1) * ( obs_vec(n_obs) / obs_vec(1) - obs_avg(n_obs))
      obs_var(1) = (obs_avg(n_obs-1) - obs_avg(1)**2)
      obs_var(2) = (obs_avg(n_obs) - obs_avg(2)**2)
      R = -obs_var(1) 
      if (R.gt.0.0_dp) then
         R = obs_avg(1)**2 / R
         obs_var(2) = obs_var(2) * R
      endif
      do i1 = 1_int32, n_blck_sizes
         blck_avg(1:n_obs-2,i1) = blck_avg(1:n_obs-2,i1) + obs_vec(1:n_obs-2)
         blck_l(i1) = blck_l(i1) + 1_int32
         if ( blck_l(i1).eq.blck_sizes(i1) ) then
            blck_var(1,i1)       = blck_var(1,i1)       + blck_avg(1,i1)**2
            blck_var(2:n_obs-2,i1) = blck_var(2:n_obs-2,i1) + blck_avg(2:n_obs-2,i1)**2 / blck_avg(1,i1)
            blck_avg(:,i1) = 0.0_dp
            blck_l(i1) = 0_int32
         endif
      enddo
      if ( present(obs_pp_vec)) then
         if (encrr_type.ge.20_int32) then
            obs_pp_avg(:,2)  = obs_pp_avg(:,2) + obs_pp_vec(:,1) / obs_pp_avg(:,1) * ( obs_pp_vec(:,2) / obs_pp_vec(:,1)-obs_pp_avg(:,2))
            obs_pp_avg(:,4)  = obs_pp_avg(:,4) + obs_pp_vec(:,1) / obs_pp_avg(:,1) * ( obs_pp_vec(:,4) / obs_pp_vec(:,1)-obs_pp_avg(:,4))
            obs_pp_var(:) = (obs_pp_avg(:,4) - obs_pp_avg(:,2)**2) * obs_pp_avg(:,1)**2/ (obs_pp_avg(:,1)**2 - obs_pp_avg(:,3))
         else if (encrr_type.ge.10_int32) then
            obs_pp_avg(:,1)  = obs_pp_avg(:,1) + obs_vec(1) / obs_avg(1) * ( obs_pp_vec(:,1) / obs_vec(1)-obs_pp_avg(:,1))
            obs_pp_avg(:,2)  = obs_pp_avg(:,2) + obs_vec(1) / obs_avg(1) * ( obs_pp_vec(:,2) / obs_vec(1)-obs_pp_avg(:,2))
            obs_pp_var(:) = obs_pp_avg(:,2) - obs_pp_avg(:,1)**2
            R = -obs_var(1) 
            if (R.gt.0.0_dp) then
               R = obs_avg(1)**2 / R
               obs_pp_var(:) = obs_pp_var(:) * R
            endif
         endif
      endif
   end subroutine accu_qmcavg
   subroutine write_qmcavg( )
      real(dp) :: q_iso, q_iso_blck_err
      real(dp) :: d_sqr, d_sqr_blck_err
      real(dp) :: C, diff, err1, err2
      character(100) :: fle_name
      integer(int32) :: correct_blck_length, n_blcks
      integer(int32) :: i1, i2
      character(10) :: variable_name
      if (mpi_rank.eq.0_int32) then
         if (trim(sysname).eq.'non' ) then
            if (qmc_mthd.eq.'vmc') then
               fle_name='vmc_eval.dat'
            else
               fle_name='dmc_eval.dat'
            endif
         else
            if (qmc_mthd.eq.'vmc') then
               fle_name='vmc_eval_'//trim(sysname)//'.dat'
            else
               fle_name='dmc_eval_'//trim(sysname)//'.dat'
            endif
         endif
         open(unit=22,file=fle_name,form='formatted',status='unknown', action='write')
      endif
      call write_separator_line(22,0,mpi_rank,2,"=")
      call write_simple_line(22,0,mpi_rank,2,"c","Evaluation of Monte Carlo averages")
      if (qmc_mthd.eq.'vmc') then
         call write_simple_line(22,0,mpi_rank,2,"c","Variational Monte Carlo")
         call write_empty_line(22,0,mpi_rank)
      else
         call write_simple_line(22,0,mpi_rank,2,"c","Diffusion Monte Carlo")
         call write_simple_line(22,0,mpi_rank,2,"c","(Mixed Averages)")
         call write_empty_line(22,0,mpi_rank)
      endif
      call write_empty_line(22,0,mpi_rank)
      call write_variable_line(22,0,mpi_rank,2,"Total number of samplings", n_tot_smpl+n_tot_smpl_old, var_name="n_tot_smpl")
      call write_variable_line(22,0,mpi_rank,2,"Maximum block size", max_blck_size, var_name="max_blck_size")
      call write_variable_line(22,0,mpi_rank,2,"Number block sizes", n_blck_sizes, var_name="n_blck_sizes")
      call write_empty_line(22,0,mpi_rank)
      call write_simple_line(22,0,mpi_rank,3,"l","- Computing total energy")
      if (dip_mom)    call write_simple_line(22,0,mpi_rank,3,"l","- Computing dipole moment")
      if (quad_mom)   call write_simple_line(22,0,mpi_rank,3,"l","- Computing quadrupole moment")
      if (expec_dist) call write_simple_line(22,0,mpi_rank,3,"l","- Computing Expected distances")
      if (expec_pote) call write_simple_line(22,0,mpi_rank,3,"l","- Computing potential energy components")
      call write_empty_line(22,0,mpi_rank)
      if(n_el.gt.0_int32) &
      & call write_variable_line(22,0,mpi_rank,2,"Final Time step for electrons", dt_e, var_name="dt_e")
      if(n_po.gt.0_int32) &
      & call write_variable_line(22,0,mpi_rank,2,"Final Time step for positrons", dt_p, var_name="dt_p")
      if(n_qdo.gt.0_int32) &
      & call write_variable_line(22,0,mpi_rank,2,"Final Time step for drudons", dt_d, var_name="dt_d")
      call write_empty_line(22,0,mpi_rank)
      do i1 = 1_int32, n_blck_sizes
         blck_var(2:n_obs-2,i1) = (blck_var(2:n_obs-2,i1) / obs_avg(1) - obs_avg(2:n_obs-2)**2  )  * obs_avg(1)**2 / ( obs_avg(1)**2-blck_var(1,i1) )
         blck_err(2:n_obs-2,i1) = sqrt( abs(blck_var(2:n_obs-2,i1)) * blck_var(1,i1)  )/ obs_avg(1)
      enddo
      call write_separator_line(22,0,mpi_rank,2,"_")
      call write_simple_line(22,0,mpi_rank,2,"c","Correlated estimation of energy and variance")
      call write_empty_line(22,0,mpi_rank)
      C = w_avg**2 / w2_avg
      if (e_var.gt.10d-10) then
         call write_variable_error_dble_line(22,0,mpi_rank,2,"Energy", e_avg, sqrt(e_var / C ),units="Eh")
         call write_variable_error_dble_line(22,0,mpi_rank,2,"Variance",e_var, e_var*sqrt(2.0_dp/(C-1.0_dp)),units="Eh^2")
      else
         call write_variable_error_dble_line(22,0,mpi_rank,2,"Energy", e_avg, 0.0_dp,units="Eh")
         call write_variable_error_dble_line(22,0,mpi_rank,2,"Variance",0.0_dp, 0.0_dp,units="Eh^2")
      endif
      call write_separator_line(22,0,mpi_rank,2,"_")
      call write_simple_line(22,0,mpi_rank,2,"c","Energy error variation per block size")
      call write_empty_line(22,0,mpi_rank)
      call write_simple_line(22,0,mpi_rank,2,"c","Block length       Error [Eh]       Stnd. error of error [Eh]")
      correct_blck_length = 0
      do i1 = 1_int32, n_blck_sizes
         C = 0.5 / (w_avg**2/blck_var(1,i1)-1.0)
         write(22,'(3X,I13,E22.7,3X,E22.7)') blck_sizes(i1), blck_err(2,i1), blck_err(2,i1) * sqrt( C )
         if (i1.gt.1_int32) then
            diff = (blck_err(2,i1) - blck_err(2,i1-1))
            err1 = blck_err(2,i1) * sqrt( C )
            err2 = blck_err(2,i1-1) * sqrt( 0.5 / (w_avg**2/blck_var(1,i1-1)-1.0) )
            if ( diff.le.sqrt(err1**2+err2**2) ) then
               correct_blck_length = i1
               exit
            endif
         endif
      enddo
      write(22,*)
      if (correct_blck_length.eq.0) correct_blck_length = n_blck_sizes
      i2 = correct_blck_length
      n_blcks = (n_tot_smpl + n_tot_smpl_old ) / blck_sizes(i2)
      i1 = 2_int32
      e_blck_err => blck_err(i1,i2) ; i1 = i1 + 1_int32
      k_blck_err => blck_err(i1,i2) ; i1 = i1 + 1_int32
      v_blck_err => blck_err(i1,i2) ; i1 = i1 + 1_int32
      if ( expec_pote ) then
         if ( n_fe.gt.0_int32 )  then
            v_ee_blck_err => blck_err(i1,i2) ; i1 = i1 + 1_int32
            v_fn_blck_err => blck_err(i1:i1+n_fe-1,i2) ; i1 = i1 + n_fe
         endif
         if ( n_qdo.gt.0_int32 ) then
            v_dd_blck_err => blck_err(i1,i2) ; i1 = i1 + 1_int32
            v_dq_blck_err => blck_err(i1:i1+n_qdo-1,i2) ; i1 = i1 + n_qdo
         endif
      endif
      if (dip_mom) then
         dip_blck_err => blck_err(i1:i1+2,i2) ; i1 = i1 + 3_int32
      endif
      if (quad_mom) then
         qud_blck_err => blck_err(i1:i1+5,i2) ; i1 = i1 + 6_int32
      endif
      if (expec_dist) then
         if (n_fe.gt.0_int32) then
            r_fe_blck_err => blck_err(i1:i1+3*n_fe-1,i2) ; i1 = i1 + 3*n_fe
         endif
         if (n_qdo.gt.0_int32) then
            r_drd_blck_err => blck_err(i1:i1+3*n_qdo-1,i2) ; i1 = i1 + 3*n_qdo
         endif
         if (n_el.ge.2_int32) then
            d_ee_blck_err => blck_err(i1,i2) ; i1 = i1 + 1_int32
         endif
         if (n_po.ge.2_int32) then
            d_pp_blck_err => blck_err(i1,i2) ; i1 = i1 + 1_int32
         endif
         if (n_el.gt.0_int32.and.n_po.gt.0_int32) then
            d_ep_blck_err => blck_err(i1,i2) ; i1 = i1 + 1_int32
         endif
         if (n_qdo.gt.1_int32) then
            d_dd_blck_err => blck_err(i1,i2) ; i1 = i1 + 1_int32
         endif
      endif
      call write_separator_line(22,0,mpi_rank,2,"_")
      write(22,'(2X,"Block len. / Num. of Blocks                   ",I10,"/",I10)') blck_sizes(i2), n_blcks
      call write_empty_line(22,0,mpi_rank)
      call write_variable_error_dble_line(22,0,mpi_rank,2,"Energy", e_avg, e_blck_err,units="Eh")
      call write_variable_error_dble_line(22,0,mpi_rank,2,"Kin.Ene.", k_avg, k_blck_err,units="Eh")
      call write_variable_error_dble_line(22,0,mpi_rank,2,"Pot.Ene.", v_avg, v_blck_err,units="Eh")
      if ( expec_pote ) then
         call write_simple_line(22,0,mpi_rank,2,"l","________________________________________")
         call write_empty_line(22,0,mpi_rank)
         if ( n_fe.gt.0_int32 )  then
            call write_variable_error_dble_line(22,0,mpi_rank,2,"v_ff", v_ee_avg, v_ee_blck_err,units="Eh")
            do i1 = 1_int32, n_fe
               write(variable_name,'("v_fn",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), v_fn_avg(i1), v_fn_blck_err(i1),units="Eh")
            enddo
         endif
         if ( n_qdo.gt.0_int32 ) then
            call write_variable_error_dble_line(22,0,mpi_rank,2,"v_dd", v_dd_avg, v_dd_blck_err,units="Eh")
            do i1 = 1_int32, n_qdo
               write(variable_name,'("v_dq",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), v_dq_avg(i1),v_dq_blck_err(i1),units="Eh")
            enddo
         endif
      endif
      if (dip_mom) then
         call write_simple_line(22,0,mpi_rank,2,"l","________________________________________")
         call write_simple_line(22,0,mpi_rank,2,"l","Dipole moment")
         call write_empty_line(22,0,mpi_rank)
         call write_variable_error_dble_line(22,0,mpi_rank,2,"D_x", dip_avg(1)*au_to_debye, dip_blck_err(1)*au_to_debye,units="D")
         call write_variable_error_dble_line(22,0,mpi_rank,2,"D_y", dip_avg(2)*au_to_debye, dip_blck_err(2)*au_to_debye,units="D")
         call write_variable_error_dble_line(22,0,mpi_rank,2,"D_z", dip_avg(3)*au_to_debye, dip_blck_err(3)*au_to_debye,units="D")
         d_sqr = dsqrt(sum(dip_avg(1:3)**2))
         if (d_sqr.gt.0.0_dp) then
            d_sqr_blck_err = dsqrt(sum((dip_avg(1:3)*dip_blck_err(1:3))**2)) / d_sqr
         else
            d_sqr_blck_err = 0.0_dp
         endif
         call write_variable_error_dble_line(22,0,mpi_rank,2,"|D|", d_sqr*au_to_debye, d_sqr_blck_err*au_to_debye,units="D")
      endif
      if (quad_mom) then
         C = 0.5_dp*au_to_debye*au_to_ang
         call write_simple_line(22,0,mpi_rank,2,"l","________________________________________")
         call write_simple_line(22,0,mpi_rank,2,"l","Non-traceless Quadrupole moment")
         call write_empty_line(22,0,mpi_rank)
         call write_variable_error_dble_line(22,0,mpi_rank,2,"Q_xx", qud_avg(1)*C, qud_blck_err(1)*C,units="D*A")
         call write_variable_error_dble_line(22,0,mpi_rank,2,"Q_yy", qud_avg(2)*C, qud_blck_err(2)*C,units="D*A")
         call write_variable_error_dble_line(22,0,mpi_rank,2,"Q_zz", qud_avg(3)*C, qud_blck_err(3)*C,units="D*A")
         call write_variable_error_dble_line(22,0,mpi_rank,2,"Q_xy", qud_avg(4)*C, qud_blck_err(4)*C,units="D*A")
         call write_variable_error_dble_line(22,0,mpi_rank,2,"Q_xz", qud_avg(5)*C, qud_blck_err(5)*C,units="D*A")
         call write_variable_error_dble_line(22,0,mpi_rank,2,"Q_yz", qud_avg(6)*C, qud_blck_err(6)*C,units="D*A")
         q_iso = dsqrt(2.0_dp / 3.0_dp * sum(qud_avg(1:3)**2) )
         if (q_iso.gt.0.0_dp) then
            q_iso_blck_err = 2.0_dp / 3.0_dp * dsqrt( sum(( qud_avg(1:3)*qud_blck_err(1:3) )**2) ) / q_iso
         else
            q_iso_blck_err = 0.0_dp
         endif
      endif
      if (expec_dist) then
         if (n_el.gt.0_int32) then
            call write_simple_line(22,0,mpi_rank,2,"l","________________________________________")
            call write_simple_line(22,0,mpi_rank,2,"l","Average electronic positions (3 components)")
            call write_empty_line(22,0,mpi_rank)
            do i1 = 1_int32, n_el
               write(variable_name,'("r_e_x",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), r_fe_avg(3*i1-2), r_fe_blck_err(3*i1-2),units="Bohr")
               write(variable_name,'("r_e_y",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), r_fe_avg(3*i1-1), r_fe_blck_err(3*i1-1),units="Bohr")
               write(variable_name,'("r_e_z",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), r_fe_avg(3*i1),   r_fe_blck_err(3*i1),units="Bohr")
            enddo
         endif
         if (n_po.gt.0_int32) then
            call write_simple_line(22,0,mpi_rank,2,"l","________________________________________")
            call write_simple_line(22,0,mpi_rank,2,"l","Average positronic positions (3 components)")
            call write_empty_line(22,0,mpi_rank)
            do i1 = 1_int32, n_po
               write(variable_name,'("r_p_x",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), r_fe_avg(3*i1-2+3*n_el), r_fe_blck_err(3*i1-2+3*n_el),units="Bohr")
               write(variable_name,'("r_p_y",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), r_fe_avg(3*i1-1+3*n_el), r_fe_blck_err(3*i1-1+3*n_el),units="Bohr")
               write(variable_name,'("r_p_z",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), r_fe_avg(3*i1+3*n_el),   r_fe_blck_err(3*i1+3*n_el),units="Bohr")
            enddo
         endif
         if (n_qdo.gt.0_int32) then
            call write_simple_line(22,0,mpi_rank,2,"l","________________________________________")
            call write_simple_line(22,0,mpi_rank,2,"l","Average drudonic positions (3 components)")
            call write_empty_line(22,0,mpi_rank)
            do i1 = 1_int32, n_qdo
               write(variable_name,'("r_d_x",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), r_drd_avg(3*i1-2), r_drd_blck_err(3*i1-2),units="Bohr")
               write(variable_name,'("r_d_y",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), r_drd_avg(3*i1-1), r_drd_blck_err(3*i1-1),units="Bohr")
               write(variable_name,'("r_d_z",I4)') i1
               call write_variable_error_dble_line(22,0,mpi_rank,2,trim(variable_name), r_drd_avg(3*i1),   r_drd_blck_err(3*i1),units="Bohr")
            enddo
         endif
         call write_simple_line(22,0,mpi_rank,2,"l","________________________________________")
         call write_simple_line(22,0,mpi_rank,2,"l","Average particle distances")
         call write_empty_line(22,0,mpi_rank)
         if (n_el.gt.1_int32) then
            call write_variable_error_dble_line(22,0,mpi_rank,2,"r_ee", d_ee_avg, d_ee_blck_err,units="Bohr")
         endif
         if ( n_el.gt.0_int32 .and. n_po.gt.0_int32) then
            call write_variable_error_dble_line(22,0,mpi_rank,2,"r_ep", d_ep_avg, d_ep_blck_err, units="Bohr")
         endif
         if ( n_po.gt.1_int32) then
            call write_variable_error_dble_line(22,0,mpi_rank,2,"r_pp", d_pp_avg, d_pp_blck_err, units="Bohr")
         endif
         if ( n_qdo.gt.1_int32 ) then
            call write_variable_error_dble_line(22,0,mpi_rank,2,"r_dd", d_dd_avg, d_dd_blck_err,units="Bohr")
         endif
      endif
      call write_empty_line(22,0,mpi_rank)
      call write_separator_line(22,0,mpi_rank,2,"=")
      close(22)
   end subroutine write_qmcavg
   subroutine dealloc_qmcsmp
      if ( allocated(obs_b) ) deallocate( obs_b )
      if ( allocated(obs_pp_b) ) deallocate( obs_pp_b )
      if ( allocated(blck_sizes) ) deallocate( blck_sizes )
      if ( allocated(blck_l) ) deallocate( blck_l )
      if ( allocated(obs_avg) ) deallocate( obs_avg )
      if ( allocated(obs_var) ) deallocate( obs_var )
      if ( allocated(obs_pp_avg) ) deallocate( obs_pp_avg )
      if ( allocated(obs_pp_var) ) deallocate( obs_pp_var )
      if ( allocated(blck_avg) ) deallocate( blck_avg )
      if ( allocated(blck_var) ) deallocate( blck_var )
      if ( allocated(blck_err) ) deallocate( blck_err )
   end subroutine dealloc_qmcsmp
end module monte_carlo_averages_m
