module mersenne_twister19937_m
   use fortran_kinds_v,       only: stdout, int32, dp
   use openmp_mpi_m
   use write_lines_m,         only: write_separator_line, write_simple_line, &
   & write_empty_line, write_variable_line, write_line_check, write_done
   implicit none
   integer(int32), intrinsic :: ishft, ior, ieor, iand
   real(dp),   parameter :: twopi  = 6.283185307179586476925286766559005768394_dp
   integer(int32),  save, public :: rndm_seed
   integer(int32), parameter :: n  = 624_int32
   integer(int32), parameter :: n1 = n + 1_int32
   integer(int32), parameter :: m      = 397_int32
   integer(int32), parameter :: mata   = -1727483681_int32
   integer(int32), parameter, dimension(0:1)   :: mag01 = [ 0 , mata ]
   integer(int32), parameter :: lmask  =  2147483647_int32
   integer(int32), parameter :: umask  = -lmask - 1_int32
   integer(int32), parameter :: tmaskb = -1658038656_int32
   integer(int32), parameter :: tmaskc = -272236544_int32
   type, public :: mt19937_t
      integer(int32) :: seed
      integer(int32) :: mti
      integer(int32), dimension(0:n-1) :: mt
      logical        :: spr_num
      real(dp)   :: spr_rnd_num
   contains
      procedure :: ini    => init_mt19937_rg
      procedure :: rndi   => grndi
      procedure :: rndc   => grndc
      procedure :: rndo   => grndo
      procedure :: rndoo  => grndoo
      procedure :: rndbm  => box_muller
      procedure :: rndexp => exp_dist
   end type mt19937_t
   class(mt19937_t), save, public, allocatable, dimension(:), target  :: rndn
   public :: init_mt19937, save_mt19937, load_mt19937
   private :: grndi, grndc, grndo, grndoo, box_muller,exp_dist, init_mt19937_rg
contains
   subroutine init_mt19937( n_wlk, n_wlk_max )
      integer(int32), intent(in) :: n_wlk
      integer(int32), intent(in) :: n_wlk_max
      integer(int32) :: iw
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","INITIALIZING MERSENNE TWYSTER PSEUDO RAND. NUM. GEN. MT19937")
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Initial random seed for master task", rndm_seed,var_name='rndm_seed')
      allocate( rndn(1:n_wlk_max) )
      rndm_seed = rndm_seed + n_wlk * mpi_rank
!$omp parallel default(shared) private(iw)
!$omp do
      do iw = 1_int32, n_wlk
         call rndn(iw)%ini( rndm_seed + iw - 1 )
      enddo
!$omp end do
!$omp end parallel
      if (n_wlk_max.gt.n_wlk) then
         rndm_seed = rndm_seed + (n_wlk_max - 2*n_wlk) * mpi_rank + n_wlk * n_mpi_tasks
!$omp parallel default(shared) private(iw)
!$omp do
         do iw = n_wlk+1_int32, n_wlk_max
            call rndn(iw)%ini( rndm_seed + iw - 1 )
         enddo
!$omp end do
!$omp end parallel
      endif
   end subroutine init_mt19937
   subroutine save_mt19937( n_wlk_max )
      integer(int32), intent(in) :: n_wlk_max
      character(50) :: svf_fle
      integer(int32) :: mti_buf
      logical        :: spr_num_buf
      real(dp)       :: spr_rnd_num_buf
      integer(int32), allocatable, dimension(:) :: mt_buf
      integer(int32) :: i_mpi_task, iw, icom
      svf_fle = 'qmc.save/random_numbers.sav'
      if (mpi_rank.eq.0_int32) then
         allocate( mt_buf(0:n-1) ) ; mt_buf = 0_int32
         open(unit=1,file=svf_fle,action='write',form='unformatted',status='unknown', access='sequential')
         do iw = 1_int32, n_wlk_max
            write(1) rndn(iw)%mti
            write(1) rndn(iw)%mt(0:n-1)
            write(1) rndn(iw)%spr_num
            write(1) rndn(iw)%spr_rnd_num
         enddo
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      if (n_mpi_tasks.gt.1_int32) then
         do i_mpi_task = 1_int32, n_mpi_tasks-1
            do iw = 1_int32, n_wlk_max
               icom = n_wlk_max * (i_mpi_task-1) + iw
               if (mpi_rank.eq.i_mpi_task) then
                  call mpi_send (rndn(iw)%mti, 1, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, mpierr)
               else if (mpi_rank.eq.0) then
                  call mpi_recv (mti_buf, 1, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
               endif
               call mpi_barrier(MPI_COMM_WORLD, mpierr)
               if (mpi_rank.eq.0_int32) write(1) mti_buf
               if (mpi_rank.eq.i_mpi_task) then
                  call mpi_send (rndn(iw)%mt(0:n-1), n, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, mpierr)
               else if (mpi_rank.eq.0) then
                  call mpi_recv (mt_buf(0:n-1), n, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
               endif
               call mpi_barrier(MPI_COMM_WORLD, mpierr)
               if (mpi_rank.eq.0_int32) write(1) mt_buf(0:n-1)
               if (mpi_rank.eq.i_mpi_task) then
                  call mpi_send (rndn(iw)%spr_num, 1, MPI_LOGICAL, 0, icom, MPI_COMM_WORLD, mpierr)
               else if (mpi_rank.eq.0) then
                  call mpi_recv (spr_num_buf, 1, MPI_LOGICAL, i_mpi_task, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
               endif
               call mpi_barrier(MPI_COMM_WORLD, mpierr)
               if (mpi_rank.eq.0_int32) write(1) spr_num_buf
               if (mpi_rank.eq.i_mpi_task) then
                  call mpi_send (rndn(iw)%spr_rnd_num, 1, MPI_DOUBLE_PRECISION, 0, icom, MPI_COMM_WORLD, mpierr)
               else if (mpi_rank.eq.0) then
                  call mpi_recv (spr_rnd_num_buf, 1, MPI_DOUBLE_PRECISION, i_mpi_task, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
               endif
               call mpi_barrier(MPI_COMM_WORLD, mpierr)
               if (mpi_rank.eq.0_int32) write(1) spr_rnd_num_buf
            enddo
         enddo
      endif
#endif
      if (mpi_rank.eq.0_int32) then
         deallocate(mt_buf)
         close(1)
      endif
   end subroutine save_mt19937
   subroutine load_mt19937(n_wlk_max )
      integer(int32), intent(in) :: n_wlk_max
      character(50) :: svf_fle
      integer(int32) :: mti_buf
      logical        :: spr_num_buf, file_present
      real(dp)       :: spr_rnd_num_buf
      integer(int32), allocatable, dimension(:) :: mt_buf
      integer(int32) :: i_mpi_task, iw, icom
      file_present = .false.
      if (mpi_rank.eq.0_int32) then
         svf_fle = 'qmc.save/random_numbers.sav'
         inquire(file=svf_fle, exist=file_present)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(file_present, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif
      if (file_present) then
         call write_line_check(stdout,0,mpi_rank,2,"Loading random number generation parameters")
         if (mpi_rank.eq.0_int32) then
            allocate( mt_buf(0:n-1) ) ; mt_buf = 0_int32
            open(unit=1,file=svf_fle,action='read',form='unformatted',status='old', access='sequential')
            do iw = 1_int32, n_wlk_max
               read(1) rndn(iw)%mti
               read(1) rndn(iw)%mt(0:n-1)
               read(1) rndn(iw)%spr_num
               read(1) rndn(iw)%spr_rnd_num
            enddo
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         if (n_mpi_tasks.gt.1_int32) then
            do i_mpi_task = 1_int32, n_mpi_tasks-1
               do iw = 1_int32, n_wlk_max
                  icom = n_wlk_max * (i_mpi_task-1) + iw
                  if (mpi_rank.eq.0) then
                     read(1) mti_buf
                     call mpi_send (mti_buf, 1, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, mpierr)
                  elseif ( mpi_rank.eq.i_mpi_task ) then
                     call mpi_recv (rndn(iw)%mti, 1, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
                  endif
                  call mpi_barrier(MPI_COMM_WORLD, mpierr)
                  if (mpi_rank.eq.0) then
                     read(1) mt_buf(0:n-1)
                     call mpi_send (mt_buf(0:n-1), n, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, mpierr)
                  elseif ( mpi_rank.eq.i_mpi_task ) then
                     call mpi_recv (rndn(iw)%mt(0:n-1), n, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
                  endif
                  call mpi_barrier(MPI_COMM_WORLD, mpierr)
                  if (mpi_rank.eq.0) then
                     read(1) spr_num_buf
                     call mpi_send (spr_num_buf, 1, MPI_LOGICAL, i_mpi_task, icom, MPI_COMM_WORLD, mpierr)
                  elseif ( mpi_rank.eq.i_mpi_task ) then
                     call mpi_recv (rndn(iw)%spr_num, 1, MPI_LOGICAL, 0, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
                  endif
                  call mpi_barrier(MPI_COMM_WORLD, mpierr)
                  if (mpi_rank.eq.0) then
                     read(1) spr_rnd_num_buf
                     call mpi_send (spr_rnd_num_buf, 1, MPI_DOUBLE_PRECISION, i_mpi_task, icom, MPI_COMM_WORLD, mpierr)
                  elseif ( mpi_rank.eq.i_mpi_task ) then
                     call mpi_recv (rndn(iw)%spr_rnd_num, 1, MPI_DOUBLE_PRECISION, 0, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
                  endif
                  call mpi_barrier(MPI_COMM_WORLD, mpierr)
               enddo
            enddo
         endif
#endif
         call write_done(stdout,0,mpi_rank)
         if (mpi_rank.eq.0_int32) then
            deallocate(mt_buf)
            close(1)
         endif
      endif
   end subroutine load_mt19937
   subroutine init_mt19937_rg( obj, rseed )
      class(mt19937_t), intent(inout) :: obj
      integer(int32),   intent(in) :: rseed
      integer(int32) :: i1
      obj%seed = rseed
      obj%mt(0) = iand(obj%seed,-1)
      do i1 = 1, n - 1
         obj%mt(i1) = iand( 69069 * obj%mt(i1-1), -1 )
      enddo
      obj%mti = n1
      obj%spr_num = .false.
   end subroutine init_mt19937_rg
   function grndi( obj ) result(y)
      class(mt19937_t), intent(inout) :: obj
      integer(int32) :: y
      integer(int32) :: i1
      if(obj%mti.ge.n) then
         do i1 = 0, n-m-1
            y      = ior(iand(obj%mt(i1),umask),iand(obj%mt(i1+1),lmask))
            obj%mt(i1) = ieor(ieor(obj%mt(i1+m),ishft(y,-1)),mag01(iand(y,1)))
         enddo
         do i1 = n-m, n-2
            y      = ior(iand(obj%mt(i1),umask),iand(obj%mt(i1+1),lmask))
            obj%mt(i1) = ieor(ieor(obj%mt(i1+(m-n)),ishft(y,-1)),mag01(iand(y,1)))
         enddo
         y        = ior(iand(obj%mt(n-1),umask),iand(obj%mt(0),lmask))
         obj%mt(n-1) = ieor(ieor(obj%mt(m-1),ishft(y,-1)),mag01(iand(y,1)))
         obj%mti = 0
      endif
      y = obj%mt(obj%mti)
      obj%mti = obj%mti + 1
      y = ieor(y,ishft(y,-11))
      y = ieor(y,iand(ishft(y,  7),tmaskb))
      y = ieor(y,iand(ishft(y, 15),tmaskc))
      y = ieor(y,ishft(y,-18))
   end function grndi
   function grndo( obj ) result(g)
      class(mt19937_t), intent(inout) :: obj
      integer(int32) :: y
      real(dp)   :: g
      y = obj%rndi()
      if(y .lt. 0) then
         g = ( dble(y) + 2.0**32 ) / ( 2.0**32 )
      else
         g = dble(y) /( 2.0**32 )
      endif
   end function grndo
   function grndoo( obj ) result(g)
      class(mt19937_t), intent(inout) :: obj
      integer(int32) :: y
      real(dp)   :: g
      y = obj%rndi()
      if(y .lt. 0) then
         g = ( dble(y) + 2.0**32 ) / ( 2.0**32 )
      elseif(y.gt.0) then
         g = dble(y) /( 2.0**32 )
      endif
   end function grndoo
   function grndc( obj ) result(g)
      class(mt19937_t), intent(inout) :: obj
      integer(int32) :: y
      real(dp)   :: g
      y = obj%rndi()
      if(y .lt. 0) then
         g = ( dble(y) + 2.0**32 ) / ( 2.0**32 - 1.0 )
      else
         g = dble(y) / ( 2.0**32 - 1.0 )
      endif
   end function grndc
   function box_muller ( obj ) result(g)
      class(mt19937_t), intent(inout) :: obj
      real(dp)   :: u, v
      real(dp)   :: g
      if ( obj%spr_num ) then
         g = obj%spr_rnd_num
         obj%spr_num = .false.
      else
         u = dsqrt( -2.0_dp * dlog( 1.0_dp-obj%rndo() ) )
         v = twopi*obj%rndo()
         g               = u * dcos( v )
         obj%spr_rnd_num = u * dsin( v )
         obj%spr_num = .true.
      endif
   end function box_muller
   function exp_dist ( obj ) result(g)
      class(mt19937_t), intent(inout) :: obj
      real(dp)   :: u, g
      u = obj%rndo()
      u = u * obj%rndo()
      u = u * obj%rndo()
      g = -dlog (u)
   end function exp_dist
end module mersenne_twister19937_m
