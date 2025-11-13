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
module dipqcfs_mod
   use fortran_kinds_v, only: dp, int32
   use qdo_wavefunction_v, only: ql_opt
   use openmp_mpi_m
   implicit none
   type, public :: dipqcfs_t
      real(dp), allocatable, dimension(:,:) :: mat_A
      real(dp), allocatable, dimension(:)   :: vec_xi
      integer(int32)                              :: n_A
   end type dipqcfs_t
   type(dipqcfs_t),  public, save, target       :: dipq_c_d
   public   :: ini_dipqcfs, upd_dipqcfs, save_dipqcfs
   private  :: read_dipqcfs
contains
   subroutine ini_dipqcfs( dipq_c, n_qdo, n_par_drd, n_par_drd_s, head_fle )
      use qdo_system_m,  only: qdos
      use quantum_monte_carlo_v,  only: restart
      type(dipqcfs_t),         intent(inout) :: dipq_c
      integer(int32),         intent(in)    :: n_qdo
      integer(int32),         intent(inout) :: n_par_drd, n_par_drd_s
      character(4), optional, intent(in)    :: head_fle
      integer(int32)  :: i1, i2
      dipq_c%n_A  = 3_int32*n_qdo*(3_int32*n_qdo+1)
      allocate( dipq_c%mat_A(1:3_int32*n_qdo,1:3_int32*n_qdo) ) ; dipq_c%mat_A = 0.0_dp
      allocate( dipq_c%vec_xi(1:3_int32*n_qdo) ) ; dipq_c%vec_xi = 0.0_dp
      if ( restart ) then
         call read_dipqcfs( dipq_c, n_qdo, head_fle )
      else
         do i1 = 1_int32, n_qdo
            do i2 = 1_int32, 3_int32
               dipq_c%mat_A((i1-1_int32)*3_int32 + i2, (i1-1_int32)*3_int32 + i2) = 0.5_dp * qdos(i1)%qdo_m *  qdos(i1)%qdo_omega
            enddo 
         enddo
      endif
      if (ql_opt) then
         n_par_drd = 3_int32*n_qdo*(3_int32*n_qdo + 3_int32)/2_int32
         n_par_drd_s = 3_int32*n_qdo*(3_int32*n_qdo + 3_int32)/2_int32
      else
         n_par_drd = 0_int32
         n_par_drd_s = 0_int32
      endif
   end subroutine ini_dipqcfs
   subroutine upd_dipqcfs( dipq_c, n_qdo, n_par_drd, vec_dipq_var )
      type(dipqcfs_t),  target, intent(inout) :: dipq_c
      integer(int32),                   intent(in) :: n_qdo
      integer(int32),                   intent(in) :: n_par_drd
      real(dp), dimension(n_par_drd), intent(in) :: vec_dipq_var
      integer(int32) :: i1, i2, i3
      i3 = 0_int32
      do i1 = 1_int32, 3_int32*n_qdo
         do i2 = i1, 3_int32*n_qdo
            i3 = i3 + 1_int32
            dipq_c%mat_A(i1,i2) = dipq_c%mat_A(i1,i2) + vec_dipq_var(i3)
            dipq_c%mat_A(i2,i1) = dipq_c%mat_A(i1,i2)
         enddo ! (n_qdo)
      enddo ! (n_qdo)
      dipq_c%vec_xi(:) = dipq_c%vec_xi(:) + vec_dipq_var(i3+1:i3+3*n_qdo)
   end subroutine upd_dipqcfs
   subroutine save_dipqcfs( dipq_c, n_qdo, head_fle )
      type(dipqcfs_t),  target, intent(inout) :: dipq_c
      integer(int32),         intent(in)    :: n_qdo
      character(4), optional, intent(in)    :: head_fle
      integer(int32) :: i1, j1, i2, j2
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# QDO dip-dip wave function coefficients")')
      do i1 = 1_int32, n_qdo
         do i2 = i1, n_qdo
            write(1,"(a,I3,a,I3)") 'QDO_i= ', i1, '   QDO_j= ', i2
            do j1 = 1_int32, 3_int32
               if( i1.eq.i2 ) then
                  do j2 = j1, 3_int32
                     write(1,*) dipq_c%mat_A((i1-1_int32)*3_int32 + j1, (i2-1_int32)*3_int32 + j2)
                  enddo 
               else
                  do j2 = 1_int32, 3_int32
                     write(1,*) dipq_c%mat_A((i1-1_int32)*3_int32 + j1, (i2-1_int32)*3_int32 + j2)
                  enddo 
               endif
            enddo 
         enddo 
      enddo 
      write(1,'("# QDO dip-dip wave function center shift")')
      do i1 = 1_int32, n_qdo
         write(1,*) dipq_c%vec_xi(3*i1-2:3*i1)
      enddo
      close(1)
   end subroutine save_dipqcfs
   subroutine read_dipqcfs( dipq_c, n_qdo, head_fle )
      type(dipqcfs_t),  target, intent(inout) :: dipq_c
      integer(int32),         intent(in)    :: n_qdo
      character(4), optional, intent(in)    :: head_fle
      integer(int32) :: i1, j1, i2, j2
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
      if( mpi_rank.eq.0_int32 ) then
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*)
         do i1 = 1_int32, n_qdo
            do i2 = i1, n_qdo
               read(1,*)
               do j1 = 1_int32, 3_int32
                  if( i1.eq.i2 ) then
                     do j2 = j1, 3_int32
                        read(1,*) dipq_c%mat_A((i1-1_int32)*3_int32 + j1, (i2-1_int32)*3_int32 + j2)
                        dipq_c%mat_A((i2-1_int32)*3_int32 + j2, (i1-1_int32)*3_int32 + j1) = dipq_c%mat_A((i1-1_int32)*3_int32 + j1, (i2-1_int32)*3_int32 + j2)
                     enddo
                  else
                     do j2 = 1_int32, 3_int32
                        read(1,*) dipq_c%mat_A((i1-1_int32)*3_int32 + j1, (i2-1_int32)*3_int32 + j2)
                        dipq_c%mat_A((i2-1_int32)*3_int32 + j2, (i1-1_int32)*3_int32 + j1) = dipq_c%mat_A((i1-1_int32)*3_int32 + j1, (i2-1_int32)*3_int32 + j2)
                     enddo 
                  endif
               enddo 
            enddo
         enddo 
         read(1,*)
         do i1 = 1_int32, n_qdo
            read(1,*) dipq_c%vec_xi(3*i1-2:3*i1)
         enddo
         close(1)
      endif 
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(dipq_c%mat_A, 3_int32*n_qdo*3_int32*n_qdo, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(dipq_c%vec_xi, 3_int32*n_qdo, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
   end subroutine read_dipqcfs
end module dipqcfs_mod
