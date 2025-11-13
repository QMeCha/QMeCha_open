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
module singlet_geminal_params_m
   use fortran_kinds_v,          only: dp, int32
   use openmp_mpi_m
   use geminal_symmetries_m,     only: lmbd_ele_t, lmbd_sym_t, cmp_nzr_fee, &
   & cmp_sgm_sym, conv_indx
   use molecular_system_v,       only: n_at, n_sys_sym
   use quantum_monte_carlo_v,    only: restart
   use fermionic_wavefunction_v, only: spin
   use fermionic_orbitals_m,     only: n_forbs, forbs
   implicit none
   type, public :: sgmcfs_t
      real(dp),         allocatable, dimension(:,:) :: mat_Lmbd
      integer(int32)                                :: n_nz_Lmbd, n_nz_Lmbd_u
      type(Lmbd_ele_t), allocatable, dimension(:)   :: nz_Lmbd, nz_Lmbd_u
      integer(int32)                                :: n_sym_Lmbd, n_sym_Lmbd_u
      type(Lmbd_sym_t), allocatable, dimension(:)   :: sym_Lmbd, sym_Lmbd_u
   end type sgmcfs_t
   type(sgmcfs_t), public, save, target :: sgm_c_e
   type(sgmcfs_t), public, save, target :: sgm_c_p
   public  :: ini_sgmcfs, save_sgmcfs_st, save_sgmcfs, sym_sgmcfs, upd_sgmcfs
   private :: cmp_sgmcfs_st, cmp_sgmcfu_st, read_sgmcfs_st, read_sgmcfs
contains
   subroutine ini_sgmcfs( sgm_c, n_fe_u, n_fe_d, dl_opt, n_par_det, n_par_det_s, head_fle )
      type(sgmcfs_t),           intent(inout) :: sgm_c
      integer(int32),           intent(in)    :: n_fe_u, n_fe_d
      logical,                  intent(in)    :: dl_opt
      integer(int32),           intent(inout) :: n_par_det, n_par_det_s
      character(4),   optional, intent(in)    :: head_fle
      integer(int32)  :: n_fe_s
      integer(int32)  :: i1
      n_fe_s = n_fe_u - n_fe_d
      allocate( sgm_c%mat_Lmbd(1:n_forbs,1:n_forbs+n_fe_s) ) ; sgm_c%mat_Lmbd = 0.0_dp
      if (restart) then
         call read_sgmcfs( sgm_c, n_fe_s, head_fle )
         call read_sgmcfs_st( sgm_c, n_fe_s, head_fle )
      else
         if (n_at.eq.1_int32) then
            do i1 = 1_int32, n_fe_d
               sgm_c%mat_Lmbd(i1,i1) = 1.0_dp
            enddo 
         else
            do i1 = 1_int32, n_forbs
               sgm_c%mat_Lmbd(i1,i1) = 1.0_dp
            enddo 
         endif
         if ( n_fe_s.gt.0_int32 ) then
            do i1 = 1_int32, n_fe_s
               sgm_c%mat_Lmbd(n_forbs+1-i1,n_forbs+i1) = 1.0_dp / n_forbs
            enddo 
         endif
         call cmp_sgmcfs_st( sgm_c, n_fe_s )
      endif
      if ( dl_opt ) then
         n_par_det = sgm_c%n_nz_Lmbd
         n_par_det_s = sgm_c%n_sym_Lmbd
         if ( n_fe_s.gt.0_int32 ) then
            n_par_det   = n_par_det   + sgm_c%n_nz_Lmbd_u
            n_par_det_s = n_par_det_s + sgm_c%n_sym_Lmbd_u
         endif
      else
         n_par_det   = 0_int32
         n_par_det_s = 0_int32
      endif
   end subroutine ini_sgmcfs
   subroutine cmp_sgmcfs_st( sgm_c, n_fe_s )
      type(sgmcfs_t), intent(inout) :: sgm_c
      integer(int32), intent(in)    :: n_fe_s
      type(lmbd_ele_t), allocatable, dimension(:) :: nz_Lmbd_tmp
      type(lmbd_sym_t), allocatable, dimension(:) :: sym_Lmbd_tmp
      integer(int32) :: i1
      sgm_c%n_nz_Lmbd = n_forbs**2
      allocate(nz_Lmbd_tmp(1:sgm_c%n_nz_Lmbd))
      sgm_c%n_sym_Lmbd = sgm_c%n_nz_Lmbd
      allocate(sym_Lmbd_tmp(1:sgm_c%n_sym_Lmbd))
      do i1 = 1_int32, sgm_c%n_sym_Lmbd
         sym_Lmbd_tmp(i1)%is_indp = .true.
         sym_Lmbd_tmp(i1)%n_c = 1_int32
         allocate( sym_Lmbd_tmp(i1)%c(1:2*n_sys_sym) )
         sym_Lmbd_tmp(i1)%c(1) = i1
      enddo
      call cmp_nzr_fee( 'f', n_forbs, forbs, sgm_c%n_nz_Lmbd, nz_Lmbd_tmp )
      call cmp_sgm_sym( .true., spin, 's', n_forbs, forbs, sgm_c%n_nz_Lmbd, nz_Lmbd_tmp, sgm_c%n_sym_Lmbd, sym_Lmbd_tmp)
      allocate( sgm_c%nz_Lmbd(1:sgm_c%n_nz_Lmbd) ) ; sgm_c%nz_Lmbd(1:sgm_c%n_nz_Lmbd) = nz_Lmbd_tmp(1:sgm_c%n_nz_Lmbd)
      deallocate(nz_Lmbd_tmp)
      allocate(sgm_c%sym_Lmbd(1:sgm_c%n_sym_Lmbd))
      do i1 = 1_int32, sgm_c%n_sym_Lmbd
         sgm_c%sym_Lmbd(i1)%n_c = sym_Lmbd_tmp(i1)%n_c
         allocate( sgm_c%sym_Lmbd(i1)%c(1:sgm_c%sym_Lmbd(i1)%n_c) )
         sgm_c%sym_Lmbd(i1)%c(1:sgm_c%sym_Lmbd(i1)%n_c) = sym_Lmbd_tmp(i1)%c(1:sgm_c%sym_Lmbd(i1)%n_c)
      enddo
      deallocate( sym_Lmbd_tmp )
      if(n_fe_s.gt.0_int32) then
         call cmp_sgmcfu_st( sgm_c, n_fe_s )
      else
         sgm_c%n_nz_Lmbd_u  = 0_int32
         sgm_c%n_sym_Lmbd_u = 0_int32
      endif
   end subroutine cmp_sgmcfs_st
   subroutine cmp_sgmcfu_st( sgm_c, n_fe_s )
      type(sgmcfs_t), intent(inout) :: sgm_c
      integer(int32), intent(in)    :: n_fe_s
      integer(int32) :: i1, i2, i3
      sgm_c%n_nz_Lmbd_u = n_fe_s * n_forbs
      allocate( sgm_c%nz_Lmbd_u(1:sgm_c%n_nz_Lmbd_u) )
      i3 = 1_int32
      do i2 = 1_int32 + n_forbs, n_fe_s + n_forbs ; do i1 = 1_int32, n_forbs
            sgm_c%nz_Lmbd_u(i3)%o1 = i1
            sgm_c%nz_Lmbd_u(i3)%o2 = i2
            sgm_c%nz_Lmbd_u(i3)%io = i1 + n_forbs * (i2 - 1)
            i3 = i3 + 1_int32
         enddo ; enddo
      sgm_c%n_sym_Lmbd_u = sgm_c%n_nz_Lmbd_u
      allocate( sgm_c%sym_Lmbd_u(1:sgm_c%n_sym_Lmbd_u) )
      do i1 = 1_int32, sgm_c%n_sym_Lmbd_u
         sgm_c%sym_Lmbd_u(i1)%n_c = 1_int32
         allocate( sgm_c%sym_Lmbd_u(i1)%c(1) )
         sgm_c%sym_Lmbd_u(i1)%c(1) = i1
      enddo
   end subroutine cmp_sgmcfu_st
   subroutine read_sgmcfs_st( sgm_c, n_fe_s, head_fle )
      type(sgmcfs_t), intent(inout) :: sgm_c
      integer(int32), intent(in)    :: n_fe_s
      character(4),   intent(in)    :: head_fle
      integer(int32) :: o1
      integer(int32) :: i1, i2
      integer(int32) :: n_c
      integer(int32), allocatable, dimension(:) :: c
      character(100) :: svf_fle
      if(mpi_rank.eq.0_int32) then
         svf_fle = 'wvfn.save/' // trim(head_fle) // '_s.sav'
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*) !'("# Table of non zero elements of paired Lambda matrix")'
         read(1,*) sgm_c%n_nz_Lmbd
         allocate( sgm_c%nz_Lmbd(1:sgm_c%n_nz_Lmbd) )
         do i1 = 1_int32, sgm_c%n_nz_Lmbd
            read(1,*) sgm_c%nz_Lmbd(i1)%o1, sgm_c%nz_Lmbd(i1)%o2
            sgm_c%nz_Lmbd(i1)%io = sgm_c%nz_Lmbd(i1)%o1 + n_forbs * ( sgm_c%nz_Lmbd(i1)%o2 -1 )
         enddo
         allocate( c(1:2*sgm_c%n_nz_Lmbd) )
         read(1,*)  !'("# Symmetry table of paired Lambda matrix")'
         read(1,*) sgm_c%n_sym_Lmbd
         allocate(sgm_c%sym_Lmbd(1:sgm_c%n_sym_Lmbd))
         do i1 = 1_int32, sgm_c%n_sym_Lmbd
            read(1,*) n_c, ( c(2*i2-1:2*i2), i2 = 1_int32, n_c )
            sgm_c%sym_Lmbd(i1)%is_indp = .true.
            sgm_c%sym_Lmbd(i1)%n_c     = n_c
            allocate( sgm_c%sym_Lmbd(i1)%c(1:n_c) )
            do i2 = 1_int32, n_c
               o1 = abs(c(2*i2-1)) + n_forbs * ( c(2*i2) - 1 )
               call conv_indx( sgm_c%n_nz_Lmbd, sgm_c%nz_Lmbd, o1 )
               sgm_c%sym_Lmbd(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
            enddo 
         enddo
         deallocate(c)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(sgm_c%n_nz_Lmbd, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(sgm_c%n_sym_Lmbd,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      if ( mpi_rank.ne.0_int32 ) allocate( sgm_c%nz_Lmbd(1:sgm_c%n_nz_Lmbd) )
      do i1 = 1_int32, sgm_c%n_nz_Lmbd
         call mpi_bcast(sgm_c%nz_Lmbd(i1)%o1, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(sgm_c%nz_Lmbd(i1)%o2, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(sgm_c%nz_Lmbd(i1)%io, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      if ( mpi_rank.ne.0_int32 ) allocate( sgm_c%sym_Lmbd(1:sgm_c%n_sym_Lmbd) )
      do i1 = 1_int32, sgm_c%n_sym_Lmbd
         call mpi_bcast(sgm_c%sym_Lmbd(i1)%n_c,     1, MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(sgm_c%sym_Lmbd(i1)%is_indp, 1, MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
         if ( mpi_rank.ne.0_int32 ) allocate( sgm_c%sym_Lmbd(i1)%c(1:sgm_c%sym_Lmbd(i1)%n_c) )
         call mpi_bcast(sgm_c%sym_Lmbd(i1)%c(:),sgm_c%sym_Lmbd(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
#endif
      sgm_c%n_nz_Lmbd_u  = 0_int32
      sgm_c%n_sym_Lmbd_u = 0_int32
      if(mpi_rank.eq.0_int32) then
         if (n_fe_s.gt.0_int32) then
            read(1,*) !'("# Table of non zero elements spin up unpaired orbital")')
            read(1,*) sgm_c%n_nz_Lmbd_u
            allocate( sgm_c%nz_Lmbd_u(1:sgm_c%n_nz_Lmbd_u) )
            do i1 = 1_int32, sgm_c%n_nz_Lmbd_u
               read(1,*) sgm_c%nz_Lmbd_u(i1)%o1, sgm_c%nz_Lmbd_u(i1)%o2
               sgm_c%nz_Lmbd_u(i1)%io = sgm_c%nz_Lmbd_u(i1)%o1 + n_forbs * ( sgm_c%nz_Lmbd_u(i1)%o2 -1 )
            enddo
            allocate( c(1:2*sgm_c%n_nz_Lmbd_u) )
            read(1,*) !'("# Symmetry table of spin up unpaired orbital")')
            read(1,*) sgm_c%n_sym_Lmbd_u
            allocate(sgm_c%sym_Lmbd_u(1:sgm_c%n_sym_Lmbd_u))
            do i1 = 1_int32, sgm_c%n_sym_Lmbd_u
               read(1,*) n_c, ( c(2*i2-1:2*i2), i2 = 1_int32, n_c )
               sgm_c%sym_Lmbd_u(i1)%is_indp = .true.
               sgm_c%sym_Lmbd_u(i1)%n_c     = n_c
               allocate( sgm_c%sym_Lmbd_u(i1)%c(1:n_c) )
               do i2 = 1_int32, n_c
                  o1 = abs(c(2*i2-1)) + n_forbs * ( c(2*i2) - 1 )
                  call conv_indx( sgm_c%n_nz_Lmbd_u, sgm_c%nz_Lmbd_u, o1 )
                  sgm_c%sym_Lmbd_u(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
               enddo 
            enddo
            deallocate( c )
         endif 
         close(1)
      endif 
#if defined _MPI || defined _MPIh || defined _MPI08
      if (n_fe_s.gt.0_int32) then
         call mpi_bcast(sgm_c%n_nz_Lmbd_u, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(sgm_c%n_sym_Lmbd_u,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         if ( mpi_rank.ne.0_int32 ) allocate( sgm_c%nz_Lmbd_u(1:sgm_c%n_nz_Lmbd_u) )
         do i1 = 1_int32, sgm_c%n_nz_Lmbd_u
            call mpi_bcast(sgm_c%nz_Lmbd_u(i1)%o1, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(sgm_c%nz_Lmbd_u(i1)%o2, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(sgm_c%nz_Lmbd_u(i1)%io, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         enddo
         if ( mpi_rank.ne.0_int32 ) allocate( sgm_c%sym_Lmbd_u(1:sgm_c%n_sym_Lmbd_u) )
         do i1 = 1_int32, sgm_c%n_sym_Lmbd_u
            call mpi_bcast(sgm_c%sym_Lmbd_u(i1)%n_c,     1, MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(sgm_c%sym_Lmbd_u(i1)%is_indp, 1, MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
            if ( mpi_rank.ne.0_int32 ) allocate( sgm_c%sym_Lmbd_u(i1)%c(1:sgm_c%sym_Lmbd_u(i1)%n_c) )
            call mpi_bcast(sgm_c%sym_Lmbd_u(i1)%c(:),sgm_c%sym_Lmbd_u(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         enddo
      endif 
#endif
   end subroutine read_sgmcfs_st
   subroutine save_sgmcfs_st( sgm_c, n_fe_s, head_fle )
      type(sgmcfs_t), target, intent(in) :: sgm_c
      integer(int32),         intent(in) :: n_fe_s
      character(4),           intent(in) :: head_fle
      character(100) :: svf_fle
      integer(int32), pointer :: o1
      integer(int32)          :: i1, i2
      svf_fle = 'wvfn.save/'//trim(head_fle)//'_s.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# Table of non zero elements of paired Lambda matrix")')
      write(1,'(I6)') sgm_c%n_nz_Lmbd
      do i1 = 1_int32, sgm_c%n_nz_Lmbd
         write(1,'(2I6)') sgm_c%nz_Lmbd(i1)%o1, sgm_c%nz_Lmbd(i1)%o2
      enddo
      write(1,'("# Symmetry table of paired Lambda matrix")')
      write(1,'(I6)') sgm_c%n_sym_Lmbd
      do i1 = 1_int32, sgm_c%n_sym_Lmbd
         o1 => sgm_c%sym_Lmbd(i1)%n_c
         write(1,'(1000I6)') o1,( sign(1,sgm_c%sym_Lmbd(i1)%c(i2)) * sgm_c%nz_Lmbd(abs(sgm_c%sym_Lmbd(i1)%c(i2)))%o1, &
         & sgm_c%nz_Lmbd(abs(sgm_c%sym_Lmbd(i1)%c(i2)))%o2, i2 = 1_int32, o1 )
      enddo
      if ( n_fe_s.gt.0_int32 ) then
         write(1,'("# Table of non zero elements spin up unpaired orbital")')
         write(1,'(I6)') sgm_c%n_nz_Lmbd_u
         do i1 = 1_int32, sgm_c%n_nz_Lmbd_u
            write(1,'(2I6)') sgm_c%nz_Lmbd_u(i1)%o1, sgm_c%nz_Lmbd_u(i1)%o2
         enddo
         write(1,'("# Symmetry table of spin up unpaired orbital")')
         write(1,'(I6)') sgm_c%n_sym_Lmbd_u
         do i1 = 1_int32, sgm_c%n_sym_Lmbd_u
            o1 => sgm_c%sym_Lmbd_u(i1)%n_c
            write(1,'(1000I6)') o1,( sign(1,sgm_c%sym_Lmbd_u(i1)%c(i2)) * sgm_c%nz_Lmbd_u(abs(sgm_c%sym_Lmbd_u(i1)%c(i2)))%o1, &
            & sgm_c%nz_Lmbd_u(abs(sgm_c%sym_Lmbd_u(i1)%c(i2)))%o2, i2 = 1_int32, o1)
         enddo
      endif
      close(1)
   end subroutine save_sgmcfs_st
   subroutine save_sgmcfs( sgm_c, n_fe_s, head_fle )
      type(sgmcfs_t), target, intent(in) :: sgm_c
      integer(int32),         intent(in) :: n_fe_s
      character(4),           intent(in) :: head_fle
      character(100) :: svf_fle
      integer(int32), pointer :: o1, o2
      integer(int32) :: i1
      svf_fle = 'wvfn.save/'//trim(head_fle)//'.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# Geminal Matrix coefficients (non null elements)")')
      write(1,'(I6)') sgm_c%n_nz_Lmbd
      do i1 = 1_int32, sgm_c%n_nz_Lmbd
         o1 => sgm_c%nz_Lmbd(i1)%o1
         o2 => sgm_c%nz_Lmbd(i1)%o2
         write(1,'(2I6,E18.9)') o1, o2, sgm_c%mat_Lmbd(o1,o2)
      enddo
      if (n_fe_s.gt.0_int32) then
         write(1,'("# Geminal coefficients of the unpaired orbitals")')
         write(1,'(I6)') sgm_c%n_nz_Lmbd_u
         do i1 = 1_int32, sgm_c%n_nz_Lmbd_u
            o1 => sgm_c%nz_Lmbd_u(i1)%o1
            o2 => sgm_c%nz_Lmbd_u(i1)%o2
            write(1,'(2I6,1E18.9)') o1, o2, sgm_c%mat_Lmbd(o1,o2)
         enddo
      endif
      close(1)
   end subroutine save_sgmcfs
   subroutine read_sgmcfs( sgm_c, n_fe_s, head_fle )
      type(sgmcfs_t), intent(inout) :: sgm_c
      integer(int32), intent(in)    :: n_fe_s
      character(4),   intent(in)    :: head_fle
      character(100) :: svf_fle
      integer(int32) :: i1, i_n, o1, o2
      svf_fle = 'wvfn.save/'//trim(head_fle)//'.sav'
      if(mpi_rank.eq.0_int32) then
         sgm_c%mat_Lmbd = 0.0_dp
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*)
         read(1,*) i_n
         do i1 = 1_int32, i_n
            read(1,*) o1, o2, sgm_c%mat_Lmbd(o1, o2)
         enddo
         if (n_fe_s.gt.0_int32) then
            read(1,*)
            read(1,*) i_n
            do i1 = 1_int32, i_n
               read(1,*) o1, o2, sgm_c%mat_Lmbd(o1, o2)
            enddo
         endif
         close(1)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(sgm_c%mat_Lmbd, n_forbs*(n_forbs+n_fe_s), MPI_DOUBLE_PRECISION,&
      & 0,MPI_COMM_WORLD,mpierr)
#endif
   end subroutine read_sgmcfs
   subroutine upd_sgmcfs( sgm_c, n_fe_s, n_par_det_s, vec_sgm_var )
      type(sgmcfs_t), target,           intent(inout) :: sgm_c
      integer(int32),                   intent(in)    :: n_fe_s
      integer(int32),                   intent(in)    :: n_par_det_s
      real(dp), dimension(n_par_det_s), intent(in)    :: vec_sgm_var
      integer(int32), pointer :: o1, o2
      integer(int32)          :: i1, i2
      integer(int32)          :: ip
      real(dp) :: max_lambd
      ip = 1_int32
      do i1 = 1_int32, sgm_c%n_sym_Lmbd
         do i2 = 1_int32, sgm_c%sym_Lmbd(i1)%n_c
            o1 => sgm_c%nz_Lmbd(abs(sgm_c%sym_Lmbd(i1)%c(i2)))%o1
            o2 => sgm_c%nz_Lmbd(abs(sgm_c%sym_Lmbd(i1)%c(i2)))%o2
            sgm_c%mat_Lmbd(o1,o2) = sgm_c%mat_Lmbd(o1,o2) + dble(sign(1,sgm_c%sym_Lmbd(i1)%c(i2))) * vec_sgm_var(ip)
         enddo
         ip = ip + 1_int32
      enddo
      if ( n_fe_s.gt.0_int32 ) then
         do i1 = 1_int32, sgm_c%n_sym_Lmbd_u
            do i2 = 1_int32, sgm_c%sym_Lmbd_u(i1)%n_c
               o1 => sgm_c%nz_Lmbd_u(abs(sgm_c%sym_Lmbd_u(i1)%c(i2)))%o1
               o2 => sgm_c%nz_Lmbd_u(abs(sgm_c%sym_Lmbd_u(i1)%c(i2)))%o2
               sgm_c%mat_Lmbd(o1,o2) = sgm_c%mat_Lmbd(o1,o2) + dble(sign(1,sgm_c%sym_Lmbd_u(i1)%c(i2))) * vec_sgm_var(ip)
            enddo
            ip = ip + 1_int32
         enddo
      endif 
      max_lambd = 0.0_dp
      do i1 = 1_int32, n_forbs
         do i2 = 1_int32, n_forbs
            if ( abs(sgm_c%mat_Lmbd(i2,i1)).gt.max_lambd) max_lambd = abs(sgm_c%mat_Lmbd(i2,i1))
         enddo
      enddo
      if ( n_fe_s.gt.0_int32 ) then
         do i1 = 1_int32, n_fe_s
            do i2 = 1_int32, n_forbs
               if ( abs(sgm_c%mat_Lmbd(i2,n_forbs+i1)).gt.max_lambd) max_lambd = abs(sgm_c%mat_Lmbd(i2,n_forbs+i1))
            enddo
         enddo
      endif
      sgm_c%mat_Lmbd(:,:) = sgm_c%mat_Lmbd(:,:) / max_lambd
   end subroutine upd_sgmcfs
   subroutine sym_sgmcfs( sgm_c, n_fe_s, n_par_det, n_par_det_s, dp_ln_det_G, dp_ln_det_G_s )
      type(sgmcfs_t),                   intent(in)    :: sgm_c
      integer(int32),                   intent(in)    :: n_fe_s
      integer(int32),                   intent(in)    :: n_par_det, n_par_det_s
      real(dp), dimension(n_par_det),   intent(in)    :: dp_ln_det_G
      real(dp), dimension(n_par_det_s), intent(inout) :: dp_ln_det_G_s
      integer(int32) :: i1, i2, i3, i4
      i3 = 1_int32
      i4 = 1_int32
      do i1 = 1_int32, sgm_c%n_sym_Lmbd
         dp_ln_det_G_s(i3) = dp_ln_det_G(abs(sgm_c%sym_Lmbd(i1)%c(1)))
         if ( sgm_c%sym_Lmbd(i1)%n_c.ge.2_int32 ) then
            do i2 = 2_int32, sgm_c%sym_Lmbd(i1)%n_c
               dp_ln_det_G_s(i3) = dp_ln_det_G_s(i3) + &
               & dble(sign(1,sgm_c%sym_Lmbd(i1)%c(i2))) * dp_ln_det_G(abs(sgm_c%sym_Lmbd(i1)%c(i2)))
            enddo
         endif
         i3 = i3 + 1_int32
      enddo
      i4 = i4 + sgm_c%n_nz_Lmbd
      if ( n_fe_s.gt.0_int32 ) then
         do i1 = 1_int32, sgm_c%n_sym_Lmbd_u
            dp_ln_det_G_s(i3) = dp_ln_det_G(i4-1+abs(sgm_c%sym_Lmbd_u(i1)%c(1)))
            if ( sgm_c%sym_Lmbd_u(i1)%n_c.ge.2_int32 ) then
               do i2 = 2_int32, sgm_c%sym_Lmbd_u(i1)%n_c
                  dp_ln_det_G_s(i3) = dp_ln_det_G_s(i3) + &
                  & dble(sign(1,sgm_c%sym_Lmbd_u(i1)%c(i2))) * dp_ln_det_G(i4-1+abs(sgm_c%sym_Lmbd_u(i1)%c(i2)))
               enddo
            endif
            i3 = i3 + 1_int32
         enddo
      endif
   end subroutine sym_sgmcfs
end module singlet_geminal_params_m
