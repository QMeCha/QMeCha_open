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
module pfaffian_params_m
   use fortran_kinds_v,          only: dp, int32
   use openmp_mpi_m
   use geminal_symmetries_m,     only: lmbd_ele_t, lmbd_sym_t
   use fermionic_wavefunction_v, only: spin
   use quantum_monte_carlo_v,    only: restart
   use fermionic_orbitals_m,     only: n_forbs, forbs
   use geminal_symmetries_m,     only: conv_indx, cmp_nzr_fee, cmp_sgm_sym
   use molecular_system_v,       only: n_sys_sym
   implicit none
   type, public :: pffcfs_t
      real(dp),         allocatable, dimension(:,:) :: mat_Lmbd
      real(dp),         allocatable, dimension(:)   :: vec_u, vec_d
      real(dp),         allocatable, dimension(:,:) :: mat_Z_u
      real(dp),         allocatable, dimension(:,:) :: mat_Z_d
      integer(int32)                                :: n_nz_Lmbd, n_nz_Z, n_nz_u, n_nz_d
      type(lmbd_ele_t), allocatable, dimension(:)   :: nz_Lmbd, nz_Z, nz_u, nz_d
      integer(int32)                                :: n_sym_Lmbd, n_sym_Z, n_sym_u, n_sym_d
      type(lmbd_sym_t), allocatable, dimension(:)   :: sym_Lmbd, sym_Z, sym_u, sym_d
   end type pffcfs_t
   type(pffcfs_t), public, save, target :: pff_c_e
   type(pffcfs_t), public, save, target :: pff_c_p
   public  :: ini_pffcfs, save_pffcfs_st, save_pffcfs, sym_pffcfs, upd_pffcfs
   private :: cmp_pffcfs_st, read_pffcfs_st, read_pffcfs
contains
   subroutine ini_pffcfs( pff_c, n_fe, n_fe_d, dl_opt, n_par_det, n_par_det_s, fle_head )
      type(pffcfs_t), target, intent(inout) :: pff_c
      integer(int32),         intent(in)    :: n_fe, n_fe_d
      logical,                intent(in)    :: dl_opt
      integer(int32),         intent(inout) :: n_par_det, n_par_det_s
      character(4), optional, intent(in)    :: fle_head
      integer(int32), pointer :: o1, o2
      integer(int32)          :: i1, i2
      allocate(pff_c%mat_Lmbd(1:n_forbs,1:n_forbs)) ; pff_c%mat_Lmbd = 0.0_dp
      if( mod(n_fe,2_int32).ne.0_int32 ) then
         allocate(pff_c%vec_u(1:n_forbs)) ; pff_c%vec_u = 0.0_dp
         allocate(pff_c%vec_d(1:n_forbs)) ; pff_c%vec_d = 0.0_dp
      endif
      allocate( pff_c%mat_Z_u(1:n_forbs,1:n_forbs) ) ; pff_c%mat_Z_u = 0.0_dp
      if ( n_fe_d.ge.2_int32 ) then
         allocate( pff_c%mat_Z_d(1:n_forbs,1:n_forbs) ) ; pff_c%mat_Z_d = 0.0_dp
      endif
      if ( restart ) then
         call read_pffcfs( pff_c, n_fe, n_fe_d, fle_head )
         call read_pffcfs_st( pff_c, fle_head )
      else
         call cmp_pffcfs_st( pff_c, n_fe )
         do i1 = 1_int32, n_forbs
            pff_c%mat_Lmbd(i1,i1) = 1.0_dp 
         enddo
         do i1 = 1_int32, pff_c%n_sym_Z ; do i2 = 1_int32, pff_c%sym_Z(i1)%n_c
               o1 => pff_c%nz_Z(abs(pff_c%sym_Z(i1)%c(i2)))%o1
               o2 => pff_c%nz_Z(abs(pff_c%sym_Z(i1)%c(i2)))%o2
               pff_c%mat_Z_u(o1,o2) = dble(sign(1,pff_c%sym_Z(i1)%c(i2))) * 0.0001_dp
               if( n_fe_d.gt.1_int32 ) pff_c%mat_Z_d(o1,o2) = pff_c%mat_Z_u(o1,o2)
            enddo ; enddo
         if (mod(n_fe,2_int32).gt.0_int32) then
            pff_c%vec_u(1) = 1.0_dp
            pff_c%vec_d(1) = 1.0_dp
         endif
      endif
      if ( dl_opt ) then
         n_par_det   = pff_c%n_nz_Lmbd + pff_c%n_nz_Z
         n_par_det_s = pff_c%n_sym_Lmbd + pff_c%n_sym_Z
         if ( mod(n_fe,2_int32).gt.0_int32 ) then
            n_par_det  = n_par_det + pff_c%n_nz_u
            n_par_det_s = n_par_det_s + pff_c%n_sym_u
            if (spin.eq.'U') then
               n_par_det  = n_par_det + pff_c%n_nz_d
               n_par_det_s = n_par_det_s + pff_c%n_sym_d
            endif
         endif
         if ( ( spin.eq.'U' ) .and. ( n_fe_d.ge.2_int32 ) ) then
            n_par_det  = n_par_det    + pff_c%n_nz_Z
            n_par_det_s = n_par_det_s + pff_c%n_sym_Z
         endif
      else
         n_par_det   = 0_int32
         n_par_det_s = 0_int32
      endif
   end subroutine ini_pffcfs
   subroutine cmp_pffcfs_st( pff_c, n_fe )
      type(pffcfs_t), intent(inout) :: pff_c      
      integer(int32), intent(in)    :: n_fe
      type(lmbd_ele_t), allocatable, dimension(:) :: nz_tmp
      type(lmbd_sym_t), allocatable, dimension(:) :: sym_tmp
      integer(int32) :: i1
      pff_c%n_nz_Lmbd = 0_int32
      pff_c%n_nz_Z    = 0_int32
      pff_c%n_nz_u    = 0_int32
      pff_c%n_nz_d    = 0_int32
      pff_c%n_sym_Lmbd = 0_int32
      pff_c%n_sym_Z    = 0_int32
      pff_c%n_sym_u    = 0_int32
      pff_c%n_sym_d    = 0_int32
      pff_c%n_nz_Lmbd = n_forbs**2
      allocate( nz_tmp(1:pff_c%n_nz_Lmbd) )
      pff_c%n_sym_Lmbd = pff_c%n_nz_Lmbd
      allocate( sym_tmp(1:pff_c%n_sym_Lmbd) )
      do i1 = 1_int32, pff_c%n_sym_Lmbd
         sym_tmp(i1)%is_indp = .true.
         sym_tmp(i1)%n_c = 1_int32
         allocate( sym_tmp(i1)%c(1:2*n_sys_sym) )
         sym_tmp(i1)%c(1) = i1
      enddo
      call cmp_nzr_fee( 'f', n_forbs, forbs, pff_c%n_nz_Lmbd, nz_tmp )
      call cmp_sgm_sym( .true., spin, 's', n_forbs, forbs, pff_c%n_nz_Lmbd, nz_tmp, pff_c%n_sym_Lmbd, sym_tmp )
      allocate( pff_c%nz_Lmbd(1:pff_c%n_nz_Lmbd) ) ; pff_c%nz_Lmbd(1:pff_c%n_nz_Lmbd) = nz_tmp(1:pff_c%n_nz_Lmbd)
      deallocate( nz_tmp )
      allocate( pff_c%sym_Lmbd(1:pff_c%n_sym_Lmbd) )
      do i1 = 1_int32, pff_c%n_sym_Lmbd
         pff_c%sym_Lmbd(i1)%n_c = sym_tmp(i1)%n_c
         allocate( pff_c%sym_Lmbd(i1)%c(1:pff_c%sym_Lmbd(i1)%n_c) )
         pff_c%sym_Lmbd(i1)%c(1:pff_c%sym_Lmbd(i1)%n_c) = sym_tmp(i1)%c(1:pff_c%sym_Lmbd(i1)%n_c)
      enddo
      deallocate( sym_tmp )
      pff_c%n_nz_Z = n_forbs**2
      allocate( nz_tmp(1:pff_c%n_nz_Z) )
      pff_c%n_sym_Z    = pff_c%n_nz_Z
      allocate( sym_tmp(1:pff_c%n_sym_Z) )
      do i1 = 1_int32, pff_c%n_sym_Z
         sym_tmp(i1)%is_indp = .true.
         sym_tmp(i1)%n_c = 1_int32
         allocate( sym_tmp(i1)%c(1:2*n_sys_sym) )
         sym_tmp(i1)%c(1) = i1
      enddo
      call cmp_nzr_fee( 'f', n_forbs, forbs, pff_c%n_nz_Z, nz_tmp )
      call cmp_sgm_sym( .true., 'R', 't', n_forbs, forbs, pff_c%n_nz_Z, nz_tmp, pff_c%n_sym_Z, sym_tmp )
      allocate( pff_c%nz_Z(1:pff_c%n_nz_Z) ) ; pff_c%nz_Z(1:pff_c%n_nz_Z) = nz_tmp(1:pff_c%n_nz_Z)
      deallocate( nz_tmp )
      allocate(pff_c%sym_Z(1:pff_c%n_sym_Z))
      do i1 = 1_int32, pff_c%n_sym_Z
         pff_c%sym_Z(i1)%n_c = sym_tmp(i1)%n_c
        allocate( pff_c%sym_Z(i1)%c(1:pff_c%sym_Z(i1)%n_c) ) 
        pff_c%sym_Z(i1)%c(1:pff_c%sym_Z(i1)%n_c) = sym_tmp(i1)%c(1:pff_c%sym_Z(i1)%n_c)
      enddo
      deallocate( sym_tmp )
      if( mod(n_fe,2_int32).gt.0_int32 ) then
         pff_c%n_nz_u = n_forbs
         allocate( pff_c%nz_u(1:pff_c%n_nz_u) )
         do i1 = 1_int32, n_forbs
            pff_c%nz_u(i1)%o1 = i1
            pff_c%nz_u(i1)%o2 = 1
            pff_c%nz_u(i1)%io = i1
         enddo
         pff_c%n_sym_u = pff_c%n_nz_u
         allocate( pff_c%sym_u(1:pff_c%n_sym_u) )
         do i1 = 1_int32, pff_c%n_sym_u
            pff_c%sym_u(i1)%n_c = 1_int32
            allocate( pff_c%sym_u(i1)%c(1) )
            pff_c%sym_u(i1)%c(1) = i1
         enddo
         if ( spin.eq.'U' ) then
            pff_c%n_nz_d = n_forbs
            allocate( pff_c%nz_d(1:pff_c%n_nz_d) )
            do i1 = 1_int32, n_forbs
               pff_c%nz_d(i1)%o1 = i1
               pff_c%nz_d(i1)%o2 = 1
               pff_c%nz_d(i1)%io = i1
            enddo
            pff_c%n_sym_d = pff_c%n_nz_d
            allocate( pff_c%sym_d(1:pff_c%n_sym_d) )
            do i1 = 1_int32, pff_c%n_sym_d
               pff_c%sym_d(i1)%n_c = 1_int32
               allocate( pff_c%sym_d(i1)%c(1) )
               pff_c%sym_d(i1)%c(1) = i1
            enddo
         endif
      endif
   end subroutine cmp_pffcfs_st
   subroutine read_pffcfs_st( pff_c, fle_head )
      type(pffcfs_t),           intent(inout) :: pff_c
      character(4),             intent(in)    :: fle_head
      integer(int32) :: o1
      integer(int32) :: i1, i2
      integer(int32) :: n_c
      integer(int32), allocatable, dimension(:) :: c
      character(100) :: svf_fle
      if(mpi_rank.eq.0_int32) then
         svf_fle = 'wvfn.save/' // trim(fle_head) // '_s.sav'
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*) !'("# Table of non zero elements of paired Lambda matrix")'
         read(1,*) pff_c%n_nz_Lmbd
         allocate( pff_c%nz_Lmbd(1:pff_c%n_nz_Lmbd) )
         do i1 = 1_int32, pff_c%n_nz_Lmbd
            read(1,*) pff_c%nz_Lmbd(i1)%o1, pff_c%nz_Lmbd(i1)%o2
            pff_c%nz_Lmbd(i1)%io = pff_c%nz_Lmbd(i1)%o1 + n_forbs * ( pff_c%nz_Lmbd(i1)%o2 -1 )
         enddo
         allocate( c(1:2*pff_c%n_nz_Lmbd) )
         read(1,*)  !'("# Symmetry table of paired Lambda matrix")'
         read(1,*) pff_c%n_sym_Lmbd
         allocate(pff_c%sym_Lmbd(1:pff_c%n_sym_Lmbd))
         do i1 = 1_int32, pff_c%n_sym_Lmbd
            read(1,*) n_c, ( c(2*i2-1:2*i2), i2 = 1_int32, n_c )
            pff_c%sym_Lmbd(i1)%is_indp = .true.
            pff_c%sym_Lmbd(i1)%n_c     = n_c
            allocate( pff_c%sym_Lmbd(i1)%c(1:n_c) )
            do i2 = 1_int32, n_c
               o1 = abs(c(2*i2-1)) + n_forbs * ( c(2*i2) - 1 )
               call conv_indx( pff_c%n_nz_Lmbd, pff_c%nz_Lmbd, o1 )
               pff_c%sym_Lmbd(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
            enddo 
         enddo
         deallocate(c)
         read(1,*) !'("# Table of non zero elements of triplet matrix")'
         read(1,*) pff_c%n_nz_Z
         allocate( pff_c%nz_Z(1:pff_c%n_nz_Z) )
         do i1 = 1_int32, pff_c%n_nz_Z
            read(1,*) pff_c%nz_Z(i1)%o1, pff_c%nz_Z(i1)%o2
            pff_c%nz_Z(i1)%io = pff_c%nz_Z(i1)%o1 + n_forbs * ( pff_c%nz_Z(i1)%o2 -1 )
         enddo
         allocate( c(1:2*pff_c%n_nz_Z) )
         read(1,*)  !'("# Symmetry table of triplet matrix")'
         read(1,*) pff_c%n_sym_Z
         allocate(pff_c%sym_Z(1:pff_c%n_sym_Z))
         do i1 = 1_int32, pff_c%n_sym_Z
            read(1,*) n_c, ( c(2*i2-1:2*i2), i2 = 1_int32, n_c )
            pff_c%sym_Z(i1)%is_indp = .true.
            pff_c%sym_Z(i1)%n_c     = n_c
            allocate( pff_c%sym_Z(i1)%c(1:n_c) )
            do i2 = 1_int32, n_c
               o1 = abs(c(2*i2-1)) + n_forbs * ( c(2*i2) - 1 )
               call conv_indx( pff_c%n_nz_Z, pff_c%nz_Z, o1 )
               pff_c%sym_Z(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
            enddo ! i2 (n_c)
         enddo
         deallocate(c)
         pff_c%n_nz_u  = 0_int32
         pff_c%n_sym_u = 0_int32
         pff_c%n_nz_d  = 0_int32
         pff_c%n_sym_d = 0_int32
         read(1,*) !'("# Table of non zero elements spin up unpaired orbital")')
         read(1,*) pff_c%n_nz_u
         allocate( pff_c%nz_u(1:pff_c%n_nz_u) )
         do i1 = 1_int32, pff_c%n_nz_u
            read(1,*) pff_c%nz_u(i1)%o1
            pff_c%nz_u(i1)%io = pff_c%nz_u(i1)%o1
         enddo
         allocate( c(1:pff_c%n_nz_u) )
         read(1,*) !'("# Symmetry table of spin up unpaired orbital")')
         read(1,*) pff_c%n_sym_u
         allocate( pff_c%sym_u(1:pff_c%n_sym_u) )
         do i1 = 1_int32, pff_c%n_sym_u
            read(1,*) n_c, c(1:n_c)
            pff_c%sym_u(i1)%is_indp = .true.
            pff_c%sym_u(i1)%n_c     = n_c
            allocate( pff_c%sym_u(i1)%c(1:n_c) )
            pff_c%sym_u(i1)%c(1:n_c) = c(1:n_c)
         enddo
         deallocate( c )
         if (spin.eq.'U') then
            read(1,*) !'("# Table of non zero elements spin down unpaired orbital")')
            read(1,*) pff_c%n_nz_d
            allocate( pff_c%nz_d(1:pff_c%n_nz_d) )
            do i1 = 1_int32, pff_c%n_nz_d
               read(1,*) pff_c%nz_d(i1)%o1
               pff_c%nz_u(i1)%io = pff_c%nz_d(i1)%o1
            enddo
            allocate( c(1:pff_c%n_nz_d) )
            read(1,*) !'("# Symmetry table of spin down unpaired orbital")')
            read(1,*) pff_c%n_sym_d
            allocate( pff_c%sym_d(1:pff_c%n_sym_d) )
            do i1 = 1_int32, pff_c%n_sym_d
               read(1,*) n_c, c(1:n_c)
               pff_c%sym_d(i1)%is_indp = .true.
               pff_c%sym_d(i1)%n_c     = n_c
               allocate( pff_c%sym_d(i1)%c(1:n_c) )
               pff_c%sym_d(i1)%c(1:n_c) = c(1:n_c)
            enddo
            deallocate( c )
         endif
         close(1)
      endif 
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(pff_c%n_nz_Lmbd, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(pff_c%n_sym_Lmbd,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      if ( mpi_rank.ne.0_int32 ) allocate( pff_c%nz_Lmbd(1:pff_c%n_nz_Lmbd) )
      do i1 = 1_int32, pff_c%n_nz_Lmbd
         call mpi_bcast(pff_c%nz_Lmbd(i1)%o1, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(pff_c%nz_Lmbd(i1)%o2, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(pff_c%nz_Lmbd(i1)%io, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      if ( mpi_rank.ne.0_int32 ) allocate( pff_c%sym_Lmbd(1:pff_c%n_sym_Lmbd) )
      do i1 = 1_int32, pff_c%n_sym_Lmbd
         call mpi_bcast(pff_c%sym_Lmbd(i1)%n_c,     1, MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(pff_c%sym_Lmbd(i1)%is_indp, 1, MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
         if ( mpi_rank.ne.0_int32 ) allocate( pff_c%sym_Lmbd(i1)%c(1:pff_c%sym_Lmbd(i1)%n_c) )
         call mpi_bcast(pff_c%sym_Lmbd(i1)%c(:),pff_c%sym_Lmbd(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      call mpi_bcast(pff_c%n_nz_Z, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(pff_c%n_sym_Z,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      if ( mpi_rank.ne.0_int32 ) allocate( pff_c%nz_Z(1:pff_c%n_nz_Z) )
      do i1 = 1_int32, pff_c%n_nz_Z
         call mpi_bcast(pff_c%nz_Z(i1)%o1, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(pff_c%nz_Z(i1)%o2, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(pff_c%nz_Z(i1)%io, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      if ( mpi_rank.ne.0_int32 ) allocate( pff_c%sym_Z(1:pff_c%n_sym_Z) )
      do i1 = 1_int32, pff_c%n_sym_Z
         call mpi_bcast(pff_c%sym_Z(i1)%n_c,     1, MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(pff_c%sym_Z(i1)%is_indp, 1, MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
         if ( mpi_rank.ne.0_int32 ) allocate( pff_c%sym_Z(i1)%c(1:pff_c%sym_Z(i1)%n_c) )
         call mpi_bcast(pff_c%sym_Z(i1)%c(:),pff_c%sym_Z(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      call mpi_bcast(pff_c%n_nz_u, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(pff_c%n_sym_u,1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(pff_c%n_nz_d, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(pff_c%n_sym_d,1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      if (pff_c%n_nz_u.gt.0_int32) then
         if ( mpi_rank.ne.0_int32 ) allocate( pff_c%nz_u(1:pff_c%n_nz_u) )
         do i1 = 1_int32, pff_c%n_nz_u
            call mpi_bcast(pff_c%nz_u(i1)%o1, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(pff_c%nz_u(i1)%o2, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(pff_c%nz_u(i1)%io, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         enddo
         if ( mpi_rank.ne.0_int32 ) allocate( pff_c%sym_u(1:pff_c%n_sym_u) )
         do i1 = 1_int32, pff_c%n_sym_u
            call mpi_bcast(pff_c%sym_u(i1)%n_c,     1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(pff_c%sym_u(i1)%is_indp, 1_int32,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
            if ( mpi_rank.ne.0_int32 ) allocate( pff_c%sym_u(i1)%c(1:pff_c%sym_u(i1)%n_c) )
            call mpi_bcast(pff_c%sym_u(i1)%c(:),pff_c%sym_u(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         enddo
      endif
      if (spin.eq.'U') then
         if ( mpi_rank.ne.0_int32 ) allocate( pff_c%nz_d(1:pff_c%n_nz_d) )
         do i1 = 1_int32, pff_c%n_nz_d
            call mpi_bcast(pff_c%nz_d(i1)%o1, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(pff_c%nz_d(i1)%o2, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(pff_c%nz_d(i1)%io, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         enddo
         if ( mpi_rank.ne.0_int32 ) allocate( pff_c%sym_d(1:pff_c%n_sym_d) )
         do i1 = 1_int32, pff_c%n_sym_d
            call mpi_bcast(pff_c%sym_d(i1)%n_c,     1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(pff_c%sym_d(i1)%is_indp, 1_int32,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
            if ( mpi_rank.ne.0_int32 ) allocate( pff_c%sym_d(i1)%c(1:pff_c%sym_d(i1)%n_c) )
            call mpi_bcast(pff_c%sym_d(i1)%c(:),pff_c%sym_d(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         enddo
      endif
#endif
   end subroutine read_pffcfs_st
   subroutine read_pffcfs( pff_c, n_fe, n_fe_d, head_fle )
      type(pffcfs_t), intent(inout) :: pff_c
      integer(int32), intent(in)    :: n_fe, n_fe_d
      character(4),   intent(in)    :: head_fle
      integer(int32) :: o1, o2
      integer(int32) :: i1, n_i
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
      if( mpi_rank.eq.0_int32 ) then
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*) !'("# Triplet Geminal coefficients (only the non null elements)")')
         read(1,*) n_i
         do i1 = 1_int32, n_i
            read(1,*) o1, o2, pff_c%mat_Z_u(o1,o2)
            pff_c%mat_Z_u(o2,o1) = - pff_c%mat_Z_u(o1,o2)
         enddo
         if ( n_fe_d.gt.1_int32 ) then
            read(1,*) !'("# Triplet Geminal coefficients (only the non null elements)")')
            read(1,*) n_i
            do i1 = 1_int32, n_i
               read(1,*) o1, o2, pff_c%mat_Z_d(o1,o2)
               pff_c%mat_Z_d(o2,o1)  = - pff_c%mat_Z_d(o1,o2)
            enddo
         endif
         read(1,*) !'("# Singlet Geminal coefficients (only the non null elements)")')
         read(1,*) n_i
         do i1 = 1_int32, n_i
            read(1,*) o1, o2, pff_c%mat_Lmbd(o1,o2)
         enddo
         if ( mod(n_fe,2_int32).gt.0_int32 ) then
            read(1,*) !'("# Unpaired orbitals coefficients spin up")')
            read(1,*) n_i
            do i1 = 1_int32, n_i
               read(1,*) o1, pff_c%vec_u(o1)
            enddo
            read(1,*) !'("# Unpaired orbitals coefficients spin do")')
            read(1,*) n_i
            do i1 = 1_int32, n_i
               read(1,*) o1, pff_c%vec_d(o1)
            enddo
         endif
         close(1)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(pff_c%mat_Lmbd, n_forbs**2, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
      if ( mod(n_fe,2_int32).gt.0_int32 ) then
         call mpi_bcast(pff_c%vec_u, n_forbs, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(pff_c%vec_d, n_forbs, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
      endif
      call mpi_bcast(pff_c%mat_Z_u, n_forbs**2, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
      if (n_fe_d.gt.1_int32) &
      & call mpi_bcast(pff_c%mat_Z_d, n_forbs**2, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
   end subroutine read_pffcfs
   subroutine save_pffcfs_st( pff_c, head_fle )      
      type(pffcfs_t), target, intent(inout) :: pff_c
      character(4),           intent(in)    :: head_fle
      integer(int32), pointer :: o1
      integer(int32)          :: i1, i2
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '_s.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# Table of non zero elements of paired Lambda matrix")')
      write(1,'(I6)') pff_c%n_nz_Lmbd
      do i1 = 1_int32, pff_c%n_nz_Lmbd
         write(1,'(2I6)') pff_c%nz_Lmbd(i1)%o1, pff_c%nz_Lmbd(i1)%o2
      enddo
      write(1,'("# Symmetry table of paired Lambda matrix")')
      write(1,'(I6)') pff_c%n_sym_Lmbd
      do i1 = 1_int32, pff_c%n_sym_Lmbd
         o1 => pff_c%sym_Lmbd(i1)%n_c
         write(1,'(1000I6)') o1,( sign(1,pff_c%sym_Lmbd(i1)%c(i2)) * pff_c%nz_Lmbd(abs(pff_c%sym_Lmbd(i1)%c(i2)))%o1, &
         & pff_c%nz_Lmbd(abs(pff_c%sym_Lmbd(i1)%c(i2)))%o2, i2 = 1_int32, o1 )
      enddo
      write(1,'("# Table of non zero elements of Zeta matrices")')
      write(1,'(I6)') pff_c%n_nz_Z
      do i1 = 1_int32, pff_c%n_nz_Z
         write(1,'(2I6)') pff_c%nz_Z(i1)%o1, pff_c%nz_Z(i1)%o2
      enddo
      write(1,'("# Symmetry table of the Zeta matrices")')
      write(1,'(I6)') pff_c%n_sym_Z
      do i1 = 1_int32, pff_c%n_sym_Z
         o1 => pff_c%sym_Z(i1)%n_c
         write(1,'(1000I6)') o1,( sign(1,pff_c%sym_Z(i1)%c(i2)) * pff_c%nz_Z(abs(pff_c%sym_Z(i1)%c(i2)))%o1, &
         & pff_c%nz_Z(abs(pff_c%sym_Z(i1)%c(i2)))%o2, i2 = 1_int32, o1 )
      enddo
      write(1,'("# Table of non zero elements of the spin up unpaired orbital")')
      write(1,'(I6)') pff_c%n_nz_u
      do i1 = 1_int32, pff_c%n_nz_u
         write(1,'(I6)') pff_c%nz_u(i1)%o1
      enddo
      write(1,'("# Symmetry table of spin up unpaired orbital")')
      write(1,'(I6)') pff_c%n_sym_u
      do i1 = 1_int32, pff_c%n_sym_u
         o1 => pff_c%sym_u(i1)%n_c
         write(1,'(1000I6)') o1,( sign(1,pff_c%sym_u(i1)%c(i2)) * pff_c%nz_u(abs(pff_c%sym_u(i1)%c(i2)))%o1, i2 = 1_int32, o1)
      enddo
      if (spin.eq.'U') then
         write(1,'("# Table of non zero elements of the spin down unpaired orbital")')
         write(1,'(I6)') pff_c%n_nz_d
         do i1 = 1_int32, pff_c%n_nz_d
            write(1,'(I6)') pff_c%nz_d(i1)%o1
         enddo
         write(1,'("# Symmetry table of spin down unpaired orbital")')
         write(1,'(I6)') pff_c%n_sym_d
         do i1 = 1_int32, pff_c%n_sym_d
            o1 => pff_c%sym_d(i1)%n_c
            write(1,'(1000I6)') o1,( sign(1,pff_c%sym_d(i1)%c(i2)) * pff_c%nz_d(abs(pff_c%sym_d(i1)%c(i2)))%o1, i2 = 1_int32, o1)
         enddo
      endif
      close(1)
   end subroutine save_pffcfs_st
   subroutine save_pffcfs( pff_c, n_fe, n_fe_d, head_fle  )      
      type(pffcfs_t), target, intent(inout) :: pff_c
      integer(int32),         intent(in)    :: n_fe
      integer(int32),         intent(in)    :: n_fe_d
      character(4),           intent(in)    :: head_fle
      integer(int32), pointer :: o1, o2
      integer(int32) :: i1
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# Triplet Geminal coefficients (only the non null elements)")')
      write(1,'(I6)') pff_c%n_nz_Z / 2_int32
      do i1 = 1_int32, pff_c%n_nz_Z
         o1 => pff_c%nz_Z(i1)%o1
         o2 => pff_c%nz_Z(i1)%o2
         if (o1.lt.o2) write(1,'(2I6,E18.9)') o1, o2, pff_c%mat_Z_u(o1,o2)
      enddo
      if ( n_fe_d.gt.1_int32 ) then
         write(1,'("# Triplet Geminal coefficients (only the non null elements)")')
         write(1,'(I6)') pff_c%n_nz_Z / 2_int32
         do i1 = 1_int32, pff_c%n_nz_Z
            o1 => pff_c%nz_Z(i1)%o1
            o2 => pff_c%nz_Z(i1)%o2
            if (o1.lt.o2) write(1,'(2I6,E18.9)') o1, o2, pff_c%mat_Z_d(o1,o2)
         enddo
      endif
      write(1,'("# Singlet Geminal coefficients (only the non null elements)")')
      write(1,'(I6)') pff_c%n_nz_Lmbd
      do i1 = 1_int32, pff_c%n_nz_Lmbd
         o1 => pff_c%nz_Lmbd(i1)%o1
         o2 => pff_c%nz_Lmbd(i1)%o2
         write(1,'(2I6,E18.9)') o1, o2, pff_c%mat_Lmbd(o1,o2)
      enddo
      if ( mod(n_fe,2_int32).gt.0_int32 ) then
         write(1,'("# Unpaired orbitals coefficients spin up")')
         write(1,'(I6)') n_forbs
         do i1 = 1_int32, n_forbs
            write(1,'(I6,E18.9)') i1, pff_c%vec_u(i1)
         enddo
         write(1,'("# Unpaired orbitals coefficients spin down")')
         write(1,'(I6)') n_forbs
         do i1 = 1_int32, n_forbs
            write(1,'(I6,E18.9)') i1, pff_c%vec_d(i1)
         enddo
      endif
      close(1)
   end subroutine save_pffcfs
   subroutine upd_pffcfs( pff_c, n_fe, n_fe_d, n_par_det_s, vec_pff_var )      
      type(pffcfs_t), target,           intent(inout) :: pff_c
      integer(int32),                   intent(in)    :: n_fe, n_fe_d
      integer(int32),                   intent(in)    :: n_par_det_s
      real(dp), dimension(n_par_det_s), intent(in)    :: vec_pff_var
      integer(int32), pointer :: o1, o2
      integer(int32)          :: i1, i2, ip
      ip = 1_int32
      do i1 = 1_int32, pff_c%n_sym_Z
         do i2 = 1_int32, pff_c%sym_Z(i1)%n_c
            o1 => pff_c%nz_Z(abs(pff_c%sym_Z(i1)%c(i2)))%o1
            o2 => pff_c%nz_Z(abs(pff_c%sym_Z(i1)%c(i2)))%o2
            pff_c%mat_Z_u(o1,o2) = pff_c%mat_Z_u(o1,o2) + dble(sign(1,pff_c%sym_Z(i1)%c(i2))) * vec_pff_var(ip)
         enddo
         ip = ip + 1_int32
      enddo
      if( n_fe_d.gt.1_int32 ) then
         if ( spin.eq.'R' ) ip = 1_int32
         do i1 = 1_int32, pff_c%n_sym_Z
            do i2 = 1_int32, pff_c%sym_Z(i1)%n_c
               o1 => pff_c%nz_Z(abs(pff_c%sym_Z(i1)%c(i2)))%o1
               o2 => pff_c%nz_Z(abs(pff_c%sym_Z(i1)%c(i2)))%o2
               pff_c%mat_Z_d(o1,o2) = pff_c%mat_Z_d(o1,o2) + dble(sign(1,pff_c%sym_Z(i1)%c(i2))) * vec_pff_var(ip)
            enddo
            ip = ip + 1_int32
         enddo
      endif
      do i1 = 1_int32, pff_c%n_sym_Lmbd
         do i2 = 1_int32, pff_c%sym_Lmbd(i1)%n_c
            o1 => pff_c%nz_Lmbd(abs(pff_c%sym_Lmbd(i1)%c(i2)))%o1
            o2 => pff_c%nz_Lmbd(abs(pff_c%sym_Lmbd(i1)%c(i2)))%o2
            pff_c%mat_Lmbd(o1,o2) = pff_c%mat_Lmbd(o1,o2) + dble(sign(1,pff_c%sym_Lmbd(i1)%c(i2))) * vec_pff_var(ip)
         enddo
         ip = ip + 1_int32
      enddo
      if ( mod(n_fe,2_int32).gt.0_int32 ) then
         do i1 = 1_int32, pff_c%n_sym_u
            do i2 = 1_int32, pff_c%sym_u(i1)%n_c
               o1 => pff_c%nz_u(abs(pff_c%sym_u(i1)%c(i2)))%o1
               pff_c%vec_u(o1) = pff_c%vec_u(o1) + dble(sign(1,pff_c%sym_u(i1)%c(i2))) * vec_pff_var(ip)
            enddo
            ip = ip + 1_int32
         enddo
         if (spin.eq.'U') then
            do i1 = 1_int32, pff_c%n_sym_d
               do i2 = 1_int32, pff_c%sym_d(i1)%n_c
                  o1 => pff_c%nz_d(abs(pff_c%sym_d(i1)%c(i2)))%o1
                  pff_c%vec_d(o1) = pff_c%vec_d(o1) + dble(sign(1,pff_c%sym_d(i1)%c(i2))) * vec_pff_var(ip)
               enddo
               ip = ip + 1_int32
            enddo
         else
            pff_c%vec_d = pff_c%vec_u
         endif
      endif
   end subroutine upd_pffcfs
   subroutine sym_pffcfs( pff_c, n_fe, n_fe_d, n_par_det, n_par_det_s, dp_ln_pff_P, dp_ln_pff_P_s )      
      type(pffcfs_t),                   intent(inout) :: pff_c
      integer(int32),                   intent(in)    :: n_fe, n_fe_d
      integer(int32),                   intent(inout) :: n_par_det, n_par_det_s
      real(dp), dimension(n_par_det),   intent(in)    :: dp_ln_pff_P
      real(dp), dimension(n_par_det_s), intent(inout) :: dp_ln_pff_P_s
      integer(int32) :: i1, i2, ip, ips
      ips = 1_int32
      ip  = 0_int32
      do i1 = 1_int32, pff_c%n_sym_Z
         dp_ln_pff_P_s(ips) = dp_ln_pff_P(ip+abs(pff_c%sym_Z(i1)%c(1)))
         if ( pff_c%sym_Z(i1)%n_c.ge.2_int32 ) then
            do i2 = 2_int32, pff_c%sym_Z(i1)%n_c
               dp_ln_pff_P_s(ips) = dp_ln_pff_P_s(ips) + dble(sign(1,pff_c%sym_Z(i1)%c(i2))) * dp_ln_pff_P(ip+abs(pff_c%sym_Z(i1)%c(i2)))
            enddo
         endif
         ips = ips + 1_int32
      enddo
      ip = ip + pff_c%n_nz_Z
      if ( (spin.eq.'U') .and. (n_fe_d.gt.1_int32) ) then
         do i1 = 1_int32, pff_c%n_sym_Z
            dp_ln_pff_P_s(ips) = dp_ln_pff_P(ip+abs(pff_c%sym_Z(i1)%c(1)) )
            if ( pff_c%sym_Z(i1)%n_c.ge.2_int32 ) then
               do i2 = 2_int32, pff_c%sym_Z(i1)%n_c
                  dp_ln_pff_P_s(ips) = dp_ln_pff_P_s(ips) + dble(sign(1,pff_c%sym_Z(i1)%c(i2))) * dp_ln_pff_P(ip+abs(pff_c%sym_Z(i1)%c(i2)))
               enddo
            endif
            ips = ips + 1_int32
         enddo
         ip = ip + pff_c%n_nz_Z
      endif
      do i1 = 1_int32, pff_c%n_sym_Lmbd
         dp_ln_pff_P_s(ips) = dp_ln_pff_P(ip+abs(pff_c%sym_Lmbd(i1)%c(1)))
         if ( pff_c%sym_Lmbd(i1)%n_c.ge.2_int32 ) then
            do i2 = 2_int32, pff_c%sym_Lmbd(i1)%n_c
               dp_ln_pff_P_s(ips) = dp_ln_pff_P_s(ips) + &
               & dble(sign(1,pff_c%sym_Lmbd(i1)%c(i2))) * dp_ln_pff_P(ip+abs(pff_c%sym_Lmbd(i1)%c(i2)))
            enddo
         endif
         ips = ips + 1_int32
      enddo
      ip = ip + pff_c%n_nz_Lmbd
      if ( mod(n_fe,2_int32).gt.0_int32 ) then
         do i1 = 1_int32, pff_c%n_sym_u
            dp_ln_pff_P_s(ips) = dp_ln_pff_P(ip+abs(pff_c%sym_u(i1)%c(1)))
            if ( pff_c%sym_u(i1)%n_c.ge.2_int32 ) then
               do i2 = 2_int32, pff_c%sym_u(i1)%n_c
                  dp_ln_pff_P_s(ips) = dp_ln_pff_P_s(ips) + &
                  & dble(sign(1,pff_c%sym_u(i1)%c(i2))) * dp_ln_pff_P(ip+abs(pff_c%sym_u(i1)%c(i2)))
               enddo
            endif
            ips = ips + 1_int32
         enddo
         ip = ip + pff_c%n_sym_u
         if (spin.eq.'U') then
            do i1 = 1_int32, pff_c%n_sym_d
               dp_ln_pff_P_s(ips) = dp_ln_pff_P(ip+abs(pff_c%sym_d(i1)%c(1)))
               if ( pff_c%sym_d(i1)%n_c.ge.2_int32 ) then
                  do i2 = 2_int32, pff_c%sym_d(i1)%n_c
                     dp_ln_pff_P_s(ips) = dp_ln_pff_P_s(ips) + &
                     & dble(sign(1,pff_c%sym_d(i1)%c(i2))) * dp_ln_pff_P(ip+abs(pff_c%sym_d(i1)%c(i2)))
                  enddo
               endif
               ips = ips + 1_int32
            enddo
         endif
      endif
   end subroutine sym_pffcfs
end module pfaffian_params_m
