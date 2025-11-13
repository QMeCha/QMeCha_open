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
module triplet_geminal_params_m
   use fortran_kinds_v,          only: dp, int32
   use openmp_mpi_m
   use molecular_system_v,       only: n_sys_sym
   use geminal_symmetries_m,     only: lmbd_ele_t, lmbd_sym_t, conv_indx, cmp_nzr_fee, cmp_sgm_sym
   use fermionic_wavefunction_v, only: spin
   use quantum_monte_carlo_v,    only: restart
   use fermionic_orbitals_m,     only: n_forbs, forbs
   implicit none
   type, public :: tgmcfs_t
      real(dp), allocatable, dimension(:,:)       :: mat_Zu
      real(dp), allocatable, dimension(:,:)       :: mat_Zd
      integer(int32)                              :: d_Zu
      integer(int32)                              :: d_Zd
      integer(int32)                              :: n_nz_Z, n_nz_Zu_u, n_nz_Zd_u
      type(Lmbd_ele_t), allocatable, dimension(:) :: nz_Z, nz_Zu_u, nz_Zd_u
      integer(int32)                              :: n_sym_Z, n_sym_Zu_u, n_sym_Zd_u
      type(Lmbd_sym_t), allocatable, dimension(:) :: sym_Z, sym_Zu_u, sym_Zd_u
   end type tgmcfs_t
   type(tgmcfs_t), public, save, target :: tgm_c_e
   type(tgmcfs_t), public, save, target :: tgm_c_p
   public  :: ini_tgmcfs, save_tgmcfs_st, save_tgmcfs, sym_tgmcfs, upd_tgmcfs
   private :: cmp_tgmcfs_st, cmp_tgmcfu_st, read_tgmcfu_st, read_tgmcfs
contains
   subroutine ini_tgmcfs( tgm_c, n_fe_u, n_fe_d, dl_opt, n_par_det, n_par_det_s, head_fle )
      type(tgmcfs_t), target, intent(inout) :: tgm_c
      integer(int32),         intent(in)    :: n_fe_u
      integer(int32),         intent(in)    :: n_fe_d
      logical,                intent(in)    :: dl_opt
      integer(int32),         intent(inout) :: n_par_det
      integer(int32),         intent(inout) :: n_par_det_s
      character(4), optional, intent(in)    :: head_fle
      integer(int32), pointer :: o1, o2
      integer(int32) :: i1, i2
      if( mod(n_fe_u,2_int32).gt.0_int32 ) then
         tgm_c%d_Zu = n_forbs + 1_int32
      else
         tgm_c%d_Zu = n_forbs
      endif
      if(n_fe_d.gt.0_int32) then
         if (n_fe_d.eq.1_int32) then
            tgm_c%d_Zd = 1_int32
         else
            if( mod(n_fe_d,2_int32).gt.0_int32 ) then
               tgm_c%d_Zd = n_forbs + 1_int32
            else
               tgm_c%d_Zd = n_forbs
            endif
         endif
      else
         tgm_c%d_Zd = 0_int32
      endif
      allocate( tgm_c%mat_Zu(1:n_forbs,1:tgm_c%d_Zu) ) ; tgm_c%mat_Zu = 0.0_dp
      if (n_fe_d.gt.0_int32) then
         allocate( tgm_c%mat_Zd(1:n_forbs,1:tgm_c%d_Zd) )  ; tgm_c%mat_Zd = 0.0_dp
      endif
      if (restart) then
         call read_tgmcfs( tgm_c, n_fe_d, head_fle )
         call read_tgmcfu_st( tgm_c, head_fle )
      else
         call cmp_tgmcfs_st( tgm_c, head_fle )
         do i1 = 1_int32, tgm_c%n_sym_Z
            do i2 = 1_int32, tgm_c%sym_Z(i1)%n_c
               o1 => tgm_c%nz_Z(abs(tgm_c%sym_Z(i1)%c(i2)))%o1
               o2 => tgm_c%nz_Z(abs(tgm_c%sym_Z(i1)%c(i2)))%o2
               tgm_c%mat_Zu(o1,o2) = dble(sign(1,tgm_c%sym_Z(i1)%c(i2))) / dble(n_forbs**2)
            enddo
         enddo
         if (tgm_c%n_nz_Zu_u.gt.0_int32) then
            do i1 = 1_int32, tgm_c%n_sym_Zu_u
               do i2 = 1_int32, tgm_c%sym_Zu_u(i1)%n_c
                  o1 => tgm_c%nz_Zu_u(abs(tgm_c%sym_Zu_u(i1)%c(i2)))%o1
                  o2 => tgm_c%nz_Zu_u(abs(tgm_c%sym_Zu_u(i1)%c(i2)))%o2
                  tgm_c%mat_Zu(o1,o2) = dble(sign(1,tgm_c%sym_Zu_u(i1)%c(i2))) / dble(n_forbs)
               enddo
            enddo
         endif
         if ( n_fe_d.gt.0_int32 ) then
            tgm_c%mat_Zd(:,1:n_forbs) = tgm_c%mat_Zu(:,1:n_forbs)
            if (tgm_c%n_nz_Zd_u.gt.0_int32) then
               do i1 = 1_int32, tgm_c%n_sym_Zd_u
                  do i2 = 1_int32, tgm_c%sym_Zd_u(i1)%n_c
                     o1 => tgm_c%nz_Zd_u(abs(tgm_c%sym_Zd_u(i1)%c(i2)))%o1
                     o2 => tgm_c%nz_Zd_u(abs(tgm_c%sym_Zd_u(i1)%c(i2)))%o2
                     tgm_c%mat_Zd(o1,o2) = dble(sign(1,tgm_c%sym_Zd_u(i1)%c(i2))) / dble(n_forbs)
                  enddo
               enddo
            endif
         endif
      endif
      if ( dl_opt ) then
         n_par_det   = tgm_c%n_nz_Z + tgm_c%n_nz_Zu_u
         n_par_det_s = tgm_c%n_sym_Z + tgm_c%n_sym_Zu_u
         if( n_fe_d.gt.0_int32 ) then
            if(spin.eq.'U') then
               n_par_det   = n_par_det   + tgm_c%n_nz_Z + tgm_c%n_nz_Zd_u
               n_par_det_s = n_par_det_s + tgm_c%n_sym_Z + tgm_c%n_sym_Zd_u
            else
               if ( tgm_c%d_Zu.eq.n_forbs .and. tgm_c%d_Zd.ne.n_forbs ) then
                  n_par_det   = n_par_det   + tgm_c%n_nz_Zd_u
                  n_par_det_s = n_par_det_s + tgm_c%n_sym_Zd_u
               endif
            endif
         endif
      else
         n_par_det   = 0_int32
         n_par_det_s = 0_int32
      endif
   end subroutine ini_tgmcfs
   subroutine cmp_tgmcfs_st( tgm_c, head_fle )
      type(tgmcfs_t),           intent(inout) :: tgm_c
      character(4),   optional, intent(in)    :: head_fle
      type(lmbd_ele_t), allocatable, dimension(:) :: nz_Z_tmp
      type(lmbd_sym_t), allocatable, dimension(:) :: sym_Z_tmp
      integer(int32) :: i1
      tgm_c%n_nz_Z = n_forbs **2
      allocate(nz_Z_tmp(1:tgm_c%n_nz_Z))
      tgm_c%n_sym_Z = tgm_c%n_nz_Z
      allocate(sym_Z_tmp(1:tgm_c%n_sym_Z))
      do i1 = 1_int32, tgm_c%n_sym_Z
         sym_Z_tmp(i1)%is_indp = .true.
         sym_Z_tmp(i1)%n_c = 1_int32
         allocate( sym_Z_tmp(i1)%c(1:2*n_sys_sym) )
         sym_Z_tmp(i1)%c(1) = i1
      enddo
      call cmp_nzr_fee( 'f', n_forbs, forbs, tgm_c%n_nz_Z, nz_Z_tmp )
      call cmp_sgm_sym( .true., 'R', 't', n_forbs, forbs, tgm_c%n_nz_Z, nz_Z_tmp, tgm_c%n_sym_Z, sym_Z_tmp)
      allocate(tgm_c%nz_Z(1:tgm_c%n_nz_Z)) ; tgm_c%nz_Z(1:tgm_c%n_nz_Z) = nz_Z_tmp(1:tgm_c%n_nz_Z)
      deallocate(nz_Z_tmp)
      allocate(tgm_c%sym_Z(1:tgm_c%n_sym_Z))
      do i1 = 1_int32, tgm_c%n_sym_Z
         tgm_c%sym_Z(i1)%n_c = sym_Z_tmp(i1)%n_c
         allocate( tgm_c%sym_Z(i1)%c(1:tgm_c%sym_Z(i1)%n_c) )
         tgm_c%sym_Z(i1)%c(1:tgm_c%sym_Z(i1)%n_c) = sym_Z_tmp(i1)%c(1:tgm_c%sym_Z(i1)%n_c)
      enddo
      deallocate( sym_Z_tmp )
      if ( tgm_c%d_Zu.gt.n_forbs .or. tgm_c%d_Zd.gt.n_forbs ) then
         if (restart) then
            call read_tgmcfu_st( tgm_c, head_fle )
         else
            call cmp_tgmcfu_st( tgm_c )
         endif
      else
         tgm_c%n_nz_Zd_u = 0_int32
         tgm_c%n_nz_Zu_u = 0_int32
         tgm_c%n_sym_Zd_u = 0_int32
         tgm_c%n_sym_Zu_u = 0_int32
      endif
   end subroutine cmp_tgmcfs_st
   subroutine cmp_tgmcfu_st( tgm_c )
      type(tgmcfs_t), intent(inout) :: tgm_c
      integer(int32) :: i1
      if ( tgm_c%d_Zu .gt. n_forbs) then
         tgm_c%n_nz_Zu_u = n_forbs
         allocate(tgm_c%nz_Zu_u(1:tgm_c%n_nz_Zu_u))
         do i1 = 1_int32, n_forbs
            tgm_c%nz_Zu_u(i1)%o1 = i1
            tgm_c%nz_Zu_u(i1)%o2 = tgm_c%d_Zu
            tgm_c%nz_Zu_u(i1)%io = i1
         enddo
         tgm_c%n_sym_Zu_u = n_forbs
         allocate(tgm_c%sym_Zu_u(1:tgm_c%n_sym_Zu_u))
         do i1 = 1_int32, tgm_c%n_sym_Zu_u
            tgm_c%sym_Zu_u(i1)%is_indp = .true.
            tgm_c%sym_Zu_u(i1)%n_c = 1_int32
            allocate( tgm_c%sym_Zu_u(i1)%c(1) )
            tgm_c%sym_Zu_u(i1)%c(1) = i1
         enddo
      else
         tgm_c%n_nz_Zu_u  = 0_int32
         tgm_c%n_sym_Zu_u = 0_int32
      endif
      if ( tgm_c%d_Zd .gt. n_forbs  ) then
         tgm_c%n_nz_Zd_u = n_forbs
         allocate(tgm_c%nz_Zd_u(1:tgm_c%n_nz_Zd_u))
         do i1 = 1_int32, n_forbs
            tgm_c%nz_Zd_u(i1)%o1 = i1
            if ( tgm_c%d_Zd.gt.n_forbs) then
               tgm_c%nz_Zd_u(i1)%o2 = tgm_c%d_Zd
               tgm_c%nz_Zd_u(i1)%io = i1 + n_forbs**2
            else
               tgm_c%nz_Zd_u(i1)%o2 = 1_int32
               tgm_c%nz_Zd_u(i1)%io = i1
            endif
         enddo
         tgm_c%n_sym_Zd_u = n_forbs
         allocate(tgm_c%sym_Zd_u(1:tgm_c%n_sym_Zd_u))
         do i1 = 1_int32, tgm_c%n_sym_Zd_u
            tgm_c%sym_Zd_u(i1)%is_indp = .true.
            tgm_c%sym_Zd_u(i1)%n_c = 1_int32
            allocate( tgm_c%sym_Zd_u(i1)%c(1) )
            tgm_c%sym_Zd_u(i1)%c(1) = i1
         enddo
      else
         tgm_c%n_nz_Zd_u  = 0_int32
         tgm_c%n_sym_Zd_u = 0_int32
      endif
   end subroutine cmp_tgmcfu_st
   subroutine read_tgmcfu_st( tgm_c, head_fle )
      type(tgmcfs_t), target, intent(inout) :: tgm_c
      character(4),           intent(in)    :: head_fle
      integer(int32) :: o1
      integer(int32) :: i1, i2
      integer(int32) :: n_c
      integer(int32), allocatable, dimension(:) :: c
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '_s.sav'
      tgm_c%n_nz_Z  = 0_int32
      tgm_c%n_sym_Z = 0_int32
      tgm_c%n_nz_Zu_u  = 0_int32
      tgm_c%n_nz_Zd_u  = 0_int32
      tgm_c%n_sym_Zu_u = 0_int32
      tgm_c%n_sym_Zd_u = 0_int32
      if(mpi_rank.eq.0_int32) then
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*) !'("# Table of non zero elements of triplet matrix")'
         read(1,*) tgm_c%n_nz_Z
         allocate( tgm_c%nz_Z(1:tgm_c%n_nz_Z) )
         do i1 = 1_int32, tgm_c%n_nz_Z
            read(1,*) tgm_c%nz_Z(i1)%o1, tgm_c%nz_Z(i1)%o2
            tgm_c%nz_Z(i1)%io = tgm_c%nz_Z(i1)%o1 + n_forbs * ( tgm_c%nz_Z(i1)%o2 -1 )
         enddo
         allocate( c(1:2*tgm_c%n_nz_Z) )
         read(1,*)  !'("# Symmetry table of triplet matrix")'
         read(1,*) tgm_c%n_sym_Z
         allocate(tgm_c%sym_Z(1:tgm_c%n_sym_Z))
         do i1 = 1_int32, tgm_c%n_sym_Z
            read(1,*) n_c, ( c(2*i2-1:2*i2), i2 = 1_int32, n_c )
            tgm_c%sym_Z(i1)%is_indp = .true.
            tgm_c%sym_Z(i1)%n_c     = n_c
            allocate( tgm_c%sym_Z(i1)%c(1:n_c) )
            do i2 = 1_int32, n_c
               o1 = abs(c(2*i2-1)) + n_forbs * ( c(2*i2) - 1 )
               call conv_indx( tgm_c%n_nz_Z, tgm_c%nz_Z, o1 )
               tgm_c%sym_Z(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
            enddo 
         enddo
         deallocate(c)
         if ( tgm_c%d_Zu.ne.n_forbs ) then
            read(1,*) !'("# Table of non zero elements spin up unpaired orbital")')
            read(1,*) tgm_c%n_nz_Zu_u
            allocate( tgm_c%nz_Zu_u(1:tgm_c%n_nz_Zu_u) )
            do i1 = 1_int32, tgm_c%n_nz_Zu_u
               read(1,*) tgm_c%nz_Zu_u(i1)%o1, tgm_c%nz_Zu_u(i1)%o2
               tgm_c%nz_Zu_u(i1)%io = tgm_c%nz_Zu_u(i1)%o1 + n_forbs * ( tgm_c%nz_Zu_u(i1)%o2 - 1 )
            enddo
            allocate( c(1:2*tgm_c%n_nz_Zu_u) ) ; c = 0_int32
            read(1,*) !'("# Symmetry table of spin up unpaired orbital")')
            read(1,*) tgm_c%n_sym_Zu_u
            allocate(tgm_c%sym_Zu_u(1:tgm_c%n_sym_Zu_u))
            do i1 = 1_int32, tgm_c%n_sym_Zu_u
               read(1,*) n_c, c(1:2*n_c)
               tgm_c%sym_Zu_u(i1)%is_indp = .true.
               tgm_c%sym_Zu_u(i1)%n_c     = n_c
               allocate( tgm_c%sym_Zu_u(i1)%c(1:n_c) )
               do i2 = 1_int32, n_c
                  o1 = abs(c(2*i2-1)) + n_forbs * ( c(2*i2) - 1 )
                  call conv_indx( tgm_c%n_nz_Zu_u, tgm_c%nz_Zu_u, o1 )
                  tgm_c%sym_Zu_u(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
               enddo 
            enddo
            deallocate( c )
         endif
         if( tgm_c%d_Zd.eq.1_int32 .or. tgm_c%d_Zd.gt.n_forbs  ) then
            read(1,*) !'("# Table of non zero elements spin up unpaired orbital")')
            read(1,*) tgm_c%n_nz_Zd_u
            allocate( tgm_c%nz_Zd_u(1:tgm_c%n_nz_Zd_u) )
            do i1 = 1_int32, tgm_c%n_nz_Zd_u
               read(1,*) tgm_c%nz_Zd_u(i1)%o1, tgm_c%nz_Zd_u(i1)%o2
               tgm_c%nz_Zd_u(i1)%io = tgm_c%nz_Zd_u(i1)%o1 + n_forbs * (tgm_c%nz_Zd_u(i1)%o2 -1 )
            enddo
            allocate( c(1:2*tgm_c%n_nz_Zd_u) )
            read(1,*) !'("# Symmetry table of spin up unpaired orbital")')
            read(1,*) tgm_c%n_sym_Zd_u
            allocate(tgm_c%sym_Zd_u(1:tgm_c%n_sym_Zd_u))
            do i1 = 1_int32, tgm_c%n_sym_Zd_u
               read(1,'(1000I6)') n_c, ( c(2*i2-1:2*i2), i2 = 1_int32, n_c )
               tgm_c%sym_Zd_u(i1)%is_indp = .true.
               tgm_c%sym_Zd_u(i1)%n_c     = n_c
               allocate( tgm_c%sym_Zd_u(i1)%c(1:n_c) )
               do i2 = 1_int32, n_c
                  o1 = abs(c(2*i2-1)) + n_forbs * ( c(2*i2) - 1 )
                  call conv_indx( tgm_c%n_nz_Zd_u, tgm_c%nz_Zd_u, o1 )
                  tgm_c%sym_Zd_u(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
               enddo 
            enddo
            deallocate( c )
         endif
         close(1)
      endif 
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(tgm_c%n_nz_Z, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(tgm_c%n_sym_Z,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      if ( mpi_rank.ne.0_int32 ) allocate( tgm_c%nz_Z(1:tgm_c%n_nz_Z) )
      do i1 = 1_int32, tgm_c%n_nz_Z
         call mpi_bcast(tgm_c%nz_Z(i1)%o1, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(tgm_c%nz_Z(i1)%o2, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(tgm_c%nz_Z(i1)%io, 1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      if ( mpi_rank.ne.0_int32 ) allocate( tgm_c%sym_Z(1:tgm_c%n_sym_Z) )
      do i1 = 1_int32, tgm_c%n_sym_Z
         call mpi_bcast(tgm_c%sym_Z(i1)%n_c,     1, MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(tgm_c%sym_Z(i1)%is_indp, 1, MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
         if ( mpi_rank.ne.0_int32 ) allocate( tgm_c%sym_Z(i1)%c(1:tgm_c%sym_Z(i1)%n_c) )
         call mpi_bcast(tgm_c%sym_Z(i1)%c(:),tgm_c%sym_Z(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      call mpi_bcast(tgm_c%n_nz_Zd_u, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(tgm_c%n_nz_Zu_u, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(tgm_c%n_sym_Zu_u,1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(tgm_c%n_sym_Zd_u,1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      if ( tgm_c%n_nz_Zu_u.gt.0_int32 ) then
         if ( mpi_rank.ne.0_int32 ) allocate( tgm_c%nz_Zu_u(1:tgm_c%n_nz_Zu_u) )
         do i1 = 1_int32, tgm_c%n_nz_Zu_u
            call mpi_bcast(tgm_c%nz_Zu_u(i1)%o1, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(tgm_c%nz_Zu_u(i1)%o2, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(tgm_c%nz_Zu_u(i1)%io, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         enddo
         if ( mpi_rank.ne.0_int32 ) allocate( tgm_c%sym_Zu_u(1:tgm_c%n_sym_Zu_u) )
         do i1 = 1_int32, tgm_c%n_sym_Zu_u
            call mpi_bcast(tgm_c%sym_Zu_u(i1)%n_c,     1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(tgm_c%sym_Zu_u(i1)%is_indp, 1_int32,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
            if ( mpi_rank.ne.0_int32 ) allocate( tgm_c%sym_Zu_u(i1)%c(1:tgm_c%sym_Zu_u(i1)%n_c) )
            call mpi_bcast(tgm_c%sym_Zu_u(i1)%c(:), tgm_c%sym_Zu_u(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         enddo
      endif
      if ( tgm_c%n_nz_Zd_u.gt.0_int32 ) then
         if ( mpi_rank.ne.0_int32 ) allocate( tgm_c%nz_Zd_u(1:tgm_c%n_nz_Zd_u) )
         do i1 = 1_int32, tgm_c%n_nz_Zd_u
            call mpi_bcast(tgm_c%nz_Zd_u(i1)%o1, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(tgm_c%nz_Zd_u(i1)%o2, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(tgm_c%nz_Zd_u(i1)%io, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         enddo
         if ( mpi_rank.ne.0_int32 ) allocate( tgm_c%sym_Zd_u(1:tgm_c%n_sym_Zd_u) )
         do i1 = 1_int32, tgm_c%n_sym_Zd_u
            call mpi_bcast(tgm_c%sym_Zd_u(i1)%n_c,     1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
            call mpi_bcast(tgm_c%sym_Zd_u(i1)%is_indp, 1_int32,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
            if ( mpi_rank.ne.0_int32 ) allocate( tgm_c%sym_Zd_u(i1)%c(1:tgm_c%sym_Zd_u(i1)%n_c) )
            call mpi_bcast(tgm_c%sym_Zd_u(i1)%c(:), tgm_c%sym_Zd_u(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         enddo
      endif
#endif
   end subroutine read_tgmcfu_st
   subroutine read_tgmcfs( tgm_c, n_fe_d, head_fle )
      type(tgmcfs_t), intent(inout) :: tgm_c
      integer(int32), intent(in)    :: n_fe_d
      character(4),   intent(in)    :: head_fle
      integer(int32) :: o1, o2
      integer(int32) :: i1, n_i
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
      if(mpi_rank.eq.0_int32) then
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*) !'("# Triplet geminal Matrix coefficients Spin up")')
         read(1,*) n_i
         do i1 = 1_int32, n_i
            read(1,*) o1, o2, tgm_c%mat_Zu(o1,o2)
            tgm_c%mat_Zu(o2,o1) = - tgm_c%mat_Zu(o1,o2)
         enddo
         if(tgm_c%d_Zu.gt.n_forbs) then
            read(1,*) !'("# Coefficients of the Spin up unpaired orbital")')
            read(1,*) n_i
            do i1 = 1_int32, n_i
               read(1,*) o1, o2, tgm_c%mat_Zu(o1,o2)
            enddo
         endif
         if( n_fe_d.gt.1_int32 ) then
            read(1,*) !'("# Triplet geminal Matrix coefficients Spin down")')
            read(1,*) n_i
            do i1 = 1_int32, n_i
               read(1,*) o1, o2, tgm_c%mat_Zd(o1,o2)
               tgm_c%mat_Zd(o2,o1) = - tgm_c%mat_Zd(o1,o2)
            enddo
         endif
         if( tgm_c%d_Zd.gt.n_forbs ) then
            read(1,*) !'("# Coefficients of the Spin do unpaired orbital")')
            read(1,*) n_i
            do i1 = 1_int32, n_i
               read(1,*) o1, o2, tgm_c%mat_Zd(o1,o2)
            enddo
         endif
         close(1)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast( tgm_c%mat_Zu, n_forbs * tgm_c%d_Zu, MPI_DOUBLE_PRECISION, &
      & 0, MPI_COMM_WORLD, mpierr)
      if ( n_fe_d.gt.0_int32 ) then
         call mpi_bcast( tgm_c%mat_Zd, n_forbs * tgm_c%d_Zd, MPI_DOUBLE_PRECISION, &
         & 0, MPI_COMM_WORLD, mpierr)
      endif
#endif
   end subroutine read_tgmcfs
   subroutine save_tgmcfs( tgm_c, head_fle )
      type(tgmcfs_t), target, intent(inout) :: tgm_c
      character(4),           intent(in)    :: head_fle
      integer(int32), pointer :: o1, o2
      integer(int32) :: i1
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# Triplet geminal Matrix coefficients Spin up")')
      write(1,'(I6)') tgm_c%n_nz_Z / 2_int32
      do i1 = 1_int32, tgm_c%n_nz_Z
         o1 => tgm_c%nz_Z(i1)%o1
         o2 => tgm_c%nz_Z(i1)%o2
         if (o1.lt.o2) then
            write(1,'(2I6,2E18.9)') o1, o2, tgm_c%mat_Zu(o1,o2)
         endif
      enddo
      if( tgm_c%n_nz_Zu_u.gt.0_int32 ) then
         write(1,'("# Coefficients of the Spin up unpaired orbital")')
         write(1,'(I6)') tgm_c%n_nz_Zu_u
         do i1 = 1_int32, tgm_c%n_nz_Zu_u
            o1 => tgm_c%nz_Zu_u(i1)%o1
            o2 => tgm_c%nz_Zu_u(i1)%o2
            write(1,'(2I6,2E18.9)') o1, o2, tgm_c%mat_Zu(o1,o2)
         enddo
      endif
      if( tgm_c%d_Zd.ge.n_forbs ) then
         write(1,'("# Triplet geminal Matrix coefficients Spin down")')
         write(1,'(I6)') tgm_c%n_nz_Z / 2_int32
         do i1 = 1_int32, tgm_c%n_nz_Z
            o1 => tgm_c%nz_Z(i1)%o1
            o2 => tgm_c%nz_Z(i1)%o2
            if (o1.lt.o2) then
               write(1,'(2I6,2E18.9)') o1, o2, tgm_c%mat_Zd(o1,o2)
            endif
         enddo
      endif
      if( tgm_c%n_nz_Zd_u.gt.0_int32 ) then
         write(1,'("# Coefficients of the Spin do unpaired orbital")')
         write(1,'(I6)') tgm_c%n_nz_Zd_u
         do i1 = 1_int32, tgm_c%n_nz_Zd_u
            o1 => tgm_c%nz_Zd_u(i1)%o1
            o2 => tgm_c%nz_Zd_u(i1)%o2
            write(1,'(2I6,2E18.9)') o1, o2, tgm_c%mat_Zd(o1,o2)
         enddo
      endif
      close(1)
   end subroutine save_tgmcfs
   subroutine save_tgmcfs_st(tgm_c, head_fle)
      type(tgmcfs_t), target, intent(inout) :: tgm_c
      character(4),           intent(in)    :: head_fle
      integer(int32), pointer :: o1
      integer(int32)          :: i1, i2
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '_s.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# Table of non zero elements of Zeta matrices")')
      write(1,'(I6)') tgm_c%n_nz_Z
      do i1 = 1_int32, tgm_c%n_nz_Z
         write(1,'(2I6)') tgm_c%nz_Z(i1)%o1, tgm_c%nz_Z(i1)%o2
      enddo
      write(1,'("# Symmetry table of the Zeta matrices")')
      write(1,'(I6)') tgm_c%n_sym_Z
      do i1 = 1_int32, tgm_c%n_sym_Z
         o1 => tgm_c%sym_Z(i1)%n_c
         write(1,'(1000I6)') o1,( sign(1,tgm_c%sym_Z(i1)%c(i2)) * tgm_c%nz_Z(abs(tgm_c%sym_Z(i1)%c(i2)))%o1, &
         & tgm_c%nz_Z(abs(tgm_c%sym_Z(i1)%c(i2)))%o2, i2 = 1_int32, o1 )
      enddo
      if ( tgm_c%n_sym_Zu_u.gt.0_int32 ) then
         write(1,'("# Table of non zero elements spin up unpaired orbital")')
         write(1,'(I6)') tgm_c%n_nz_Zu_u
         do i1 = 1_int32, tgm_c%n_nz_Zu_u
            write(1,'(2I6)') tgm_c%nz_Zu_u(i1)%o1, tgm_c%nz_Zu_u(i1)%o2
         enddo
         write(1,'("# Symmetry table of spin up unpaired orbital")')
         write(1,'(I6)') tgm_c%n_sym_Zu_u
         do i1 = 1_int32, tgm_c%n_sym_Zu_u
            o1 => tgm_c%sym_Zu_u(i1)%n_c
            write(1,'(1000I6)') o1,( sign(1,tgm_c%sym_Zu_u(i1)%c(i2)) * tgm_c%nz_Zu_u(abs(tgm_c%sym_Zu_u(i1)%c(i2)))%o1, &
            & tgm_c%nz_Zu_u(abs(tgm_c%sym_Zu_u(i1)%c(i2)))%o2, i2 = 1_int32, o1)
         enddo
      endif
      if( tgm_c%n_sym_Zd_u.gt.0_int32 ) then
         write(1,'("# Table of non zero elements spin up unpaired orbital")')
         write(1,'(I6)') tgm_c%n_nz_Zd_u
         do i1 = 1_int32, tgm_c%n_nz_Zd_u
            write(1,'(2I6)') tgm_c%nz_Zd_u(i1)%o1, tgm_c%nz_Zd_u(i1)%o2
         enddo
         write(1,'("# Symmetry table of spin do unpaired orbital")')
         write(1,'(I6)') tgm_c%n_sym_Zd_u
         do i1 = 1_int32, tgm_c%n_sym_Zd_u
            o1 => tgm_c%sym_Zd_u(i1)%n_c
            write(1,'(1000I6)') o1,( sign(1,tgm_c%sym_Zd_u(i1)%c(i2)) * tgm_c%nz_Zd_u(abs(tgm_c%sym_Zd_u(i1)%c(i2)))%o1, &
            & tgm_c%nz_Zd_u(abs(tgm_c%sym_Zd_u(i1)%c(i2)))%o2, i2 = 1_int32, o1)
         enddo
      endif
      close(1)
   end subroutine save_tgmcfs_st
   subroutine upd_tgmcfs( tgm_c, n_par_det_s, vec_tgm_var )
      type(tgmcfs_t), target,           intent(inout) :: tgm_c
      integer(int32),                   intent(in)    :: n_par_det_s
      real(dp), dimension(n_par_det_s), intent(in)    :: vec_tgm_var
      integer(int32), pointer :: o1, o2
      integer(int32)          :: i1, i2
      integer(int32)          :: ip
      ip = 1_int32
      if ( spin.eq.'R' ) then
         do i1 = 1_int32, tgm_c%n_sym_Z
            do i2 = 1_int32, tgm_c%sym_Z(i1)%n_c
               o1 => tgm_c%nz_Z(abs(tgm_c%sym_Z(i1)%c(i2)))%o1
               o2 => tgm_c%nz_Z(abs(tgm_c%sym_Z(i1)%c(i2)))%o2
               tgm_c%mat_Zu(o1,o2) = tgm_c%mat_Zu(o1,o2) + dble(sign(1,tgm_c%sym_Z(i1)%c(i2))) * vec_tgm_var(ip)
               if( tgm_c%d_Zd.ge.n_forbs ) tgm_c%mat_Zd(o1,o2) = tgm_c%mat_Zu(o1,o2)
            enddo
            ip = ip + 1_int32
         enddo
         if ( tgm_c%n_sym_Zu_u.gt.0_int32 ) then
            do i1 = 1_int32, tgm_c%n_sym_Zu_u
               do i2 = 1_int32, tgm_c%sym_Zu_u(i1)%n_c
                  o1 => tgm_c%nz_Zu_u(abs(tgm_c%sym_Zu_u(i1)%c(i2)))%o1
                  o2 => tgm_c%nz_Zu_u(abs(tgm_c%sym_Zu_u(i1)%c(i2)))%o2
                  tgm_c%mat_Zu(o1,o2) = tgm_c%mat_Zu(o1,o2) + dble(sign(1,tgm_c%sym_Zu_u(i1)%c(i2))) * vec_tgm_var(ip)
                  if ( tgm_c%n_sym_Zd_u.gt.0_int32  ) tgm_c%mat_Zd(o1,tgm_c%d_Zd) = tgm_c%mat_Zu(o1,o2)
               enddo
               ip = ip + 1_int32
            enddo
         elseif ( tgm_c%n_sym_Zd_u.gt.0_int32 ) then
            do i1 = 1_int32, tgm_c%n_sym_Zd_u
               do i2 = 1_int32, tgm_c%sym_Zd_u(i1)%n_c
                  o1 => tgm_c%nz_Zd_u(abs(tgm_c%sym_Zd_u(i1)%c(i2)))%o1
                  o2 => tgm_c%nz_Zd_u(abs(tgm_c%sym_Zd_u(i1)%c(i2)))%o2
                  tgm_c%mat_Zd(o1,o2) = tgm_c%mat_Zd(o1,o2) + dble(sign(1,tgm_c%sym_Zd_u(i1)%c(i2))) * vec_tgm_var(ip)
               enddo
               ip = ip + 1_int32
            enddo
         endif
      else 
         do i1 = 1_int32, tgm_c%n_sym_Z
            do i2 = 1_int32, tgm_c%sym_Z(i1)%n_c
               o1 => tgm_c%nz_Z(abs(tgm_c%sym_Z(i1)%c(i2)))%o1
               o2 => tgm_c%nz_Z(abs(tgm_c%sym_Z(i1)%c(i2)))%o2
               tgm_c%mat_Zu(o1,o2) = tgm_c%mat_Zu(o1,o2) + dble(sign(1,tgm_c%sym_Z(i1)%c(i2))) * vec_tgm_var(ip)
            enddo
            ip = ip + 1_int32
         enddo
         if ( tgm_c%n_sym_Zu_u.gt.0_int32 ) then
            do i1 = 1_int32, tgm_c%n_sym_Zu_u
               do i2 = 1_int32, tgm_c%sym_Zu_u(i1)%n_c
                  o1 => tgm_c%nz_Zu_u(abs(tgm_c%sym_Zu_u(i1)%c(i2)))%o1
                  o2 => tgm_c%nz_Zu_u(abs(tgm_c%sym_Zu_u(i1)%c(i2)))%o2
                  tgm_c%mat_Zu(o1,o2) = tgm_c%mat_Zu(o1,o2) + dble(sign(1,tgm_c%sym_Zu_u(i1)%c(i2))) * vec_tgm_var(ip)
               enddo
               ip = ip + 1_int32
            enddo
         endif
         if ( tgm_c%d_Zd .ge. n_forbs ) then
            do i1 = 1_int32, tgm_c%n_sym_Z
               do i2 = 1_int32, tgm_c%sym_Z(i1)%n_c
                  o1 => tgm_c%nz_Z(abs(tgm_c%sym_Z(i1)%c(i2)))%o1
                  o2 => tgm_c%nz_Z(abs(tgm_c%sym_Z(i1)%c(i2)))%o2
                  tgm_c%mat_Zd(o1,o2) = tgm_c%mat_Zd(o1,o2) + dble(sign(1,tgm_c%sym_Z(i1)%c(i2))) * vec_tgm_var(ip)
               enddo
               ip = ip + 1_int32
            enddo
         endif
         if ( tgm_c%n_sym_Zd_u.gt.0_int32 ) then
            do i1 = 1_int32, tgm_c%n_sym_Zd_u
               do i2 = 1_int32, tgm_c%sym_Zd_u(i1)%n_c
                  o1 => tgm_c%nz_Zd_u(abs(tgm_c%sym_Zd_u(i1)%c(i2)))%o1
                  o2 => tgm_c%nz_Zd_u(abs(tgm_c%sym_Zd_u(i1)%c(i2)))%o2
                  tgm_c%mat_Zd(o1,o2) = tgm_c%mat_Zd(o1,o2) + dble(sign(1,tgm_c%sym_Zd_u(i1)%c(i2))) * vec_tgm_var(ip)
               enddo
               ip = ip + 1_int32
            enddo
         endif
      endif
   end subroutine upd_tgmcfs
   subroutine sym_tgmcfs( tgm_c, n_par_det, n_par_det_s, dp_ln_pff_T, dp_ln_pff_T_s )
      type(tgmcfs_t), target,           intent(in)    :: tgm_c
      integer(int32),                   intent(in)    :: n_par_det
      integer(int32),                   intent(in)    :: n_par_det_s
      real(dp), dimension(n_par_det),   intent(in)    :: dp_ln_pff_T
      real(dp), dimension(n_par_det_s), intent(inout) :: dp_ln_pff_T_s
      integer(int32) :: i1, i2, i3, i4
      i3 = 1_int32
      i4 = 1_int32
      if (spin.eq.'R') then
         do i1 = 1_int32, tgm_c%n_sym_Z
            dp_ln_pff_T_s(i3) = dp_ln_pff_T(abs(tgm_c%sym_Z(i1)%c(1)))
            if (tgm_c%sym_Z(i1)%n_c.gt.1_int32) then
               do i2 = 2_int32, tgm_c%sym_Z(i1)%n_c
                  dp_ln_pff_T_s(i3) = dp_ln_pff_T_s(i3) &
                  & + dble(sign(1,tgm_c%sym_Z(i1)%c(i2))) * dp_ln_pff_T(abs(tgm_c%sym_Z(i1)%c(i2)))
               enddo
            endif
            i3 = i3 + 1_int32
         enddo
         i4 = i4 + tgm_c%n_nz_Z
         if( tgm_c%n_sym_Zu_u.gt.0_int32 ) then
            do i1 = 1_int32, tgm_c%n_sym_Zu_u
               dp_ln_pff_T_s(i3) = dp_ln_pff_T(i4-1+abs(tgm_c%sym_Zu_u(i1)%c(1)))
               if (tgm_c%sym_Zu_u(i1)%n_c.gt.1_int32) then
                  do i2 = 2_int32, tgm_c%sym_Zu_u(i1)%n_c
                     dp_ln_pff_T_s(i3) = dp_ln_pff_T_s(i3) &
                     & + dble(sign(1,tgm_c%sym_Zu_u(i1)%c(i2))) * dp_ln_pff_T(i4-1+abs(tgm_c%sym_Zu_u(i1)%c(i2)))
                  enddo
               endif
               i3 = i3 + 1_int32
            enddo
         elseif ( tgm_c%n_sym_Zd_u.gt.0_int32 ) then
            do i1 = 1_int32, tgm_c%n_sym_Zd_u
               dp_ln_pff_T_s(i3) = dp_ln_pff_T(i4-1+abs(tgm_c%sym_Zd_u(i1)%c(1)))
               if (tgm_c%sym_Zd_u(i1)%n_c.gt.1_int32) then
                  do i2 = 2_int32, tgm_c%sym_Zd_u(i1)%n_c
                     dp_ln_pff_T_s(i3) = dp_ln_pff_T_s(i3) &
                     & + dble(sign(1,tgm_c%sym_Zd_u(i1)%c(i2))) * dp_ln_pff_T(i4-1+abs(tgm_c%sym_Zd_u(i1)%c(i2)))
                  enddo
               endif
               i3 = i3 + 1_int32
            enddo
         endif
      else ! if (spin.eq.'U')
         do i1 = 1_int32, tgm_c%n_sym_Z
            dp_ln_pff_T_s(i3) = dp_ln_pff_T(abs(tgm_c%sym_Z(i1)%c(1)))
            if (tgm_c%sym_Z(i1)%n_c.gt.1_int32) then
               do i2 = 2_int32, tgm_c%sym_Z(i1)%n_c
                  dp_ln_pff_T_s(i3) = dp_ln_pff_T_s(i3) &
                  & + dble(sign(1,tgm_c%sym_Z(i1)%c(i2))) * dp_ln_pff_T(abs(tgm_c%sym_Z(i1)%c(i2)))
               enddo
            endif
            i3 = i3 + 1_int32
         enddo
         i4 = i4 + tgm_c%n_nz_Z
         if ( tgm_c%n_sym_Zu_u.gt.0_int32 ) then
            do i1 = 1_int32, tgm_c%n_sym_Zu_u
               dp_ln_pff_T_s(i3) = dp_ln_pff_T(i4-1+abs(tgm_c%sym_Zu_u(i1)%c(1)))
               if (tgm_c%sym_Zu_u(i1)%n_c.gt.1_int32) then
                  do i2 = 2_int32, tgm_c%sym_Zu_u(i1)%n_c
                     dp_ln_pff_T_s(i3) = dp_ln_pff_T_s(i3) &
                     & + dble(sign(1,tgm_c%sym_Zu_u(i1)%c(i2))) * dp_ln_pff_T(i4-1+abs(tgm_c%sym_Zu_u(i1)%c(i2)))
                  enddo
               endif
               i3 = i3 + 1_int32
            enddo
            i4 = i4 + tgm_c%n_nz_Zu_u
         endif
         if( tgm_c%d_Zd.ge.n_forbs ) then
            do i1 = 1_int32, tgm_c%n_sym_Z
               dp_ln_pff_T_s(i3) = dp_ln_pff_T(i4-1+abs(tgm_c%sym_Z(i1)%c(1)))
               if (tgm_c%sym_Z(i1)%n_c.gt.1_int32) then
                  do i2 = 2_int32, tgm_c%sym_Z(i1)%n_c
                     dp_ln_pff_T_s(i3) = dp_ln_pff_T_s(i3) &
                     & + dble(sign(1,tgm_c%sym_Z(i1)%c(i2))) * dp_ln_pff_T(i4-1+abs(tgm_c%sym_Z(i1)%c(i2)))
                  enddo
               endif
               i3 = i3 + 1_int32
            enddo
            i4 = i4 + tgm_c%n_nz_Z
         endif
         if( tgm_c%n_sym_Zd_u.gt.0_int32 ) then
            do i1 = 1_int32, tgm_c%n_sym_Zd_u
               dp_ln_pff_T_s(i3) = dp_ln_pff_T(i4-1+abs(tgm_c%sym_Zd_u(i1)%c(1)))
               if (tgm_c%sym_Zd_u(i1)%n_c.gt.1_int32) then
                  do i2 = 2_int32, tgm_c%sym_Zd_u(i1)%n_c
                     dp_ln_pff_T_s(i3) = dp_ln_pff_T_s(i3) &
                     & + dble(sign(1,tgm_c%sym_Zd_u(i1)%c(i2))) * dp_ln_pff_T(i4-1+abs(tgm_c%sym_Zd_u(i1)%c(i2)))
                  enddo
               endif
               i3 = i3 + 1_int32
            enddo
         endif
      endif 
   end subroutine sym_tgmcfs
end module triplet_geminal_params_m
