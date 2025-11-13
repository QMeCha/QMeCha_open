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
module slater_determinant_params_m
   use fortran_kinds_v,          only: dp, int32
   use fermionic_wavefunction_v, only: spin
   use write_lines_m
   use geminal_symmetries_m,     only: lmbd_ele_t, lmbd_sym_t, conv_indx
   use openmp_mpi_m
   use fermionic_orbitals_m,     only: n_forbs
   use quantum_monte_carlo_v,    only: restart
   implicit none
   type, public :: sldcfs_t
      real(dp),         allocatable, dimension(:,:) :: mat_L_u, mat_L_d
      integer(int32)                                :: n_nz_Lu, n_nz_Ld
      type(Lmbd_ele_t), allocatable, dimension(:)   :: nz_Lu, nz_Ld
      integer(int32)                                :: n_sym_Lu, n_sym_Ld
      type(Lmbd_sym_t), allocatable, dimension(:)   :: sym_Lu, sym_Ld
   end type sldcfs_t
   type(sldcfs_t), public, save, target :: sld_c_e
   type(sldcfs_t), public, save, target :: sld_c_p
   public  :: ini_sldcfs, save_sldcfs_st, upd_sldcfs, save_sldcfs, sym_sldcfs
   private :: cmp_sldcfs_st, read_sldcfs_st, read_sldcfs
contains
   subroutine ini_sldcfs( sldcfs_in, n_fe_u, n_fe_d, dl_opt, n_par_det, n_par_det_s, head_fle )
      type(sldcfs_t),         intent(inout) :: sldcfs_in
      integer(int32),         intent(in)    :: n_fe_u, n_fe_d
      logical,                intent(in)    :: dl_opt
      integer(int32),         intent(inout) :: n_par_det
      integer(int32),         intent(inout) :: n_par_det_s
      character(4), optional, intent(in)    :: head_fle
      integer(int32)  :: i1
      allocate( sldcfs_in%mat_L_u(1:n_forbs, 1:n_fe_u) ) ; sldcfs_in%mat_L_u = 0.0_dp
      do i1 = 1_int32, n_fe_u
         sldcfs_in%mat_L_u(i1,i1) = 1.0_dp
      enddo 
      if (n_fe_d.gt.0_int32) then
         allocate( sldcfs_in%mat_L_d(1:n_forbs, 1:n_fe_d) ) ; sldcfs_in%mat_L_d = sldcfs_in%mat_L_u(:, 1:n_fe_d)
      endif ! n_fe_d.gt.0_int32
      if ( restart ) then
         call read_sldcfs( sldcfs_in, n_fe_u, n_fe_d, head_fle )
         call read_sldcfs_st( sldcfs_in, n_fe_d, head_fle )
      else
         call cmp_sldcfs_st( sldcfs_in, n_fe_u, n_fe_d )
      endif
      if ( dl_opt ) then
         if (spin.eq.'R') then
            n_par_det   = sldcfs_in%n_nz_Lu
            n_par_det_s = sldcfs_in%n_sym_Lu
         else
            n_par_det   = sldcfs_in%n_nz_Lu  + sldcfs_in%n_nz_Ld
            n_par_det_s = sldcfs_in%n_sym_Lu + sldcfs_in%n_sym_Ld
         endif
      else
         n_par_det   = 0_int32
         n_par_det_s = 0_int32
      endif
   end subroutine ini_sldcfs
   subroutine cmp_sldcfs_st( sldcfs_in, n_fe_u, n_fe_d )
      type(sldcfs_t), intent(inout) :: sldcfs_in
      integer(int32) , intent(in) :: n_fe_u, n_fe_d
      integer(int32) :: i1, i2, i3
      sldcfs_in%n_nz_Lu  = n_forbs * n_fe_u
      sldcfs_in%n_nz_Ld  = n_forbs * n_fe_d
      sldcfs_in%n_sym_Lu = n_forbs * n_fe_u
      sldcfs_in%n_sym_Ld = n_forbs * n_fe_d
      allocate( sldcfs_in%nz_Lu(1:sldcfs_in%n_nz_Lu), sldcfs_in%nz_Ld(1:sldcfs_in%n_nz_Ld) )
      i3 = 1_int32
      do i2 = 1_int32, n_fe_u ; do i1 = 1_int32, n_forbs
            sldcfs_in%nz_Lu(i3)%o1 = i1
            sldcfs_in%nz_Lu(i3)%o2 = i2
            sldcfs_in%nz_Lu(i3)%io = i1 + n_forbs * (i2 - 1)
            if (i2.le.n_fe_d) then
               sldcfs_in%nz_Ld(i3)%o1 = i1
               sldcfs_in%nz_Ld(i3)%o2 = i2
               sldcfs_in%nz_Ld(i3)%io = i1 + n_forbs * (i2 - 1)
            endif
            i3 = i3 + 1_int32
         enddo ; enddo
      allocate( sldcfs_in%sym_Lu(1:sldcfs_in%n_sym_Lu), sldcfs_in%sym_Ld(1:sldcfs_in%n_sym_Ld) )
      do i1 = 1_int32, sldcfs_in%n_sym_Lu
         sldcfs_in%sym_Lu(i1)%is_indp = .true.
         sldcfs_in%sym_Lu(i1)%n_c = 1_int32
         allocate( sldcfs_in%sym_Lu(i1)%c(1) ) ; sldcfs_in%sym_Lu(i1)%c(1) = i1
         if ( i1.le.sldcfs_in%n_sym_Ld) then
            sldcfs_in%sym_Ld(i1)%is_indp = .true.
            sldcfs_in%sym_Ld(i1)%n_c = 1_int32
            allocate( sldcfs_in%sym_Ld(i1)%c(1) ) ; sldcfs_in%sym_Ld(i1)%c(1) = i1
         endif
      enddo
   end subroutine cmp_sldcfs_st
   subroutine save_sldcfs_st( sldcfs_in, n_fe_d, head_fle )
      type(sldcfs_t), target , intent(inout) :: sldcfs_in
      integer(int32),         intent(in) :: n_fe_d
      character(4), optional, intent(in) :: head_fle
      character(100) :: svf_fle
      integer(int32), pointer :: o1
      integer(int32)          :: i1, i2
      svf_fle = 'wvfn.save/'//trim(head_fle)//'_s.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# Table of non zero elements spin up orbitals matrix")')
      write(1,'(I6)') sldcfs_in%n_nz_Lu
      do i1 = 1_int32, sldcfs_in%n_nz_Lu
         write(1,'(2I6)') sldcfs_in%nz_Lu(i1)%o1, sldcfs_in%nz_Lu(i1)%o2
      enddo
      write(1,'("# Symmetry table of spin up orbitals matrix")')
      write(1,'(I6)') sldcfs_in%n_sym_Lu
      do i1 = 1_int32, sldcfs_in%n_sym_Lu
         o1 => sldcfs_in%sym_Lu(i1)%n_c
         write(1,'(1000I6)') o1,( sign(1,sldcfs_in%sym_Lu(i1)%c(i2)) * sldcfs_in%nz_Lu(abs(sldcfs_in%sym_Lu(i1)%c(i2)))%o1, &
         & sldcfs_in%nz_Lu(abs(sldcfs_in%sym_Lu(i1)%c(i2)))%o2, i2 = 1_int32, o1)
      enddo
      if (n_fe_d.gt.0_int32) then
         write(1,'("# Table of non zero elements spin do orbitals matrix")')
         write(1,'(I6)') sldcfs_in%n_nz_Ld
         do i1 = 1_int32, sldcfs_in%n_nz_Ld
            write(1,'(2I6)') sldcfs_in%nz_Ld(i1)%o1, sldcfs_in%nz_Ld(i1)%o2
         enddo
         write(1,'("# Symmetry table of spin do orbitals matrix")')
         write(1,'(I6)') sldcfs_in%n_sym_Ld
         do i1 = 1_int32, sldcfs_in%n_sym_Ld
            o1 => sldcfs_in%sym_Ld(i1)%n_c
            write(1,'(1000I6)') o1, ( sign(1,sldcfs_in%sym_Ld(i1)%c(i2)) * sldcfs_in%nz_Ld(abs(sldcfs_in%sym_Ld(i1)%c(i2)))%o1, &
            & sldcfs_in%nz_Ld(abs(sldcfs_in%sym_Ld(i1)%c(i2)))%o2, i2 = 1_int32, o1)
         enddo
      endif
      close(1)
   end subroutine save_sldcfs_st
   subroutine read_sldcfs_st( sldcfs_in, n_fe_d, head_fle )
      type(sldcfs_t), intent(inout) :: sldcfs_in
      integer(int32) , intent(in) :: n_fe_d
      character(4),    intent(in) :: head_fle
      character(100) :: svf_fle
      integer(int32) :: o1
      integer(int32) :: i1, i2
      integer(int32) :: n_c
      integer(int32), allocatable, dimension(:) :: c
      svf_fle = 'wvfn.save/' // trim(head_fle) // '_s.sav'
      if(mpi_rank.eq.0_int32) then
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*) !'("# Table of non zero elements spin up unpaired orbital")')
         read(1,*) sldcfs_in%n_nz_Lu
         allocate( sldcfs_in%nz_Lu(1:sldcfs_in%n_nz_Lu) )
         do i1 = 1_int32, sldcfs_in%n_nz_Lu
            read(1,*) sldcfs_in%nz_Lu(i1)%o1, sldcfs_in%nz_Lu(i1)%o2
            sldcfs_in%nz_Lu(i1)%io = sldcfs_in%nz_Lu(i1)%o1 + n_forbs * ( sldcfs_in%nz_Lu(i1)%o2 -1 )
         enddo
         allocate( c(1:2*sldcfs_in%n_nz_Lu) )
         read(1,*) !'("# Symmetry table of spin up unpaired orbital")')
         read(1,*) sldcfs_in%n_sym_Lu
         allocate( sldcfs_in%sym_Lu(1:sldcfs_in%n_sym_Lu) )
         do i1 = 1_int32, sldcfs_in%n_sym_Lu
            read(1,*) n_c, ( c(2*i2-1:2*i2), i2 = 1_int32, n_c )
            sldcfs_in%sym_Lu(i1)%is_indp = .true.
            sldcfs_in%sym_Lu(i1)%n_c     = n_c
            allocate( sldcfs_in%sym_Lu(i1)%c(1:n_c) )
            do i2 = 1_int32, n_c
               o1 = abs(c(2*i2-1)) + n_forbs * ( c(2*i2) - 1 )
               call conv_indx( sldcfs_in%n_nz_Lu, sldcfs_in%nz_Lu, o1 )
               sldcfs_in%sym_Lu(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
            enddo ! i2 (n_c)
         enddo
         deallocate( c )
         if( n_fe_d.gt.0_int32) then
            read(1,*) !'("# Table of non zero elements spin up unpaired orbital")')
            read(1,*) sldcfs_in%n_nz_Ld
            allocate( sldcfs_in%nz_Ld(1:sldcfs_in%n_nz_Ld) )
            do i1 = 1_int32, sldcfs_in%n_nz_Ld
               read(1,*) sldcfs_in%nz_Ld(i1)%o1, sldcfs_in%nz_Ld(i1)%o2
               sldcfs_in%nz_Ld(i1)%io = sldcfs_in%nz_Ld(i1)%o1 + n_forbs * ( sldcfs_in%nz_Ld(i1)%o2 -1 )
            enddo
            allocate( c(1:2*sldcfs_in%n_nz_Ld) )
            read(1,*) !'("# Symmetry table of spin up unpaired orbital")')
            read(1,*) sldcfs_in%n_sym_Ld
            allocate( sldcfs_in%sym_Ld(1:sldcfs_in%n_sym_Ld) )
            do i1 = 1_int32, sldcfs_in%n_sym_Ld
               read(1,*) n_c, ( c(2*i2-1:2*i2), i2 = 1_int32, n_c )
               sldcfs_in%sym_Ld(i1)%is_indp = .true.
               sldcfs_in%sym_Ld(i1)%n_c     = n_c
               allocate( sldcfs_in%sym_Ld(i1)%c(1:n_c) )
               do i2 = 1_int32, n_c
                  o1 = abs(c(2*i2-1)) + n_forbs * ( c(2*i2) - 1 )
                  call conv_indx( sldcfs_in%n_nz_Ld, sldcfs_in%nz_Ld, o1 )
                  sldcfs_in%sym_Ld(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
               enddo 
            enddo
            deallocate( c )
         endif
         close(1)
      endif 
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(sldcfs_in%n_nz_Lu, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(sldcfs_in%n_sym_Lu,1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(sldcfs_in%n_nz_Ld, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(sldcfs_in%n_sym_Ld,1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      if ( mpi_rank.ne.0_int32 ) allocate( sldcfs_in%nz_Lu(1:sldcfs_in%n_nz_Lu), sldcfs_in%nz_Ld(1:sldcfs_in%n_nz_Ld))
      do i1 = 1_int32, sldcfs_in%n_nz_Lu
         call mpi_bcast(sldcfs_in%nz_Lu(i1)%o1, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(sldcfs_in%nz_Lu(i1)%o2, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(sldcfs_in%nz_Lu(i1)%io, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      do i1 = 1_int32, sldcfs_in%n_nz_Ld
         call mpi_bcast(sldcfs_in%nz_Ld(i1)%o1, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(sldcfs_in%nz_Ld(i1)%o2, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(sldcfs_in%nz_Ld(i1)%io, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      if ( mpi_rank.ne.0_int32 ) allocate( sldcfs_in%sym_Lu(1:sldcfs_in%n_sym_Lu), sldcfs_in%sym_Ld(1:sldcfs_in%n_sym_Ld) )
      do i1 = 1_int32, sldcfs_in%n_sym_Lu
         call mpi_bcast(sldcfs_in%sym_Lu(i1)%n_c,     1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(sldcfs_in%sym_Lu(i1)%is_indp, 1_int32,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
         if ( mpi_rank.ne.0_int32 ) allocate( sldcfs_in%sym_Lu(i1)%c(1:sldcfs_in%sym_Lu(i1)%n_c) )
         call mpi_bcast(sldcfs_in%sym_Lu(i1)%c(:),sldcfs_in%sym_Lu(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      do i1 = 1_int32, sldcfs_in%n_sym_Ld
         call mpi_bcast(sldcfs_in%sym_Ld(i1)%n_c,     1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(sldcfs_in%sym_Ld(i1)%is_indp, 1_int32,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
         if ( mpi_rank.ne.0_int32 ) allocate( sldcfs_in%sym_Ld(i1)%c(1:sldcfs_in%sym_Ld(i1)%n_c) )
         call mpi_bcast(sldcfs_in%sym_Ld(i1)%c(:),sldcfs_in%sym_Ld(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
#endif
   end subroutine read_sldcfs_st
   subroutine upd_sldcfs( sldcfs_in, n_fe_u, n_fe_d, n_par_det_s, vec_sld_var )
      type(sldcfs_t), target , intent(inout) :: sldcfs_in
      integer(int32) , intent(in) :: n_fe_u
      integer(int32) , intent(in) :: n_fe_d
      integer(int32) , intent(in) :: n_par_det_s
      real(dp), dimension(n_par_det_s), intent(in) :: vec_sld_var
      integer(int32), pointer :: o1, o2
      integer(int32)          :: i1, i2
      integer(int32)          :: ip
      real(dp) :: max_coef
      ip = 1_int32
      do i1 = 1_int32, sldcfs_in%n_sym_Lu
         do i2 = 1_int32, sldcfs_in%sym_Lu(i1)%n_c
            o1 => sldcfs_in%nz_Lu(abs(sldcfs_in%sym_Lu(i1)%c(i2)))%o1
            o2 => sldcfs_in%nz_Lu(abs(sldcfs_in%sym_Lu(i1)%c(i2)))%o2
            sldcfs_in%mat_L_u(o1,o2) = sldcfs_in%mat_L_u(o1,o2) + dble(sign(1,sldcfs_in%sym_Lu(i1)%c(i2))) * vec_sld_var(ip)
         enddo
         ip = ip + 1_int32
      enddo
      if (n_fe_d.gt.0_int32) then
         if (spin.eq.'R') then
            sldcfs_in%mat_L_d(:,:) = sldcfs_in%mat_L_u(:,1:n_fe_d)
         else
            do i1 = 1_int32, sldcfs_in%n_sym_Ld
               do i2 = 1_int32, sldcfs_in%sym_Ld(i1)%n_c
                  o1 => sldcfs_in%nz_Ld(abs(sldcfs_in%sym_Ld(i1)%c(i2)))%o1
                  o2 => sldcfs_in%nz_Ld(abs(sldcfs_in%sym_Ld(i1)%c(i2)))%o2
                  sldcfs_in%mat_L_d(o1,o2) = sldcfs_in%mat_L_d(o1,o2) + dble(sign(1,sldcfs_in%sym_Ld(i1)%c(i2))) * vec_sld_var(ip)
               enddo
               ip = ip + 1_int32
            enddo
         endif
      endif
      max_coef = 0.0_dp
      do i1 = 1_int32, n_fe_u
         do i2 = 1_int32, n_forbs
            if ( abs(sldcfs_in%mat_L_u(i2,i1)).gt.max_coef) max_coef = abs(sldcfs_in%mat_L_u(i2,i1))
         enddo
      enddo
      sldcfs_in%mat_L_u(:,:) = sldcfs_in%mat_L_u(:,:) / max_coef
      if (n_fe_d.gt.0_int32) sldcfs_in%mat_L_d(:,:) = sldcfs_in%mat_L_d(:,:) / max_coef
   end subroutine upd_sldcfs
   subroutine save_sldcfs( sldcfs_in, n_fe_u, n_fe_d, head_fle )
      type(sldcfs_t), intent(inout) :: sldcfs_in
      integer(int32) , intent(in) :: n_fe_u, n_fe_d
      character(4),    intent(in) :: head_fle
      integer(int32) :: i1, i2
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# Slater coefficients for Spin up matrix")')
      i1 = 1_int32
      if (n_fe_u.le.6_int32) then
         do i2 = 1_int32, n_forbs
            write(1,'(6E18.9)') sldcfs_in%mat_L_u(i2,1:n_fe_u)
         enddo
      else
         do i1 = 1_int32, int(n_fe_u/6_int32)
            do i2 = 1_int32, n_forbs
               write(1,'(6E18.9)') sldcfs_in%mat_L_u(i2,6*(i1-1)+1:6*i1)
            enddo
         enddo
         if( mod(n_fe_u, 6_int32).gt.0_int32 ) then
            do i2 = 1_int32, n_forbs
               write(1,'(6E18.9)') sldcfs_in%mat_L_u(i2,n_fe_u-mod(n_fe_u, 6_int32)+1:n_fe_u)
            enddo
         endif
      endif
      if (n_fe_d.gt.0_int32) then
         write(1,'("# Slater coefficients for Spin down matrix")')
         i1 = 1_int32
         if (n_fe_d.le.6_int32) then
            do i2 = 1_int32, n_forbs
               write(1,'(6E18.9)') sldcfs_in%mat_L_d(i2,1:n_fe_d)
            enddo
         else
            do i1 = 1_int32, int(n_fe_d/6_int32)
               do i2 = 1_int32, n_forbs
                  write(1,'(6E18.9)') sldcfs_in%mat_L_d(i2,6*(i1-1)+1:6*i1)
               enddo
            enddo
            if( mod(n_fe_d, 6_int32).gt.0_int32 ) then
               do i2 = 1_int32, n_forbs
                  write(1,'(6E18.9)') sldcfs_in%mat_L_d(i2,n_fe_d-mod(n_fe_d, 6_int32)+1:n_fe_d)
               enddo
            endif
         endif
      endif
      close(1)
   end subroutine save_sldcfs
   subroutine read_sldcfs( sldcfs_in, n_fe_u, n_fe_d, head_fle )
      type(sldcfs_t), intent(inout) :: sldcfs_in
      integer(int32) , intent(in) :: n_fe_u, n_fe_d
      character(4),    intent(in) :: head_fle
      integer(int32) :: i1, i2
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
      call write_simple_line(stdout,0,mpi_rank,2,"l","Reading Slater determinant coefficients save file.")
      if( mpi_rank.eq.0_int32 ) then
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*)
         i1 = 1_int32
         if (n_fe_u.le.6_int32) then
            do i2 = 1_int32, n_forbs
               read(1,*) sldcfs_in%mat_L_u(i2,1:n_fe_u)
            enddo
         else
            do i1 = 1_int32, int(n_fe_u/6_int32)
               do i2 = 1_int32, n_forbs
                  read(1,*) sldcfs_in%mat_L_u(i2,6*(i1-1)+1:6*i1)
               enddo
            enddo
            if( mod(n_fe_u, 6_int32).gt.0_int32 ) then
               do i2 = 1_int32, n_forbs
                  read(1,*) sldcfs_in%mat_L_u(i2,n_fe_u-mod(n_fe_u, 6_int32)+1:n_fe_u)
               enddo
            endif
         endif
         if (n_fe_d.gt.0_int32) then
            read(1,*)
            i1 = 1_int32
            if (n_fe_d.le.6_int32) then
               do i2 = 1_int32, n_forbs
                  read(1,*) sldcfs_in%mat_L_d(i2,1:n_fe_d)
               enddo
            else
               do i1 = 1_int32, int(n_fe_d/6_int32)
                  do i2 = 1_int32, n_forbs
                     read(1,*) sldcfs_in%mat_L_d(i2,6*(i1-1)+1:6*i1)
                  enddo
               enddo
               if( mod(n_fe_d, 6_int32).gt.0_int32 ) then
                  do i2 = 1_int32, n_forbs
                     read(1,*) sldcfs_in%mat_L_d(i2,n_fe_d-mod(n_fe_d, 6_int32)+1:n_fe_d)
                  enddo
               endif
            endif
         endif
         close(1)
         if (spin.eq.'R'.and.n_fe_d.gt.0_int32) then
            sldcfs_in%mat_L_d(:,1:n_fe_d) = sldcfs_in%mat_L_u(:,1:n_fe_d)
         endif
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(sldcfs_in%mat_L_u, n_forbs*n_fe_u, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
      if (n_fe_d.gt.0_int32) then
         call mpi_bcast(sldcfs_in%mat_L_d, n_forbs*n_fe_d, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
      endif ! n_fe_d.gt.0_int32
#endif
   end subroutine read_sldcfs
   subroutine sym_sldcfs( sldcfs_in, n_fe_d, n_par_det, n_par_det_s, dp_ln_det_S, dp_ln_det_S_s )
      type(sldcfs_t), intent(in) :: sldcfs_in
      integer(int32),                   intent(in)    :: n_fe_d
      integer(int32),                   intent(in)    :: n_par_det, n_par_det_s
      real(dp), dimension(n_par_det),   intent(in)    :: dp_ln_det_S
      real(dp), dimension(n_par_det_s), intent(inout) :: dp_ln_det_S_s
      integer(int32) :: i1, i2, i3, i4
      i3 = 1_int32
      i4 = 1_int32
      do i1 = 1_int32, sldcfs_in%n_sym_Lu
         dp_ln_det_S_s(i3) = dp_ln_det_S(abs(sldcfs_in%sym_Lu(i1)%c(1)))
         if ( sldcfs_in%sym_Lu(i1)%n_c.gt.1_int32 ) then
            do i2 = 2_int32, sldcfs_in%sym_Lu(i1)%n_c
               dp_ln_det_S_s(i3) = dp_ln_det_S_s(i3) + &
               & dble(sign(1,sldcfs_in%sym_Lu(i1)%c(i2))) * dp_ln_det_S(abs(sldcfs_in%sym_Lu(i1)%c(i2)))
            enddo
         endif
         i3 = i3 + 1_int32
      enddo
      i4 = i4 + sldcfs_in%n_nz_Lu
      if(n_fe_d.gt.0_int32 .and. spin.eq.'U') then
         do i1 = 1_int32, sldcfs_in%n_sym_Ld
            dp_ln_det_S_s(i3) = dp_ln_det_S(i4-1+abs(sldcfs_in%sym_Ld(i1)%c(1)))
            if ( sldcfs_in%sym_Ld(i1)%n_c.gt.1_int32 ) then
               do i2 = 2_int32, sldcfs_in%sym_Ld(i1)%n_c
                  dp_ln_det_S_s(i3) = dp_ln_det_S_s(i3) + &
                  & dble(sign(1,sldcfs_in%sym_Ld(i1)%c(i2))) * dp_ln_det_S(i4-1+abs(sldcfs_in%sym_Ld(i1)%c(i2)))
               enddo
            endif
            i3 = i3 + 1_int32
         enddo
      endif
   end subroutine sym_sldcfs
end module slater_determinant_params_m
