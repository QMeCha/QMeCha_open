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
module eppcfs_mod
   use fortran_kinds_v, only: dp, int32
   use openmp_mpi_m
   use fermionic_wavefunction_v, only: spin
   use molecular_system_v, only: n_el, n_po
   use geminal_symmetries_m, only: lmbd_ele_t, lmbd_sym_t
   use quantum_monte_carlo_v, only: restart
   use feporb_mod, only: n_eporbs
   implicit none
   type, public :: eppcfs_t
      real(dp),         allocatable, dimension(:,:) :: mat_L
      integer(int32)                                :: n_nz_L
      type(Lmbd_ele_t), allocatable, dimension(:)   :: nz_L
      integer(int32)                                :: n_sym_L
      type(Lmbd_sym_t), allocatable, dimension(:)   :: sym_L
   end type eppcfs_t
   type(eppcfs_t), public, save, target :: epo_c
   public  :: ini_eppcfs, save_eppcfs_st, upd_eppcfs, save_eppcfs, sym_eppcfs
   private :: cmp_eppcfs_st, read_eppcfs_st, read_eppcfs
contains
   subroutine ini_eppcfs( dl_opt, n_par_gem, n_par_gem_s, head_fle )
      logical,                  intent(in)    :: dl_opt
      integer(int32),           intent(inout) :: n_par_gem, n_par_gem_s
      character(4),   optional, intent(in)    :: head_fle
      allocate( epo_c%mat_L(1:n_eporbs, 1:n_el*n_po) ) ; epo_c%mat_L = 0.0_dp
      epo_c%mat_L(1,:) = 1.0_dp
      if ( restart ) then
         call read_eppcfs( head_fle )
         call read_eppcfs_st( head_fle )
      else
         call cmp_eppcfs_st(  )
      endif
      if ( dl_opt ) then
         n_par_gem   = epo_c%n_nz_L
         n_par_gem_s = epo_c%n_sym_L
      else
         n_par_gem   = 0_int32
         n_par_gem_s = 0_int32
      endif
   end subroutine ini_eppcfs
   subroutine cmp_eppcfs_st( )
      integer(int32) :: i1, i2, i3
      epo_c%n_nz_L  = n_eporbs * n_el * n_po
      epo_c%n_sym_L = n_eporbs * n_el * n_po
      allocate( epo_c%nz_L(1:epo_c%n_nz_L) )
      i3 = 1_int32
      do i2 = 1_int32, n_el*n_po ; do i1 = 1_int32, n_eporbs
            epo_c%nz_L(i3)%o1 = i1
            epo_c%nz_L(i3)%o2 = i2
            epo_c%nz_L(i3)%io = i1 + n_eporbs * (i2 - 1)
            i3 = i3 + 1_int32
         enddo ; enddo
      allocate( epo_c%sym_L(1:epo_c%n_sym_L) )
      do i1 = 1_int32, epo_c%n_sym_L
         epo_c%sym_L(i1)%is_indp = .true.
         epo_c%sym_L(i1)%n_c = 1_int32
         allocate( epo_c%sym_L(i1)%c(1) ) ; epo_c%sym_L(i1)%c(1) = i1
      enddo
   end subroutine cmp_eppcfs_st
   subroutine save_eppcfs_st( head_fle )
      character(4),  optional, intent(in)    :: head_fle
      character(100) :: svf_fle
      integer(int32), pointer :: p1
      integer(int32)          :: i1, i2
      svf_fle = 'wvfn.save/'//trim(head_fle)//'_s.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# Table of non zero elements of the geminal matrix")')
      write(1,'(I6)') epo_c%n_nz_L
      do i1 = 1_int32, epo_c%n_nz_L
         write(1,'(2I6)') epo_c%nz_L(i1)%o1, epo_c%nz_L(i1)%o2
      enddo
      write(1,'("# Symmetry table of the geminal matrix")')
      write(1,'(I6)') epo_c%n_sym_L
      do i1 = 1_int32, epo_c%n_sym_L
         p1 => epo_c%sym_L(i1)%n_c
         write(1,'(1000I6)') p1,( sign(1,epo_c%sym_L(i1)%c(i2)) * epo_c%nz_L(abs(epo_c%sym_L(i1)%c(i2)))%o1, &
         & epo_c%nz_L(abs(epo_c%sym_L(i1)%c(i2)))%o2, i2 = 1_int32, p1)
      enddo
      close(1)
   end subroutine save_eppcfs_st
   subroutine read_eppcfs_st( head_fle )
      use openmp_mpi_m
      use geminal_symmetries_m, only: conv_indx
      use feporb_mod, only: n_eporbs
      character(4),  optional, intent(in)    :: head_fle
      character(100) :: svf_fle
      integer(int32) :: n_c
      integer(int32), allocatable, dimension(:) :: c
      integer(int32) :: o1
      integer(int32) :: i1, i2
      svf_fle = 'wvfn.save/' // trim(head_fle) // '_s.sav'
      if(mpi_rank.eq.0_int32) then
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*) !'("# Table of non zero elements of geminal")')
         read(1,*) epo_c%n_nz_L
         allocate( epo_c%nz_L(1:epo_c%n_nz_L) )
         do i1 = 1_int32, epo_c%n_nz_L
            read(1,*) epo_c%nz_L(i1)%o1, epo_c%nz_L(i1)%o2
            epo_c%nz_L(i1)%io = epo_c%nz_L(i1)%o1 + n_eporbs * ( epo_c%nz_L(i1)%o2 -1 )
         enddo
         allocate( c(1:2*epo_c%n_nz_L) )
         read(1,*) !'("# Symmetry table of geminal elements")')
         read(1,*) epo_c%n_sym_L
         allocate( epo_c%sym_L(1:epo_c%n_sym_L) )
         do i1 = 1_int32, epo_c%n_sym_L
            read(1,*) n_c, ( c(2*i2-1:2*i2), i2 = 1_int32, n_c )
            epo_c%sym_L(i1)%is_indp = .true.
            epo_c%sym_L(i1)%n_c     = n_c
            allocate( epo_c%sym_L(i1)%c(1:n_c) )
            do i2 = 1_int32, n_c
               o1 = abs(c(2*i2-1)) + n_eporbs * ( c(2*i2) - 1 )
               call conv_indx( epo_c%n_nz_L, epo_c%nz_L, o1 )
               epo_c%sym_L(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
            enddo 
         enddo
         deallocate( c )
         close(1)
      endif 
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(epo_c%n_nz_L, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      call mpi_bcast(epo_c%n_sym_L,1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      if ( mpi_rank.ne.0_int32 ) allocate( epo_c%nz_L(1:epo_c%n_nz_L))
      do i1 = 1_int32, epo_c%n_nz_L
         call mpi_bcast(epo_c%nz_L(i1)%o1, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(epo_c%nz_L(i1)%o2, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(epo_c%nz_L(i1)%io, 1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
      if ( mpi_rank.ne.0_int32 ) allocate( epo_c%sym_L(1:epo_c%n_sym_L) )
      do i1 = 1_int32, epo_c%n_sym_L
         call mpi_bcast(epo_c%sym_L(i1)%n_c,     1_int32,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
         call mpi_bcast(epo_c%sym_L(i1)%is_indp, 1_int32,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
         if ( mpi_rank.ne.0_int32 ) allocate( epo_c%sym_L(i1)%c(1:epo_c%sym_L(i1)%n_c) )
         call mpi_bcast(epo_c%sym_L(i1)%c(:),epo_c%sym_L(i1)%n_c,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
      enddo
#endif
   end subroutine read_eppcfs_st
   subroutine upd_eppcfs( n_par_gem_s, vec_epo_var )
      integer(int32) , intent(in) :: n_par_gem_s
      real(dp), dimension(n_par_gem_s), intent(in) :: vec_epo_var
      integer(int32), pointer :: o1, o2
      integer(int32)          :: i1, i2
      integer(int32)          :: ip
      ip = 1_int32
      do i1 = 1_int32, epo_c%n_sym_L
         do i2 = 1_int32, epo_c%sym_L(i1)%n_c
            o1 => epo_c%nz_L(abs(epo_c%sym_L(i1)%c(i2)))%o1
            o2 => epo_c%nz_L(abs(epo_c%sym_L(i1)%c(i2)))%o2
            epo_c%mat_L(o1,o2) = epo_c%mat_L(o1,o2) + dble(sign(1,epo_c%sym_L(i1)%c(i2))) * vec_epo_var(ip)
         enddo
         ip = ip + 1_int32
      enddo
   end subroutine upd_eppcfs
   subroutine save_eppcfs( head_fle )
      character(4),  optional, intent(in)    :: head_fle
      character(100)                         :: svf_fle
      integer(int32)                         :: i1, i2, i3
      svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# Ele-pos geminal coefficients")')
      i1 = 1_int32
      if (n_el*n_po.le.6_int32) then
         write(1,'(6I18)') (i2, i2 = 1_int32, n_el*n_po)
         do i2 = 1_int32, n_eporbs
            write(1,'(6E18.9)') epo_c%mat_L(i2,1:n_el*n_po)
         enddo
      else
         i3 = 1_int32
         do i1 = 1_int32, int(n_el*n_po/6_int32)
            write(1,'(6I18)') (i2, i2 = i3,i3+5_int32)
            i3 = i3 + 6_int32
            do i2 = 1_int32, n_eporbs
               write(1,'(6E18.9)') epo_c%mat_L(i2,6*(i1-1)+1:6*i1)
            enddo
         enddo
         if( mod(n_el*n_po, 6_int32).gt.0_int32 ) then
            write(1,'(6I18)') (i2, i2=i3,i3+mod(n_el*n_po, 6_int32)-1)
            do i2 = 1_int32, n_eporbs
               write(1,'(6E18.9)') epo_c%mat_L(i2,n_el*n_po-mod(n_el*n_po, 6_int32)+1:n_el*n_po)
            enddo
         endif
      endif
      close(1)
   end subroutine save_eppcfs
   subroutine read_eppcfs( head_fle )
      character(4),  optional,  intent(in)    :: head_fle
      character(100)                         :: svf_fle
      integer(int32)                         :: i1, i2
      svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
      if( mpi_rank.eq.0_int32 ) then
         open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
         read(1,*) ! Title
         i1 = 1_int32
         if (n_el*n_po.le.6_int32) then
            read(1,*) ! Numbers
            do i2 = 1_int32, n_eporbs
               read(1,*) epo_c%mat_L(i2,1:n_el*n_po)
            enddo
         else
            do i1 = 1_int32, int(n_el*n_po/6_int32)
               read(1,*) ! Numbers
               do i2 = 1_int32, n_eporbs
                  read(1,*) epo_c%mat_L(i2,6*(i1-1)+1:6*i1)
               enddo
            enddo
            if( mod(n_el*n_po, 6_int32).gt.0_int32 ) then
               read(1,*) ! Numbers
               do i2 = 1_int32, n_eporbs
                  read(1,*) epo_c%mat_L(i2,n_el*n_po-mod(n_el*n_po, 6_int32)+1:n_el*n_po)
               enddo
            endif
         endif
         close(1)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(epo_c%mat_L, n_eporbs*n_el*n_po, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
   end subroutine read_eppcfs
   subroutine sym_eppcfs( n_par_gem, n_par_gem_s, dl_ln_gem, dl_ln_gem_s )
      integer(int32),                    intent(in)    :: n_par_gem
      integer(int32),                    intent(in)    :: n_par_gem_s
      real(dp), dimension(n_par_gem),   intent(in)    :: dl_ln_gem
      real(dp), dimension(n_par_gem_s), intent(inout) :: dl_ln_gem_s
      integer(int32)                                    :: i1, i2, i3, i4
      i3 = 1_int32
      i4 = 1_int32
      do i1 = 1_int32, epo_c%n_sym_L
         dl_ln_gem_s(i3) = dl_ln_gem(abs(epo_c%sym_L(i1)%c(1)))
         if ( epo_c%sym_L(i1)%n_c.ge.2_int32 ) then
            do i2 = 2_int32, epo_c%sym_L(i1)%n_c
               dl_ln_gem_s(i3) = dl_ln_gem_s(i3) + &
               & dble(sign(1,epo_c%sym_L(i1)%c(i2))) * dl_ln_gem(abs(epo_c%sym_L(i1)%c(i2)))
            enddo
         endif
         i3 = i3 + 1_int32
      enddo
      i4 = i4 + epo_c%n_nz_L
   end subroutine sym_eppcfs
end module eppcfs_mod
