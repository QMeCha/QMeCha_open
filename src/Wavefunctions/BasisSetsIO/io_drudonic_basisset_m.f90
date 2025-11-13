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
module io_drudonic_basisset_m
   use fortran_kinds_v, only: int8, int32, dp
   use openmp_mpi_m
   use atomic_orbital_c, only: orb_t
   use qdo_wavefunction_v, only: wf_type_d
   use basisset_v, only: n_qbs, drd_bsst, fle_basis_qdo
   use jstqepar_var, only: jstqe_prs, jstqe_type
   use qdo_system_v, only: n_qdo, qdos
   use cubic_harmonics_m, only: l_max
   use primorbs_operations_m, only: cmp_orb_type, cmp_orb_ang, ini_orb_nrm
   use qdo_jastrow_m, only: b_dd, b_qd
   implicit none
   integer(int32),  external :: wvfn_cod
   public  :: read_qbssst_fle
#if defined _MPI || defined _MPIh || defined _MPI08
   private :: bcst_prdbssst_var
#endif
contains
   subroutine read_qbssst_fle()
      logical        :: diff_basis
      character(1)   :: orb_l
      character(2)   :: orb_t
      character(3)             :: wf_name
      integer(int32) :: i1, i2, i3, i4
      if (mpi_rank.eq.0_int32) then
         open(unit=1,file=fle_basis_qdo,action='read',form='formatted')
         read(1,*) ! Comment line
         read(1,*) wf_name
         read(1,*) ! Comment line
         read(1,*) jstqe_type
         read(1,*) ! Comment line
         read(1,*) b_qd%typ, b_qd%ord
         read(1,*) ! Comment line
         read(1,*) b_dd%typ, b_dd%ord
         wf_type_d = wvfn_cod( wf_name )
         if(jstqe_type.ne.0) jstqe_prs = .true.
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(wf_type_d,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jstqe_type,  1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(jstqe_prs,   1, MPI_LOGICAL,  0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(b_qd%typ,    1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(b_dd%typ,    1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(b_qd%ord,    1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(b_dd%ord,    1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
      if ( abs(wf_type_d) .eq.5 ) then
         n_qbs = 1_int32
         do i1 = 2_int32, n_qdo
            diff_basis = .true.
            do i2 = 1_int32, (i1 - 1_int32)
               if(trim(qdos(i1)%qdo_name).eq.trim(qdos(i2)%qdo_name)) then
                  diff_basis = .false.
                  cycle
               endif
            enddo
            if (diff_basis) n_qbs = n_qbs + 1_int32
         enddo
         allocate( drd_bsst(1:n_qbs) )
         if (mpi_rank.eq.0_int32) then
            do i1 = 1_int32, n_qbs
               read(1,*)  ! Comment line
               read(1,*) drd_bsst(i1)%atm_name, drd_bsst(i1)%n_shlls
               drd_bsst(i1)%n_orbs = 0_int32
               allocate( drd_bsst(i1)%shlls(1:drd_bsst(i1)%n_shlls) )
               do i2 = 1, drd_bsst(i1)%n_shlls
                  read(1,*) orb_l, drd_bsst(i1)%shlls(i2)%n
                  drd_bsst(i1)%shlls(i2)%l = cmp_orb_ang( orb_l )
                  drd_bsst(i1)%n_orbs      = drd_bsst(i1)%n_orbs + 2_int32*drd_bsst(i1)%shlls(i2)%l + 1_int32
                  allocate(drd_bsst(i1)%shlls(i2)%c(1:drd_bsst(i1)%shlls(i2)%n))
                  allocate(drd_bsst(i1)%shlls(i2)%z(1:drd_bsst(i1)%shlls(i2)%n))
                  allocate(drd_bsst(i1)%shlls(i2)%prm_n(1:drd_bsst(i1)%shlls(i2)%n))
                  allocate(drd_bsst(i1)%shlls(i2)%prm_t(1:drd_bsst(i1)%shlls(i2)%n))
                  do i3 = 1_int32, drd_bsst(i1)%shlls(i2)%n
                     read(1,*) drd_bsst(i1)%shlls(i2)%z(i3), drd_bsst(i1)%shlls(i2)%c(i3), orb_t
                     call cmp_orb_type(orb_t, .true., drd_bsst(i1)%shlls(i2)%prm_t(i3), drd_bsst(i1)%shlls(i2)%prm_n(i3))
                  enddo
               enddo
            enddo
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         call bcst_prdbssst_var()
#endif
         do i1 = 1_int32, n_qdo ; do i2 = 1_int32, n_qbs
               if(trim(drd_bsst(i2)%atm_name).eq.trim(qdos(i1)%qdo_name)) then
                  do i3 = 1_int32, drd_bsst(i2)%n_shlls
                     if(drd_bsst(i2)%shlls(i3)%l.gt.l_max) l_max = drd_bsst(i2)%shlls(i3)%l
                     do i4 = 1_int32, drd_bsst(i2)%shlls(i3)%n
                        call ini_orb_nrm( drd_bsst(i2)%shlls(i3)%prm_t(i4) )
                     enddo
                  enddo
               endif
            enddo ; enddo
      endif
      if (mpi_rank.eq.0_int32) close(1)
   end subroutine read_qbssst_fle
#if defined _MPI || defined _MPIh || defined _MPI08
   subroutine bcst_prdbssst_var( )
      integer(int32) :: n
      integer(int32) :: i1, i2
      do i1 = 1_int32, n_qbs
         call mpi_bcast(drd_bsst(i1)%atm_name, 4, MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(drd_bsst(i1)%n_shlls,  1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(drd_bsst(i1)%n_orbs,   1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
         if(mpi_rank.ne.0_int32) allocate( drd_bsst(i1)%shlls(1:drd_bsst(i1)%n_shlls) )
         do i2 = 1_int32, drd_bsst(i1)%n_shlls
            call mpi_bcast(drd_bsst(i1)%shlls(i2)%l,     1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(drd_bsst(i1)%shlls(i2)%n,     1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            n = int(drd_bsst(i1)%shlls(i2)%n)
            if( mpi_rank.ne.0_int32 ) then
               allocate( drd_bsst(i1)%shlls(i2)%z(1:n),&
               & drd_bsst(i1)%shlls(i2)%c(1:n), &
               & drd_bsst(i1)%shlls(i2)%prm_n(1:n), &
               & drd_bsst(i1)%shlls(i2)%prm_t(1:n) )
            endif
            call mpi_bcast(drd_bsst(i1)%shlls(i2)%z, n, MPI_DOUBLE_PRECISION,&
            & 0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(drd_bsst(i1)%shlls(i2)%c, n, MPI_DOUBLE_PRECISION,&
            & 0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(drd_bsst(i1)%shlls(i2)%prm_n, n, MPI_INTEGER,&
            & 0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(drd_bsst(i1)%shlls(i2)%prm_t, n, MPI_INTEGER,&
            & 0,MPI_COMM_WORLD,mpierr)
         enddo
      enddo 
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
   end subroutine bcst_prdbssst_var
#endif
end module io_drudonic_basisset_m
