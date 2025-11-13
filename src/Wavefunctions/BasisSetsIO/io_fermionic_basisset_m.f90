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
module io_fermionic_basisset_m
   use fortran_kinds_v, only: int32, dp
   use openmp_mpi_m
   use atomic_orbital_c, only: orb_t
   use fermionic_wavefunction_v, only: wf_type_e, wf_type_p, wf_type_ep
   use basisset_v, only: n_fbs, jst_bsst, frm_bsst, fle_basis
   use molecular_system_v, only: n_at, atoms, n_po
   use cubic_harmonics_m, only: l_max
   use primorbs_operations_m, only: cmp_orb_type, cmp_orb_ang, ini_orb_nrm, prnt_orb_type, prnt_orb_ang
   use cusps_functions_params_m, only: b_en, b_ee, b_pp, b_pn, b_ep
   implicit none
   public  :: read_bssst_fle, dall_bssst_var, prnt_bssst_fle
   private :: bcst_bssst_var
contains
   subroutine read_bssst_fle()
      character(1)   :: orb_l
      character(2)   :: orb_t
      character(3)   :: wf_name_e, wf_name_p
      character(4)   :: tmp_atm_name
      integer(int32), external :: wvfn_cod
      integer(int32), external :: comp_atomic_charge
      integer(int32) :: na_tmp, nj_tmp, nc_tmp
      integer(int32) :: i1, i2, i3, i4, io
      if (trim(fle_basis).ne.'none' ) then
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_line(stdout,0,mpi_rank,2,"Reading basis-set file", fle_basis,var_name="fle_basis")
         n_fbs = 0_int32
         if (mpi_rank.eq.0_int32) then
            open(unit=1,file=fle_basis,action='read',form='formatted')
            do i1 = 1_int32, 6_int32
               read(1,*)
            enddo
            do
               read(1,*,iostat=io)
               if (io.lt.0_int32) exit
               read(1,*,iostat=io) tmp_atm_name, na_tmp, nj_tmp
               if (io.lt.0_int32) exit
               n_fbs = n_fbs + 1_int32
               do i2 = 1_int32, na_tmp
                  read(1,*) orb_l, nc_tmp
                  do i3 = 1_int32, nc_tmp
                     read(1,*)
                  enddo
               enddo
               if ( nj_tmp .gt.0_int32) then
                  read(1,*)
                  do i2 = 1_int32, nj_tmp
                     read(1,*) orb_l, nc_tmp
                     do i3 = 1_int32, nc_tmp
                        read(1,*)
                     enddo
                  enddo
               endif
            enddo
            rewind(1)
            read(1,*)
            if ( n_po.gt.0_int32 ) then
               read(1,*) wf_name_e, wf_name_p
            else
               read(1,*) wf_name_e
            endif
            wf_type_e = wvfn_cod( trim(wf_name_e(1:3)) )
            if ( n_po.gt.0_int32 ) then
               wf_type_p = wvfn_cod( trim(wf_name_p(1:3)) )
            else
               wf_type_p = 0_int32
            endif
            if (abs(wf_type_p).ge.10_int32) then
               wf_type_ep = wf_type_p
               wf_type_p = 0_int32
            endif
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_bcast(n_fbs,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(wf_type_e,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(wf_type_p,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(wf_type_ep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
         allocate( frm_bsst(1:n_fbs), jst_bsst(1:n_fbs) )
         if (mpi_rank.eq.0_int32) then
            if ( n_po.eq.0_int32 ) then
               read(1,*) ! Read 1Body parameter
               read(1,*) b_en%typ, b_en%ord
               read(1,*) ! Read 2Body parameter
               read(1,*) b_ee%typ, b_ee%ord
            else
               read(1,*) ! Read 1Body parameter
               read(1,*) b_en%typ, b_en%ord, b_pn%typ, b_pn%ord
               read(1,*) ! Read 2Body parameter
               read(1,*) b_ee%typ, b_ee%ord, b_ep%typ, b_ep%ord
               b_pp%typ = b_ee%typ
               b_pp%ord = b_ee%ord
            endif
            if (n_fbs.gt.0_int32) then
               do i1 = 1_int32, n_fbs
                  read(1,*)  ! Comment line
                  read(1,*) frm_bsst(i1)%atm_name, frm_bsst(i1)%n_shlls, jst_bsst(i1)%n_shlls
                  jst_bsst(i1)%atm_name = frm_bsst(i1)%atm_name
                  frm_bsst(i1)%n_orbs = 0_int32
                  jst_bsst(i1)%n_orbs = 0_int32
                  if (frm_bsst(i1)%n_shlls.gt.0_int32) then
                     allocate( frm_bsst(i1)%shlls(1:frm_bsst(i1)%n_shlls) )
                     do i2 = 1, frm_bsst(i1)%n_shlls
                        read(1,*) orb_l, frm_bsst(i1)%shlls(i2)%n
                        frm_bsst(i1)%shlls(i2)%l = cmp_orb_ang( orb_l )
                        frm_bsst(i1)%n_orbs      = frm_bsst(i1)%n_orbs + 2_int32*frm_bsst(i1)%shlls(i2)%l + 1_int32
                        allocate(frm_bsst(i1)%shlls(i2)%c(1:frm_bsst(i1)%shlls(i2)%n))
                        allocate(frm_bsst(i1)%shlls(i2)%z(1:frm_bsst(i1)%shlls(i2)%n))
                        allocate(frm_bsst(i1)%shlls(i2)%prm_n(1:frm_bsst(i1)%shlls(i2)%n))
                        allocate(frm_bsst(i1)%shlls(i2)%prm_t(1:frm_bsst(i1)%shlls(i2)%n))
                        do i3 = 1_int32, frm_bsst(i1)%shlls(i2)%n
                           read(1,*) frm_bsst(i1)%shlls(i2)%z(i3), frm_bsst(i1)%shlls(i2)%c(i3), orb_t
                           call cmp_orb_type(orb_t, .true., frm_bsst(i1)%shlls(i2)%prm_t(i3), frm_bsst(i1)%shlls(i2)%prm_n(i3))
                        enddo
                     enddo
                  endif
                  if (jst_bsst(i1)%n_shlls.gt.0_int32) then
                     read(1,*) ! Jastrow Factor.
                     allocate( jst_bsst(i1)%shlls(1:jst_bsst(i1)%n_shlls) )
                     do i2 = 1, jst_bsst(i1)%n_shlls
                        read(1,*) orb_l, jst_bsst(i1)%shlls(i2)%n
                        jst_bsst(i1)%shlls(i2)%l = cmp_orb_ang( orb_l )
                        jst_bsst(i1)%n_orbs      = jst_bsst(i1)%n_orbs + 2_int32*jst_bsst(i1)%shlls(i2)%l + 1_int32
                        allocate(jst_bsst(i1)%shlls(i2)%c(1:jst_bsst(i1)%shlls(i2)%n))
                        allocate(jst_bsst(i1)%shlls(i2)%z(1:jst_bsst(i1)%shlls(i2)%n))
                        allocate(jst_bsst(i1)%shlls(i2)%prm_n(1:jst_bsst(i1)%shlls(i2)%n))
                        allocate(jst_bsst(i1)%shlls(i2)%prm_t(1:jst_bsst(i1)%shlls(i2)%n))
                        do i3 = 1_int32, jst_bsst(i1)%shlls(i2)%n
                           read(1,*) jst_bsst(i1)%shlls(i2)%z(i3), jst_bsst(i1)%shlls(i2)%c(i3), orb_t
                           call cmp_orb_type(orb_t, .false., jst_bsst(i1)%shlls(i2)%prm_t(i3), jst_bsst(i1)%shlls(i2)%prm_n(i3))
                        enddo
                     enddo
                  endif
               enddo
            else
               n_fbs = 1_int32
               l_max = 0_int32
               wf_type_e = 0
               wf_type_p = 0
               wf_type_ep = 0
               do i1 = 1_int32, n_fbs
                  frm_bsst(i1)%n_shlls = 0_int32 ; frm_bsst(i1)%n_orbs = 0_int32
                  jst_bsst(i1)%n_shlls = 0_int32 ; jst_bsst(i1)%n_orbs = 0_int32
               enddo
            endif
            close(1)
         endif
         call bcst_bssst_var()
         do i2 = 1_int32, n_fbs
            do i3 = 1_int32, frm_bsst(i2)%n_shlls
               if(frm_bsst(i2)%shlls(i3)%l.gt.l_max) l_max = frm_bsst(i2)%shlls(i3)%l
               do i4 = 1_int32, frm_bsst(i2)%shlls(i3)%n
                  call ini_orb_nrm( frm_bsst(i2)%shlls(i3)%prm_t(i4) )
               enddo
            enddo
            do i3 = 1_int32, jst_bsst(i2)%n_shlls
               if(jst_bsst(i2)%shlls(i3)%l.gt.l_max) l_max = jst_bsst(i2)%shlls(i3)%l
            enddo
         enddo
      else
         n_fbs = 1_int32
         l_max = 0_int32
         wf_type_e = 0
         wf_type_p = 0
         do i1 = 1_int32, n_fbs
            frm_bsst(i1)%n_shlls = 0_int32 ; frm_bsst(i1)%n_orbs = 0_int32
            jst_bsst(i1)%n_shlls = 0_int32 ; jst_bsst(i1)%n_orbs = 0_int32
         enddo
      endif
   end subroutine read_bssst_fle
   subroutine bcst_bssst_var( )
      integer(int32) :: n
      integer(int32) :: i1, i2
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast( b_en%typ, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast( b_en%ord, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast( b_ee%typ, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast( b_ee%ord, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast( b_pp%typ, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast( b_pp%ord, 1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast( b_pn%ord,    1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast( b_pn%typ,    1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast( b_ep%ord,    1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast( b_ep%typ,    1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      if (n_fbs.gt.0_int32) then
         do i1 = 1_int32, n_fbs
            call mpi_bcast(jst_bsst(i1)%atm_name, 4, MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(jst_bsst(i1)%n_shlls,  1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(jst_bsst(i1)%n_orbs,   1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(frm_bsst(i1)%atm_name, 4, MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(frm_bsst(i1)%n_shlls,  1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(frm_bsst(i1)%n_orbs,   1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            if(mpi_rank.ne.0_int32) then
               allocate( frm_bsst(i1)%shlls(1:frm_bsst(i1)%n_shlls) )
               allocate( jst_bsst(i1)%shlls(1:jst_bsst(i1)%n_shlls) )
            endif
            do i2 = 1_int32, frm_bsst(i1)%n_shlls
               call mpi_bcast(frm_bsst(i1)%shlls(i2)%l,     1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
               call mpi_bcast(frm_bsst(i1)%shlls(i2)%n,     1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
               n = int(frm_bsst(i1)%shlls(i2)%n)
               if( mpi_rank.ne.0_int32 ) then
                  allocate( frm_bsst(i1)%shlls(i2)%z(1:n),&
                  & frm_bsst(i1)%shlls(i2)%c(1:n), &
                  & frm_bsst(i1)%shlls(i2)%prm_n(1:n), &
                  & frm_bsst(i1)%shlls(i2)%prm_t(1:n) )
               endif
               call mpi_bcast(frm_bsst(i1)%shlls(i2)%z, n, MPI_DOUBLE_PRECISION,&
               & 0,MPI_COMM_WORLD,mpierr)
               call mpi_bcast(frm_bsst(i1)%shlls(i2)%c, n, MPI_DOUBLE_PRECISION,&
               & 0,MPI_COMM_WORLD,mpierr)
               call mpi_bcast(frm_bsst(i1)%shlls(i2)%prm_n, n, MPI_INTEGER,&
               & 0,MPI_COMM_WORLD,mpierr)
               call mpi_bcast(frm_bsst(i1)%shlls(i2)%prm_t, n, MPI_INTEGER,&
               & 0,MPI_COMM_WORLD,mpierr)
            enddo
            do i2 = 1_int32, jst_bsst(i1)%n_shlls
               call mpi_bcast(jst_bsst(i1)%shlls(i2)%l,     1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
               call mpi_bcast(jst_bsst(i1)%shlls(i2)%n,     1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
               n = int(jst_bsst(i1)%shlls(i2)%n)
               if( mpi_rank.ne.0_int32 ) then
                  allocate( jst_bsst(i1)%shlls(i2)%z(1:n), &
                  & jst_bsst(i1)%shlls(i2)%c(1:n), &
                  & jst_bsst(i1)%shlls(i2)%prm_n(1:n), &
                  & jst_bsst(i1)%shlls(i2)%prm_t(1:n) )
               endif
               call mpi_bcast(jst_bsst(i1)%shlls(i2)%z, n, MPI_DOUBLE_PRECISION,&
               & 0,MPI_COMM_WORLD,mpierr)
               call mpi_bcast(jst_bsst(i1)%shlls(i2)%c, n, MPI_DOUBLE_PRECISION,&
               & 0,MPI_COMM_WORLD,mpierr)
               call mpi_bcast(jst_bsst(i1)%shlls(i2)%prm_n, n, MPI_INTEGER,&
               & 0,MPI_COMM_WORLD,mpierr)
               call mpi_bcast(jst_bsst(i1)%shlls(i2)%prm_t, n, MPI_INTEGER,&
               & 0,MPI_COMM_WORLD,mpierr)
            enddo
         enddo
      endif
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
   end subroutine bcst_bssst_var
   subroutine dall_bssst_var()
      deallocate(frm_bsst,jst_bsst)
   end subroutine dall_bssst_var
   subroutine prnt_bssst_fle()
      character(2)   :: orb_t
      character(3),  external :: wvfn_name
      integer(int32) :: i1, i2, i3
      if (mpi_rank.eq.0_int32) then
         open(unit=1,file=fle_basis,action='write',form='formatted')
         write(1,'("# Wave function")')
         write(1,*) wvfn_name ( wf_type_e )
         write(1,'("# 1body type")') ! Read 1Body parameter
         write(1,'(2I4)') b_en%typ, b_en%ord
         write(1,'("# 2body type")') ! Read 2Body parameter
         write(1,'(2I4)') b_ee%typ, b_ee%ord
         do i1 = 1_int32, n_fbs
            write(1,'("# Fermionic basis set")')
            write(1,'(A3,2I4)') frm_bsst(i1)%atm_name, frm_bsst(i1)%n_shlls, jst_bsst(i1)%n_shlls
            do i2 = 1, frm_bsst(i1)%n_shlls
               write(1,'(A1,I4)') prnt_orb_ang( frm_bsst(i1)%shlls(i2)%l ), frm_bsst(i1)%shlls(i2)%n
               do i3 = 1_int32, frm_bsst(i1)%shlls(i2)%n
                  call prnt_orb_type( frm_bsst(i1)%shlls(i2)%prm_t(i3), frm_bsst(i1)%shlls(i2)%prm_n(i3), orb_t )
                  write(1,'(2E16.7,A3)') frm_bsst(i1)%shlls(i2)%z(i3), frm_bsst(i1)%shlls(i2)%c(i3), orb_t
               enddo
            enddo 
            if (jst_bsst(i1)%n_shlls.gt.0_int32) then
               write(1,'("# Jastrow basis set")')
               do i2 = 1, jst_bsst(i1)%n_shlls
                  write(1,'(A3,I4)') prnt_orb_ang( jst_bsst(i1)%shlls(i2)%l ), jst_bsst(i1)%shlls(i2)%n
                  do i3 = 1_int32, jst_bsst(i1)%shlls(i2)%n
                     call prnt_orb_type( jst_bsst(i1)%shlls(i2)%prm_t(i3), jst_bsst(i1)%shlls(i2)%prm_n(i3), orb_t )
                     write(1,'(2E16.7,A3)') jst_bsst(i1)%shlls(i2)%z(i3), jst_bsst(i1)%shlls(i2)%c(i3), orb_t
                  enddo
               enddo 
            endif
         enddo 
         close(1)
      endif 
   end subroutine prnt_bssst_fle
end module io_fermionic_basisset_m
