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
subroutine read_forbs_par( )
   use fortran_kinds_v, only: int32
   use openmp_mpi_m
   use write_lines_m
   use molecular_system_v, only: n_at, atoms
   use fermionic_orbitals_m, only: forbs
   use jastrow_orbitals_m, only: jorbs
   use primorbs_operations_m, only: cmp_orb_type, cmp_orb_ang
   implicit none
   integer(int32), external :: comp_atomic_charge
   character(2)   :: prm_type
   character(1)   :: orb_ang_mom
   integer(int32) :: n
   integer(int32) :: i1, i2, i3
   integer(int32) :: n_forbs, n_jorbs
   call write_simple_line(stdout,0,mpi_rank,2,"l","Reading atomic orbitals save file")
   if( mpi_rank.eq.0_int32 ) then
      open(unit=1,file='wvfn.save/frmbss.sav',action='read',form='formatted',status='old')
      read(1,*)
      do i1 = 1_int32, n_at
         read(1,*) atoms(i1)%atm_name, n_forbs, n_jorbs
         do i2 = 1_int32, n_forbs
            read(1,*) orb_ang_mom, forbs(i1)%orb(i2)%n
            forbs(i1)%orb(i2)%l = cmp_orb_ang (trim(orb_ang_mom))
            do i3 = 1_int32, forbs(i1)%orb(i2)%n
               read(1,*) forbs(i1)%orb(i2)%z(i3), forbs(i1)%orb(i2)%c(i3), prm_type
               call cmp_orb_type(prm_type, .true., forbs(i1)%orb(i2)%prm_t(i3), forbs(i1)%orb(i2)%prm_n(i3))
            enddo
         enddo
         if ( n_jorbs.gt.0_int32 ) then
            read(1,*)
            do i2 = 1_int32, n_jorbs
               read(1,*) orb_ang_mom, jorbs(i1)%orb(i2)%n
               jorbs(i1)%orb(i2)%l = cmp_orb_ang (trim(orb_ang_mom))
               do i3 = 1_int32, jorbs(i1)%orb(i2)%n
                  read(1,*) jorbs(i1)%orb(i2)%z(i3), jorbs(i1)%orb(i2)%c(i3), prm_type
                  call cmp_orb_type(prm_type, .false., jorbs(i1)%orb(i2)%prm_t(i3), jorbs(i1)%orb(i2)%prm_n(i3))
               enddo
            enddo
         endif
      enddo
      close(1)
   endif
#if defined _MPI || defined _MPIh || defined _MPI08
   do i1 = 1_int32, n_at
      do i2 = 1_int32, forbs(i1)%n_orbs
         n = int(forbs(i1)%orb(i2)%n)
         call mpi_bcast(forbs(i1)%orb(i2)%z, n, MPI_DOUBLE_PRECISION,&
         & 0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(forbs(i1)%orb(i2)%c, n, MPI_DOUBLE_PRECISION,&
         & 0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(forbs(i1)%orb(i2)%prm_n, n, MPI_INTEGER,&
         & 0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(forbs(i1)%orb(i2)%prm_t, n, MPI_INTEGER,&
         & 0,MPI_COMM_WORLD,mpierr)
      enddo
      if (jorbs(i1)%n_orbs.gt.0_int32) then
         do i2 = 1_int32, jorbs(i1)%n_orbs
            n = int(jorbs(i1)%orb(i2)%n)
            call mpi_bcast(jorbs(i1)%orb(i2)%z, n, MPI_DOUBLE_PRECISION,&
            & 0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(jorbs(i1)%orb(i2)%c, n, MPI_DOUBLE_PRECISION,&
            & 0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(jorbs(i1)%orb(i2)%prm_n, n, MPI_INTEGER,&
            & 0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(jorbs(i1)%orb(i2)%prm_t, n, MPI_INTEGER,&
            & 0,MPI_COMM_WORLD,mpierr)
         enddo
      endif
   enddo
   call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
end subroutine read_forbs_par
subroutine read_eporbs_par( )
   use fortran_kinds_v, only: int32
   use openmp_mpi_m
   use write_lines_m
   use feporb_mod, only: eporbs, eporbs_s, n_eporbs
   use jeporb_mod, only: jepos, jepos_s, n_jepos
   use primorbs_operations_m, only: cmp_orb_type, cmp_orb_ang
   implicit none
   character(2)   :: prm_type
   character(1)   :: orb_ang_mom
   integer(int32) :: n
   integer(int32) :: i1, i2
   integer(int32) :: n_forbs, n_jorbs
   logical        :: file_exists
   call write_simple_line(stdout,0,mpi_rank,2,"l","Reading positronic orbitals save file")
   if( mpi_rank.eq.0_int32 ) then
      inquire(file='wvfn.save/epobss.sav', exist=file_exists)
      if ( file_exists .and. (n_eporbs.gt.0_int32.or.n_jepos.gt.0_int32) ) then
         open(unit=1,file='wvfn.save/epobss.sav',action='read',form='formatted',status='old')
         read(1,*)
         read(1,*) n_forbs, n_jorbs
         do i1 = 1_int32, n_forbs
            read(1,*) orb_ang_mom, eporbs(i1)%n
            eporbs(i1)%l = cmp_orb_ang (trim(orb_ang_mom))
            do i2 = 1_int32, eporbs(i1)%n
               read(1,*) eporbs(i1)%z(i2), eporbs(i1)%c(i2), prm_type
               call cmp_orb_type(prm_type, .true., eporbs(i1)%prm_t(i2), eporbs(i1)%prm_n(i2))
            enddo
         enddo
         read(1,*)
         do i1 = 1_int32, n_jorbs
            read(1,*) orb_ang_mom, jepos(i1)%n
            jepos(i1)%l = cmp_orb_ang (trim(orb_ang_mom))
            do i2 = 1_int32, jepos(i1)%n
               read(1,*) jepos(i1)%z(i2), jepos(i1)%c(i2), prm_type
               call cmp_orb_type(prm_type, .false., jepos(i1)%prm_t(i2), jepos(i1)%prm_n(i2))
            enddo
         enddo
         close(1)
      else
         n_forbs = 0_int32
         n_jorbs = 0_int32
      endif
   endif
#if defined _MPI || defined _MPIh || defined _MPI08
   call mpi_bcast(n_forbs, 1, MPI_INTEGER, 0,MPI_COMM_WORLD,mpierr)
   call mpi_bcast(n_jorbs, 1, MPI_INTEGER, 0,MPI_COMM_WORLD,mpierr)
   if (n_forbs.gt.0_int32) then
      do i1 = 1_int32, eporbs_s%n_orbs
         n = int(eporbs(i1)%n)
         call mpi_bcast(eporbs(i1)%z, n, MPI_DOUBLE_PRECISION,&
         & 0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(eporbs(i1)%c, n, MPI_DOUBLE_PRECISION,&
         & 0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(eporbs(i1)%prm_n, n, MPI_INTEGER,&
         & 0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(eporbs(i1)%prm_t, n, MPI_INTEGER,&
         & 0,MPI_COMM_WORLD,mpierr)
      enddo
   endif
   if ( n_jorbs.gt.0_int32 ) then
      do i1 = 1_int32, jepos_s%n_orbs
         n = int(jepos(i1)%n)
         call mpi_bcast(jepos(i1)%z, n, MPI_DOUBLE_PRECISION,&
         & 0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(jepos(i1)%c, n, MPI_DOUBLE_PRECISION,&
         & 0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(jepos(i1)%prm_n, n, MPI_INTEGER,&
         & 0,MPI_COMM_WORLD,mpierr)
         call mpi_bcast(jepos(i1)%prm_t, n, MPI_INTEGER,&
         & 0,MPI_COMM_WORLD,mpierr)
      enddo
   endif
   call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
end subroutine read_eporbs_par
subroutine read_dorbs_par( )
   use fortran_kinds_v, only: int32
   use openmp_mpi_m
   use qdo_system_v, only: n_qdo, qdos
   use qdo_wavefunction_v, only: wf_type_d
   use drudonic_orbitals_m, only: dorbs
   use primorbs_operations_m, only: cmp_orb_type, cmp_orb_ang
   implicit none
   character(2)   :: prm_type
   character(1)   :: orb_ang_mom
   integer(int32) :: n
   integer(int32) :: i1, i2, i3
   if( mpi_rank.eq.0_int32 .and. wf_type_d.eq.5_int32) then
      open(unit=1,file='wvfn.save/qdobss.sav',action='read',form='formatted',status='old')
      read(1,*)
      do i1 = 1_int32, n_qdo
         read(1,*) qdos(i1)%qdo_name, n
         do i2 = 1_int32, dorbs(i1)%n_orbs
            read(1,*) orb_ang_mom, dorbs(i1)%orb(i2)%n
            dorbs(i1)%orb(i2)%l = cmp_orb_ang (trim(orb_ang_mom))
            do i3 = 1_int32, dorbs(i1)%orb(i2)%n
               read(1,*) dorbs(i1)%orb(i2)%z(i3), dorbs(i1)%orb(i2)%c(i3), prm_type
               call cmp_orb_type(prm_type, .true., dorbs(i1)%orb(i2)%prm_t(i3), dorbs(i1)%orb(i2)%prm_n(i3))
            enddo
         enddo 
      enddo 
      close(1)
   endif 
#if defined _MPI || defined _MPIh || defined _MPI08
   if (wf_type_d.eq.5_int32) then
      do i1 = 1_int32, n_qdo
         do i2 = 1_int32, dorbs(i1)%n_orbs
            n = int(dorbs(i1)%orb(i2)%n)
            call mpi_bcast(dorbs(i1)%orb(i2)%z, n, MPI_DOUBLE_PRECISION,&
            & 0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(dorbs(i1)%orb(i2)%c, n, MPI_DOUBLE_PRECISION,&
            & 0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(dorbs(i1)%orb(i2)%prm_n, n, MPI_INTEGER,&
            & 0,MPI_COMM_WORLD,mpierr)
            call mpi_bcast(dorbs(i1)%orb(i2)%prm_t, n, MPI_INTEGER,&
            & 0,MPI_COMM_WORLD,mpierr)
         enddo
      enddo 
   endif
   call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
end subroutine read_dorbs_par
subroutine save_forbs_par( )
   use fortran_kinds_v, only: int32
   use molecular_system_v, only: n_at, atoms
   use fermionic_orbitals_m, only: forbs
   use jastrow_orbitals_m, only: jorbs
   use primorbs_operations_m, only: prnt_orb_type, prnt_orb_ang
   implicit none
   character(2)   :: prm_type
   integer(int32) :: i1, i2, i3
   open(unit=1,file='wvfn.save/frmbss.sav',action='write',form='formatted',status='unknown')
   write(1,'("# Optimize atomic basis set")')
   do i1 = 1_int32, n_at
      write(1,'(A3,2I8)') atoms(i1)%atm_name, forbs(i1)%n_orbs, jorbs(i1)%n_orbs
      do i2 = 1_int32, forbs(i1)%n_orbs
         write(1,'(A1,I8)') prnt_orb_ang(forbs(i1)%orb(i2)%l), forbs(i1)%orb(i2)%n
         do i3 = 1_int32, forbs(i1)%orb(i2)%n
            call prnt_orb_type( forbs(i1)%orb(i2)%prm_t(i3), forbs(i1)%orb(i2)%prm_n(i3), prm_type )
            write(1,'(2E19.8,1X,A2)') forbs(i1)%orb(i2)%z(i3), forbs(i1)%orb(i2)%c(i3), prm_type
         enddo
      enddo
      if ( jorbs(i1)%n_orbs.gt.0_int32 ) then
         write(1,'("# Jastrow factor of the atom")')
         do i2 = 1_int32, jorbs(i1)%n_orbs
            write(1,'(A1,I8)') prnt_orb_ang(jorbs(i1)%orb(i2)%l), jorbs(i1)%orb(i2)%n
            do i3 = 1_int32, jorbs(i1)%orb(i2)%n
               call prnt_orb_type( jorbs(i1)%orb(i2)%prm_t(i3), jorbs(i1)%orb(i2)%prm_n(i3), prm_type )
               write(1,'(2E19.8,1X,A2)') jorbs(i1)%orb(i2)%z(i3), jorbs(i1)%orb(i2)%c(i3), prm_type
            enddo
         enddo
      endif
   enddo
   close(1)
end subroutine save_forbs_par
subroutine save_dorbs_par( )
   use fortran_kinds_v, only: int32
   use qdo_system_v, only: n_qdo, qdos
   use qdo_wavefunction_v, only: wf_type_d
   use drudonic_orbitals_m, only: dorbs
   use primorbs_operations_m, only: prnt_orb_type, prnt_orb_ang
   implicit none
   character(2)   :: prm_type
   integer(int32) :: i1, i2, i3
   if (wf_type_d.eq.5_int32) then
      open(unit=1,file='wvfn.save/qdobss.sav',action='write',form='formatted',status='unknown')
      write(1,'("# Optimized basis set of the QDO product wave function")')
      do i1 = 1_int32, n_qdo
         write(1,'(A3,2I8)') qdos(i1)%qdo_name, dorbs(i1)%n_orbs
         do i2 = 1_int32, dorbs(i1)%n_orbs
            write(1,'(A1,I8)') prnt_orb_ang(dorbs(i1)%orb(i2)%l), dorbs(i1)%orb(i2)%n
            do i3 = 1_int32, dorbs(i1)%orb(i2)%n
               call prnt_orb_type( dorbs(i1)%orb(i2)%prm_t(i3), dorbs(i1)%orb(i2)%prm_n(i3), prm_type )
               write(1,'(2E19.8,1X,A2)') dorbs(i1)%orb(i2)%z(i3), dorbs(i1)%orb(i2)%c(i3), prm_type
            enddo ! i3 (dorbs(i1)%orb(i2)%n )
         enddo ! i2 (dorbs(i1)%n_orbs)
      enddo ! i1 (n_qdo)
      close(1)
   endif ! (wf_type_d.eq.5_int32)
end subroutine save_dorbs_par
subroutine save_eporbs_par( )
   use fortran_kinds_v, only: int32
   use feporb_mod, only: eporbs, eporbs_s
   use jeporb_mod, only: jepos, jepos_s
   use primorbs_operations_m, only: prnt_orb_type, prnt_orb_ang
   implicit none
   character(2)   :: prm_type
   integer(int32) :: i1, i2
   if ( eporbs_s%n_orbs.gt.0 .or. jepos_s%n_orbs.gt.0 ) then
      open(unit=1,file='wvfn.save/epobss.sav',action='write',form='formatted',status='unknown')
      write(1,'("# Postronic basis set")')
      write(1,'(2I8)') eporbs_s%n_orbs, jepos_s%n_orbs
      do i1 = 1_int32, eporbs_s%n_orbs
         write(1,'(A1,I8)') prnt_orb_ang(eporbs(i1)%l), eporbs(i1)%n
         do i2 = 1_int32, eporbs(i1)%n
            call prnt_orb_type( eporbs(i1)%prm_t(i2), eporbs(i1)%prm_n(i2), prm_type )
            write(1,'(2E19.8,1X,A2)') eporbs(i1)%z(i2), eporbs(i1)%c(i2), prm_type
         enddo
      enddo
      write(1,'("# Jastrow Postronic basis set")')
      do i1 = 1_int32, jepos_s%n_orbs
         write(1,'(A1,I8)') prnt_orb_ang(jepos(i1)%l), jepos(i1)%n
         do i2 = 1_int32, jepos(i1)%n
            call prnt_orb_type( jepos(i1)%prm_t(i2), jepos(i1)%prm_n(i2), prm_type )
            write(1,'(2E19.8,1X,A2)') jepos(i1)%z(i2), jepos(i1)%c(i2), prm_type
         enddo
      enddo
      close(1)
   endif
end subroutine save_eporbs_par
pure function wvfn_cod( wf_name )
   use fortran_kinds_v, only: int32
   implicit none
   character(3),  intent(in)  :: wf_name
   integer(int32)             :: wvfn_cod
   select case ( wf_name )
    case('rsd')
      wvfn_cod =  1_int32
    case('usd')
      wvfn_cod = -1_int32
    case('rsg')
      wvfn_cod =  2_int32
    case('usg')
      wvfn_cod = -2_int32
    case('rtg')
      wvfn_cod =  3_int32
    case('utg')
      wvfn_cod = -3_int32
    case('rpf')
      wvfn_cod =  4_int32
    case('upf')
      wvfn_cod = -4_int32
    case('prd')
      wvfn_cod =  5_int32
    case('dip')
      wvfn_cod =  6_int32
    case('psg')
      wvfn_cod =  11_int32
    case default
      wvfn_cod =  0_int32
   end select ! wvfn_name
end function wvfn_cod
pure function wvfn_name( wf_code )
   use fortran_kinds_v, only: int32
   implicit none
   integer(int32), intent(in)  :: wf_code
   character(3)                :: wvfn_name
   select case ( wf_code )
    case(1_int32)
      wvfn_name = 'rsd'
    case(-1_int32)
      wvfn_name = 'usd'
    case(2_int32)
      wvfn_name = 'rsg'
    case(-2_int32)
      wvfn_name = 'usg'
    case(3_int32)
      wvfn_name = 'rtg'
    case(-3_int32)
      wvfn_name = 'utg'
    case(4_int32)
      wvfn_name = 'rpf'
    case(-4_int32)
      wvfn_name = 'upf'
    case(5_int32)
      wvfn_name = 'prd'
    case(6_int32)
      wvfn_name = 'dip'
    case(11_int32)
      wvfn_name = 'psg'
    case default
      wvfn_name = 'non'
   end select ! wf_code
end function wvfn_name
