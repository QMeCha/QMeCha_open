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
module pseudopotentials_v
   use fortran_kinds_v,       only: dp, int32, stdout
   use openmp_mpi_m
   use write_lines_m,         only: write_separator_line, write_simple_line, write_empty_line
   use quantum_monte_carlo_v, only: qmc_mthd
   use nl_pseudo_grid_m,      only: n_psd_grd_pnts
   use molecular_system_v,    only: n_at, atoms
   implicit none
   type, public :: psdpot_t
      character(2)   :: atm_name
      integer(int32) :: n_psd_c
      integer(int32) :: Z_rep
      integer(int32) :: tot_psd_ele
      real(dp)       :: r_cut
      integer(int32), allocatable, dimension(:) :: n_psd_ele
      real(dp),       allocatable, dimension(:) :: c
      real(dp),       allocatable, dimension(:) :: e
      integer(int32), allocatable, dimension(:) :: n
   end type psdpot_t
   logical,        save, public                            :: psd_prs
   character(20),  save, public                            :: fle_pspot
   integer(int32), save, public                            :: n_psdpot
   type(psdpot_t), save, public, allocatable, dimension(:) :: psdpot
   real(dp),       save, public                            :: ene_psd_cut
   integer(int32), save, public                            :: psd_int_mthd
   public  :: read_psdpot_fle , prnt_psdpot
   private :: count_pseudo_atoms
contains
   subroutine read_psdpot_fle()
      integer(int32) :: ierr
      integer(int32) :: n_pseudo
      integer(int32) :: i1, i2, i3
      logical        :: add_pseudo, repeated_pseudo
      character(3)   :: atm_name_t
      integer(int32) :: n_psd_c_t
      integer(int32) :: Z_rep_t 
      integer(int32) :: n_psd_ele_t(1:7)
      namelist /pseudo/ ene_psd_cut, n_pseudo, n_psd_grd_pnts, psd_int_mthd
      call count_pseudo_atoms()
      n_pseudo = n_psdpot
      ene_psd_cut = 1.0d-5
      n_psd_grd_pnts = 6_int32
      psd_int_mthd = 0_int32
      ierr = 0_int32
      if ( n_psdpot.gt.0_int32 .and. trim(fle_pspot).ne.'none' ) then
         allocate(psdpot(1:n_psdpot))
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_line(stdout,0,mpi_rank,2,"Reading pseudopotential file", trim(fle_pspot),var_name="fle_pspot")
         if(mpi_rank.eq.0_int32) then
            open(unit=1,file=fle_pspot,action='read',form='formatted')
            read(1,nml=pseudo,IOSTAT=ierr)
            call check_namelist_mandatory( ierr, "&pseudo" )
            if (n_pseudo.lt.n_psdpot) then
               call write_simple_line(stdout,0,mpi_rank,2,"l","ERROR!!! Number of pseudos do not correspond to those in the molecular system.")
               ierr = -1
            endif
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_bcast(ierr,  1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(n_psd_grd_pnts,  1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(ene_psd_cut,     1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(psd_int_mthd,    1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
#endif
         if (ierr.eq.0_int32) then
            i2 = 1_int32 
            if(mpi_rank.eq.0_int32) then
               do i1 = 1_int32, n_pseudo
                  read(1,*) atm_name_t, n_psd_c_t, Z_rep_t 
                  call chck_pseudo_atom( atm_name_t, i2, add_pseudo )
                  if ( add_pseudo ) then 
                     psdpot(i2)%atm_name = atm_name_t
                     psdpot(i2)%n_psd_c = n_psd_c_t
                     psdpot(i2)%Z_rep = Z_rep_t
                     allocate( psdpot(i2)%n_psd_ele(1:psdpot(i2)%n_psd_c) )
                     read(1,*) psdpot(i2)%n_psd_ele(1:psdpot(i2)%n_psd_c)
                     psdpot(i2)%tot_psd_ele = sum(psdpot(i2)%n_psd_ele(1:psdpot(i2)%n_psd_c))
                     allocate( psdpot(i2)%c(1:psdpot(i2)%tot_psd_ele) )
                     allocate( psdpot(i2)%e(1:psdpot(i2)%tot_psd_ele) )
                     allocate( psdpot(i2)%n(1:psdpot(i2)%tot_psd_ele) )
                     do i3 = 1_int32, psdpot(i2)%tot_psd_ele
                        read(1,*) psdpot(i2)%n(i3), psdpot(i2)%e(i3), psdpot(i2)%c(i3)
                     enddo 
                     i2 = i2 + 1_int32
                  else 
                     read(1,*) n_psd_ele_t(1:n_psd_c_t) 
                     do i3 = 1_int32, sum( n_psd_ele_t(1:n_psd_c_t) )
                        read(1,*) 
                     enddo 
                  endif
               enddo 
            endif 
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         do i1 = 1_int32, n_psdpot
            call mpi_bcast(psdpot(i1)%Z_rep,       1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast(psdpot(i1)%n_psd_c,     1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast(psdpot(i1)%atm_name,    2, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast(psdpot(i1)%tot_psd_ele, 1, MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
            if ( mpi_rank.ne.0_int32 ) allocate( psdpot(i1)%n_psd_ele(1:psdpot(i1)%n_psd_c) )
            call mpi_bcast( psdpot(i1)%n_psd_ele(1:psdpot(i1)%n_psd_c),psdpot(i1)%n_psd_c, &
            & MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
            if ( mpi_rank.ne.0_int32 ) then
               allocate( psdpot(i1)%c(1:psdpot(i1)%tot_psd_ele) )
               allocate( psdpot(i1)%e(1:psdpot(i1)%tot_psd_ele) )
               allocate( psdpot(i1)%n(1:psdpot(i1)%tot_psd_ele) )
            endif
            call mpi_bcast( psdpot(i1)%c(1:psdpot(i1)%tot_psd_ele), psdpot(i1)%tot_psd_ele,&
            & MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast( psdpot(i1)%e(1:psdpot(i1)%tot_psd_ele), psdpot(i1)%tot_psd_ele,&
            & MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast( psdpot(i1)%n(1:psdpot(i1)%tot_psd_ele), psdpot(i1)%tot_psd_ele,&
            & MPI_INTEGER,  0, MPI_COMM_WORLD, mpierr)
         enddo
#endif
         if (mpi_rank.eq.0_int32) close(1)
      endif
      if ( ierr .ne. 0_int32 ) then
         call fnlz_ompmpi_env()
         stop
      endif
   end subroutine read_psdpot_fle
   subroutine prnt_psdpot( )
      integer(int32) :: i1
      call write_simple_line(stdout,0,mpi_rank,2,"c","________________________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","List of pseudopotentials")
      call write_empty_line(stdout,0,mpi_rank)
      if (mpi_rank.eq.0) then
         write(*,'(5X,"Atom       N° Non-loc. c.     Charge Replaced      R_cut [bohr] ")')
         do i1 = 1_int32, n_psdpot
            write(*,'(5X,A2,I17,I19,8X,F14.4)') psdpot(i1)%atm_name, psdpot(i1)%n_psd_c-1, psdpot(i1)%Z_rep, psdpot(i1)%r_cut
         enddo
      endif
   end subroutine prnt_psdpot
   subroutine count_pseudo_atoms()
      logical        :: diff_pspot
      integer(int32) :: i1, i2
      psd_prs = .false.
      do i1 = 1_int32, n_at
         if (atoms(i1)%i_pp.eq.-1) psd_prs = .true.
      enddo 
      if (.not.psd_prs) then
         n_psdpot = 0_int32
      else
         if (atoms(1)%i_pp.eq.-1_int32) then
            n_psdpot = 1_int32
         else
            n_psdpot = 0_int32
         endif
         do i1 = 2_int32, n_at
            if (atoms(i1)%i_pp.eq.-1_int32) then
               diff_pspot = .true.
               do i2 = 1_int32, (i1 - 1_int32)
                  if (atoms(i2)%i_pp.eq.-1_int32) then
                     if(trim(atoms(i1)%atm_name(2:3)).eq.trim(atoms(i2)%atm_name(2:3))) then
                        diff_pspot = .false.
                        cycle
                     endif
                  endif
               enddo
               if (diff_pspot) n_psdpot = n_psdpot + 1_int32
            endif
         enddo 
      endif
   end subroutine count_pseudo_atoms
   subroutine chck_pseudo_atom( atom_name, n_tmp_pseudos, add_pseudo )
      character(2),    intent(in)  :: atom_name
      integer(int32),  intent(in)  :: n_tmp_pseudos   
      logical,         intent(out) :: add_pseudo   
      integer(int32) :: i1
      logical :: repeated_pseudo
      repeated_pseudo = .false.
      add_pseudo      = .false.
      if (n_tmp_pseudos.gt.1_int32) then 
         do i1 = 1_int32, n_tmp_pseudos-1
            if ( psdpot(i1)%atm_name .eq. atom_name ) then
               repeated_pseudo =.true. 
               exit 
            endif 
         enddo
      endif 
      if (repeated_pseudo) then
         call write_simple_line(stdout,0,mpi_rank,2,"l","WARNING!!! Repeated pseudo for atom "//trim(atom_name)//" check pseudo.dat file.")
         return 
      endif
      do i1 = 1_int32, n_at
         if (atoms(i1)%i_pp.eq.-1_int32) then
            if(trim(atoms(i1)%atm_name(2:3)).eq.atom_name) then
               add_pseudo = .true.
               exit
            endif
         endif
      enddo            
   end subroutine chck_pseudo_atom
end module pseudopotentials_v
