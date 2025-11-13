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
module molecular_system_m
   use fortran_kinds_v,       only: int32, dp
   use openmp_mpi_m
   use physical_constants_v,  only: ang_to_au
   use write_lines_m,         only: write_separator_line, write_simple_line, &
   & write_empty_line, write_line_check, write_done
   use molecular_system_v
   use fermionic_config_c,    only: frm_cnf
   use quantum_monte_carlo_v, only: n_wlk_max
   use fermionic_positions_m, only: read_frmdst, deallocate_frmdst, frmdst
   use merge_sort_m,          only: mrg_srt_int
   use molec_symmetry_m,      only: atm_trs, n_sym_op, sym_op
   use external_field_m,      only: ext_E_field
   implicit none
   public  :: init_molec_sys, read_coord_fle, save_fermionic_confs, load_fermionic_confs
   private :: cmp_atms_dst, cmp_syst_sym, cmp_atms_sym
contains
   subroutine init_molec_sys
      integer(int32) :: i1, ierr
      call cmp_syst_sym()
      call cmp_atms_sym()
      nuc_chrg = 0_int32
      do i1 = 1_int32, n_at
         nuc_chrg = nuc_chrg + atoms(i1)%atm_z
      enddo ! i1 (n_at)
      ierr = 0
      if ( n_el.eq.0_int32.and.n_po.eq.0_int32 ) then
         n_el = nuc_chrg - chrg
         n_po = 0_int32
         if (mult.eq.-1_int32) then
            mult = 2_int32 * int(spin_e) + 1_int32
         else
            spin_e = dble(mult - 1_int32) / 2.0_dp
         endif
      else
         chrg = nuc_chrg - n_el
         mult = 2_int32 * int(spin_e) + 1_int32
      endif
      if ( spin_e .le. 10d-13 ) then
         if ( mod(dble(n_el), 2.0_dp) .gt. 10d-13 ) then
            call write_simple_line(stdout,0,mpi_rank,2,"l","ERROR!!! Number of electrons is odd while spin is zero")
            ierr = -1
         else
            n_el_u = n_el / 2_int32
            n_el_d = n_el_u
         endif
      else
         if ( abs(mod(dble(n_el) - 2.0_dp * spin_e, 2.0_dp))  .gt. 10d-13 ) then
            call write_simple_line(stdout,0,mpi_rank,2,"l","ERROR!!! Number of electrons is inconsistent with the multiplicity")
            ierr = -1
         else
            n_el_u = (n_el - int(2.0_dp * spin_e)) / 2_int32 + int(2.0_dp * spin_e)
            n_el_d = n_el - n_el_u
         endif
      endif
      n_el_s = int(2.0_dp * spin_e)
      if ( n_po.gt.0_int32 ) then
         if ( spin_p .le. 10d-13 ) then
            if ( mod(dble(n_po), 2.0_dp) .gt. 10d-13 ) then
               call write_simple_line(stdout,0,mpi_rank,2,"l","ERROR!!! Number of positrons is odd while spin is zero")
               ierr = -1
            else
               n_po_u = n_po / 2_int32
               n_po_d = n_po_u
            endif
         else
            if ( abs(mod(dble(n_po) - 2.0_dp * spin_p, 2.0_dp)) .gt. 10d-13 ) then
               call write_simple_line(stdout,0,mpi_rank,2,"l","ERROR!!! Number of positrons is inconsistent with the multiplicity")
               ierr = -1
            else
               n_po_u = ( n_po - int(2.0_dp * spin_p)) / 2_int32 + int(2.0_dp * spin_p)
               n_po_d = n_po - n_po_u
            endif
         endif
         n_po_s = int(2.0_dp * spin_p)
      endif
      n_fe = n_el + n_po
      if ( ierr .ne. 0_int32 ) then
         call fnlz_ompmpi_env()
         stop
      endif
      call write_separator_line(stdout,0,mpi_rank,2,"_")
      call write_simple_line(stdout,0,mpi_rank,2,"c","INITIALIZING MOLECULAR SYSTEM")
      call write_simple_line(stdout,0,mpi_rank,2,"c","________________________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","Molecular structure [Bohr]")
      call write_empty_line(stdout,0,mpi_rank)
      if(mpi_rank.eq.0) then
         write(stdout,'(5X,"Atom              x                   y                   z     ")')
         do i1 = 1_int32, n_at
            if (atoms(i1)%atm_name(1:1).eq.'*') then
               write(stdout,'(5X,A3,1X,3F20.8)') atoms(i1)%atm_name(1:3), r_at(1:3,i1)
            else
               write(stdout,'(6X,A2,2X,3F20.8)') atoms(i1)%atm_name(1:2), r_at(1:3,i1)
            endif
         enddo
      endif ! (mpi_rank.eq.0)
      call write_simple_line(stdout,0,mpi_rank,2,"c","________________________________")
      call write_simple_line(stdout,0,mpi_rank,2,"c","Electronic variables")
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Total charge of the system", nuc_chrg-n_el+n_po)
      call write_variable_line(stdout,0,mpi_rank,2,"Number of atoms in the system", n_at,var_name="n_at")
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of electrons", n_el,var_name="n_el")
      call write_variable_line(stdout,0,mpi_rank,2,"Number of spin up electrons", n_el_u,var_name="n_el_u")
      call write_variable_line(stdout,0,mpi_rank,2,"Number of spin down electrons", n_el_d,var_name="n_el_d")
      call write_variable_line(stdout,0,mpi_rank,2,"Total spin of the electrons", spin_e,var_name="spin_e")
      if(n_po.gt.0_int32) then
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_line(stdout,0,mpi_rank,2,"Total number of positrons", n_po,var_name="n_po")
         call write_variable_line(stdout,0,mpi_rank,2,"Number of spin up positrons", n_po_u,var_name="n_po_u")
         call write_variable_line(stdout,0,mpi_rank,2,"Number of spin down positrons", n_po_d,var_name="n_po_d")
         call write_variable_line(stdout,0,mpi_rank,2,"Total spin of the positrons", spin_p,var_name="spin_p")
      endif
      allocate( f_chrg(1:n_fe) ) ; f_chrg = 0.0
      allocate( f_spin(1:n_fe) ) ; f_spin = 0_int32
      do i1 = 1_int32, n_fe
         if ( i1.le.n_el ) then
            f_chrg(i1) = -1.0
            if ( i1.le.n_el_u ) then
               f_spin(i1) =  1_int32
            else
               f_spin(i1) = -1_int32
            endif
         else
            f_chrg(i1) =  1.0
            if ( i1-n_el.le.n_po_u ) then
               f_spin(i1) =  1_int32
            else
               f_spin(i1) = -1_int32
            endif
         endif ! ( i1.le.n_el )
      enddo ! ( n_fe )
      call cmp_atms_dst()
      if (n_at.gt.0_int32) call read_frmdst( )
      allocate( frm_cnf(1:n_wlk_max) )
      do i1 = 1_int32, n_wlk_max
         call frm_cnf(i1)%ini( i1 )
      enddo
      if (n_at.gt.0_int32) call deallocate_frmdst()
      mu_ep = m_e * m_p / ( m_e + m_p )
   end subroutine init_molec_sys
   subroutine load_fermionic_confs( )
      logical       :: file_present
      character(50) :: svf_fle
      real(dp),       allocatable, dimension(:,:) :: r_buf
      integer(int32), allocatable, dimension(:)   :: i_buf
      integer(int32) :: i_mpi_task, iw, icom
      file_present = .false.
      if (mpi_rank.eq.0_int32) then
         svf_fle = 'qmc.save/frm_conf.sav'
         inquire(file=svf_fle, exist=file_present)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(file_present, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif
      if (file_present) then
         call write_line_check(stdout,0,mpi_rank,2,"Loading fermionic configurations")
         if (mpi_rank.eq.0_int32) then
            open(unit=1,file=svf_fle,action='read',form='unformatted',status='old',access='sequential')
            do iw = 1_int32, n_wlk_max
               read(1) frm_cnf(iw)%r_fe(1:3,1:n_fe)
               read(1) frm_cnf(iw)%n_fe_tmv
               read(1) frm_cnf(iw)%i_fe_tmv(1:n_fe)
            enddo
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         if (n_mpi_tasks.gt.1_int32) then
            if (mpi_rank.eq.0_int32) then
               allocate( r_buf(1:3,1:n_fe) ) ; r_buf = 0.0_dp
               allocate( i_buf(1:n_fe) )     ; i_buf = 0_int32
            endif
            do i_mpi_task = 1_int32, n_mpi_tasks-1
               do iw = 1_int32, n_wlk_max
                  icom = n_wlk_max * (i_mpi_task-1) + iw
                  if (mpi_rank.eq.0_int32) then
                     read(1) r_buf(1:3,1:n_fe)
                     call mpi_send (r_buf(1:3,1:n_fe), 3*n_fe, MPI_DOUBLE_PRECISION, i_mpi_task, icom, MPI_COMM_WORLD, mpierr)
                  else if (mpi_rank.eq.i_mpi_task) then
                     call mpi_recv (frm_cnf(iw)%r_fe(1:3,1:n_fe), 3*n_fe, MPI_DOUBLE_PRECISION, 0, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
                  endif
                  call mpi_barrier(MPI_COMM_WORLD,mpierr)
                  if (mpi_rank.eq.0_int32) then
                     read(1) i_buf(1)
                     call mpi_send (i_buf(1), 1, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, mpierr)
                  else if (mpi_rank.eq.i_mpi_task) then
                     call mpi_recv (frm_cnf(iw)%n_fe_tmv, 1, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
                  endif
                  call mpi_barrier(MPI_COMM_WORLD,mpierr)
                  if (mpi_rank.eq.0_int32) then
                     read(1) i_buf(1:n_fe)
                     call mpi_send (i_buf(1:n_fe), n_fe, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, mpierr)
                  else if (mpi_rank.eq.i_mpi_task) then
                     call mpi_recv (frm_cnf(iw)%i_fe_tmv(1:n_fe), n_fe, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
                  endif
                  call mpi_barrier(MPI_COMM_WORLD,mpierr)
               enddo
            enddo
            if (mpi_rank.eq.0_int32) deallocate(r_buf,i_buf)
         endif
#endif
         call write_done(stdout,0,mpi_rank)
         if (mpi_rank.eq.0_int32) then
            close(1)
         endif
      endif
   end subroutine load_fermionic_confs
   subroutine save_fermionic_confs( )
      character(50)  :: svf_fle
      real(dp),       allocatable, dimension(:,:) :: r_buf
      integer(int32), allocatable, dimension(:)   :: i_buf
      integer(int32) :: i_mpi_task, iw, icom
      if (mpi_rank.eq.0_int32) then
         svf_fle = 'qmc.save/frm_conf.sav'
         open(unit=1,file=svf_fle,action='write',form='unformatted',status='unknown',access='sequential')
         do iw = 1_int32, n_wlk_max
            write(1) frm_cnf(iw)%r_fe(1:3,1:n_fe)
            write(1) frm_cnf(iw)%n_fe_tmv
            write(1) frm_cnf(iw)%i_fe_tmv(1:n_fe)
         enddo
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      if (n_mpi_tasks.gt.1_int32) then
         if (mpi_rank.eq.0_int32) then
            allocate( r_buf(1:3,1:n_fe) ) ; r_buf = 0.0_dp
            allocate( i_buf(1:n_fe) )     ; i_buf = 0_int32
         endif
         do i_mpi_task = 1_int32, n_mpi_tasks-1
            do iw = 1_int32, n_wlk_max
               icom = n_wlk_max * (i_mpi_task-1) + iw
               if (mpi_rank.eq.i_mpi_task) then
                  call mpi_send (frm_cnf(iw)%r_fe(1:3,1:n_fe), 3*n_fe, MPI_DOUBLE_PRECISION, 0, icom, MPI_COMM_WORLD, mpierr)
               else if (mpi_rank.eq.0) then
                  call mpi_recv (r_buf, 3*n_fe, MPI_DOUBLE_PRECISION, i_mpi_task, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
               endif
               call mpi_barrier(MPI_COMM_WORLD, mpierr)
               if (mpi_rank.eq.0_int32)  write(1) r_buf(1:3,1:n_fe)
               if (mpi_rank.eq.i_mpi_task) then
                  call mpi_send (frm_cnf(iw)%n_fe_tmv, 1, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, mpierr)
               else if (mpi_rank.eq.0) then
                  call mpi_recv (i_buf(1), 1, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
               endif
               call mpi_barrier(MPI_COMM_WORLD, mpierr)
               if (mpi_rank.eq.0_int32)  write(1) i_buf(1)
               if (mpi_rank.eq.i_mpi_task) then
                  call mpi_send (frm_cnf(iw)%i_fe_tmv(1:n_fe), n_fe, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, mpierr)
               else if (mpi_rank.eq.0) then
                  call mpi_recv (i_buf, n_fe, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
               endif
               call mpi_barrier(MPI_COMM_WORLD, mpierr)
               if (mpi_rank.eq.0_int32)  write(1) i_buf(1:n_fe)
            enddo
         enddo
         if (mpi_rank.eq.0_int32) deallocate(r_buf,i_buf)
      endif
#endif
      if ( mpi_rank.eq.0_int32 ) then
         close(1)
      endif
   end subroutine save_fermionic_confs
   subroutine read_coord_fle
      integer(int32)    :: ierr
      external          :: init_atomic_properties
      character(4)      :: units
      integer(int32) :: i1
      namelist /molsys/ n_at, n_el, n_po, spin_e, spin_p, m_p, m_e, chrg, mult,&
      & ext_E_field, units
      units  = 'bohr'
      chrg   =  0_int32
      mult   = -1_int32
      m_e    = 1.0_dp
      m_p    = 1.0_dp
      spin_e = 0.0_dp
      spin_p = 0.0_dp
      n_el   = 0_int32
      n_po   = 0_int32
      ext_E_field = 0.0_dp
      ierr = 0_int32
      if (trim(fle_coord).eq.'none') then
         n_fe = 0_int32
         n_el = 0_int32
         n_el_u = 0_int32
         n_el_d = 0_int32
         n_po = 0_int32
         n_po_u = 0_int32
         n_po_d = 0_int32
         molsys_prs = .false.
         return
      endif
      call write_empty_line(stdout,0,mpi_rank)
      call write_variable_line(stdout,0,mpi_rank,2,"Reading coordinates from file", fle_coord,var_name="file_coord")
      molsys_prs = .true.
      if(mpi_rank.eq.0_int32) then
         open(unit=1,file=fle_coord,action='read',form='formatted')
         read(1,nml=molsys,IOSTAT=ierr)
         call check_namelist_mandatory( ierr, "&molsys" )
         allocate( atoms(1:n_at), r_at(1:3,1:n_at) )
         do i1 = 1_int32, n_at
            atoms(i1)%atm_name = ""
            read(1,*) atoms(i1)%atm_name, r_at(1:3,i1)
         enddo
         if (units.eq.'angs') r_at(:,:) = r_at(:,:) * ang_to_au
         close(1)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_at, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(spin_e, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(spin_p, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(m_e,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(m_p,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_el, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(n_po, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(mult, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(chrg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(ext_E_field, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      if(mpi_rank.ne.0_int32) allocate( atoms(1:n_at), r_at(1:3,1:n_at) )
      do i1 = 1_int32, n_at
         call mpi_bcast(atoms(i1)%atm_name, 4, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      enddo
      call mpi_bcast(r_at, 3 * n_at, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
#endif
      call init_atomic_properties()
      if ( ierr .ne. 0_int32 ) then
         call fnlz_ompmpi_env()
         stop
      endif
   end subroutine read_coord_fle
   subroutine cmp_atms_dst
      integer(int32) :: i1, i2 ,i3
      allocate( d_nn(1:n_at*(n_at-1)/2) )
      if( n_at.ge.2_int32 ) then
         i3 = 0_int32
         do i2 = 1_int32, n_at ; do i1 = i2 + 1_int32, n_at
               i3 = i3 + 1_int32
               d_nn(i3)%v = r_at(:,i2) - r_at(:,i1)
               d_nn(i3)%m    = dsqrt(sum(d_nn(i3)%v(:)**2))
            enddo ; enddo
      endif
   end subroutine cmp_atms_dst
   subroutine cmp_syst_sym
      real(dp),        dimension(3,n_at)       :: r_ro
      real(dp),        dimension(3)              :: r_cc
      integer(int32)                               :: n_cnct_at
      integer(int32), dimension(n_at)            :: vec_cnct_at
      integer(int32), dimension(n_sym_op,0:n_at) :: sys_sym_tab_tmp
      logical,        dimension(n_at)            :: vec_indp_at
      integer(int32) :: i1, i2, is
      if( atm_trs ) then
         r_cc(:) = 0.0_dp
         do i1 = 1_int32, n_at
            if (atoms(i1)%atm_z.gt.0.0_dp) r_cc(:) = r_cc(:) + r_at(:,i1) * atoms(i1)%atm_z
         enddo
         if (sum(abs(r_cc(:))).gt.0.0_dp ) then
            r_cc(:) = r_cc(:) / dble(sum(atoms(:)%atm_z))
            do i1 = 1_int32, n_at
               r_at(:,i1) = r_at(:,i1) - r_cc(:)
            enddo
         endif
      endif
      if (sum(abs(ext_E_field)).gt.0.0_dp) then
         do i1 = 1_int32, n_at
            r_at(:,i1) = r_at(:,i1) + ext_E_field(:)
         enddo
      endif
      n_sys_sym = 0_int32
      do is = 1_int32, n_sym_op
         call dgemm('N', 'N', 3_int32, n_at, 3_int32, 1.0_dp, dble(sym_op(:,:,is)), 3_int32, &
         & r_at, 3_int32, 0.0_dp, r_ro, 3_int32)
         n_cnct_at   = 0_int32
         do i1 = 1_int32, n_at
            do i2 = 1_int32, n_at
               if( atoms(i1)%atm_z .eq. atoms(i2)%atm_z ) then
                  if( sum( abs(r_ro(:,i1)-r_at(:,i2)) ) .lt. 10.0d-13  ) then
                     vec_cnct_at(i1) = i2
                     n_cnct_at = n_cnct_at + 1_int32
                  endif
               endif
            enddo
         enddo
         if( n_cnct_at.eq.n_at ) then
            n_sys_sym = n_sys_sym + 1_int32
            sys_sym_tab_tmp(n_sys_sym,0) = is
            sys_sym_tab_tmp(n_sys_sym,1:n_at) = vec_cnct_at(:)
         endif
      enddo 
      vec_indp_at(:) = .true.
      do i1 = 1_int32, n_at
         if (vec_indp_at(i1) .eqv. .true.) then
            do i2 = 2_int32, n_sys_sym
               if ( sys_sym_tab_tmp(1,i1) .ne. sys_sym_tab_tmp(i2,i1) ) then
                  vec_indp_at(sys_sym_tab_tmp(i2,i1))=.false.
               endif
            enddo
         endif
      enddo
      n_indp_at = 0_int32
      do i1 = 1_int32, n_at
         if (vec_indp_at(i1) .eqv. .true.) n_indp_at = n_indp_at + 1_int32
      enddo
      if(.not.allocated(sys_sym_tab)) allocate( sys_sym_tab(1:n_sys_sym,0:n_at) )
      sys_sym_tab(1:n_sys_sym,0) = sys_sym_tab_tmp(1:n_sys_sym,0)
      i2 = 1_int32
      do i1 = 1_int32, n_at
         if (vec_indp_at(i1) .eqv. .true.) then
            sys_sym_tab(1:n_sys_sym,i2) = sys_sym_tab_tmp(1:n_sys_sym,i1)
            i2 = i2 + 1_int32
         endif
      enddo
      do i1 = 1_int32, n_at
         if (vec_indp_at(i1) .eqv. .false.) then
            sys_sym_tab(1:n_sys_sym,i2) = sys_sym_tab_tmp(1:n_sys_sym,i1)
            i2 = i2 + 1_int32
         endif
      enddo
      if (sum(abs(ext_E_field)).gt.0.0_dp) then
         do i1 = 1_int32, n_at
            r_at(:,i1) = r_at(:,i1) - ext_E_field(:)
         enddo 
      endif
   end subroutine cmp_syst_sym
   subroutine cmp_atms_sym
      integer(int32) :: n_tmp_cnct
      integer(int32), allocatable, dimension(:) :: tmp_cnct
      logical :: cnct_prs
      integer(int32) :: i1, i2, i3
      allocate( tmp_cnct(1:n_at) ) ; tmp_cnct = 0_int32
      allocate( atm_sym(1:n_indp_at) )
      do i1 = 1_int32, n_indp_at
         tmp_cnct(1) = sys_sym_tab(1,i1)
         n_tmp_cnct = 1_int32
         do i2 = 2_int32, n_sys_sym
            cnct_prs = .false.
            do i3 = 1_int32, n_tmp_cnct
               if ( sys_sym_tab(i2,i1) .eq. tmp_cnct(i3) ) then
                  cnct_prs = .true.
                  exit
               endif
            enddo 
            if (.not.cnct_prs) then
               n_tmp_cnct = n_tmp_cnct + 1_int32
               tmp_cnct(n_tmp_cnct) = sys_sym_tab(i2,i1)
            endif
         enddo 
         atm_sym(i1)%n_cnct_at = n_tmp_cnct
         allocate( atm_sym(i1)%cnct(1:n_tmp_cnct) )
         atm_sym(i1)%cnct(1:n_tmp_cnct) = tmp_cnct(1:n_tmp_cnct)
      enddo 
      do i1 = 1_int32, n_indp_at
         call mrg_srt_int( atm_sym(i1)%n_cnct_at, atm_sym(i1)%cnct(1:atm_sym(i1)%n_cnct_at) )
      enddo
      deallocate(tmp_cnct)
   end subroutine cmp_atms_sym
end module molecular_system_m
