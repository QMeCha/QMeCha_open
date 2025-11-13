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
module qdo_system_m
   use fortran_kinds_v,       only: int8, int32, dp
   use physical_constants_v,  only: ang_to_au
   use openmp_mpi_m
   use external_field_m,      only: ext_E_field
   use write_lines_m,         only: write_separator_line, write_simple_line, &
   & write_empty_line, write_line_check, write_done
   use qdo_system_v
   use drudonic_config_c,     only: drd_cnf
   use molecular_system_v,    only: n_at, r_at
   use quantum_monte_carlo_v, only: n_wlk_max
   use molecular_system_v,    only: molsys_prs, n_fe
   implicit none
   public  :: ini_qdo_sys, read_qcoord_fle, save_drudonic_confs, load_drudonic_confs
   private :: cmp_qdo_dist, cmp_qdonuc_dist, cmp_pchrg_dist, cmp_qdopchrg_dist, cmp_nucpchrg_dist
contains
   subroutine ini_qdo_sys
      integer(int32) :: iw, i1
      if(mpi_rank.eq.0) then
         write(*,'(3X,"================================================================")')
         write(*,'(3X,"                       QDOs SYSTEM                         ")')
         write(*,'(3X,"               ________________________________                 ")')
         write(*,'(3X,"                  QDOs structure [Bohr]                    ")')
         call write_empty_line(stdout,0,mpi_rank)
         write(*,'(3X,"Reading QDO coordinates from file   (file_coord_qdo) : ",A30)') fle_coord_qdo
         call write_empty_line(stdout,0,mpi_rank)
         write(*,'(3X,"QDO              x                   y                   z     ")')
         do i1 = 1_int32, n_qdo
            write(*,'(3X,A2,2X,3F20.8)') qdos(i1)%qdo_name(1:2), r_qdo(1:3,i1)
         enddo
         call write_empty_line(stdout,0,mpi_rank)
         write(*,'(3X,"QDO              q                   omega               m     ")')
         do i1 = 1_int32, n_qdo
            write(*,'(3X,A2,2X,3F20.8)') qdos(i1)%qdo_name(1:2), qdos(i1)%qdo_q, qdos(i1)%qdo_omega, qdos(i1)%qdo_m
         enddo
         call write_empty_line(stdout,0,mpi_rank)
         write(*,'(3X,"QDO              sigma_c             sigma_d")')
         do i1 = 1_int32, n_qdo
            write(*,'(3X,A2,2X,3F20.8)') qdos(i1)%qdo_name(1:2), qdos(i1)%qdo_sigma_c, qdos(i1)%qdo_sigma_d
         enddo
         call write_empty_line(stdout,0,mpi_rank)
         if( n_fe.gt.0_int32 ) then
            write(*,'(3X,"Sigma for electrons                                 : ",F10.4)') el_sigma
         endif
         write(*,'(3X,"Total number of QDOs/Drudons                        : ",I10)') n_qdo
         if( n_pchrg.gt.0_int32 ) then
            call write_empty_line(stdout,0,mpi_rank)
            write(*,'(3X,"================================================================")')
            write(*,'(3X,"                       POINT CHARGES                         ")')
            write(*,'(3X,"               ________________________________                 ")')
            write(*,'(3X,"              point charges structure [Bohr]                    ")')
            call write_empty_line(stdout,0,mpi_rank)
            write(*,'(3X,"pchrg            x                   y                   z     ")')
            do i1 = 1_int32, n_pchrg
               write(*,'(3X,A2,2X,3F20.8)') pchrgs(i1)%pchrg_name, r_pchrg(1:3,i1)
            enddo
            call write_empty_line(stdout,0,mpi_rank)
            write(*,'(3X,"pchrg            q                   sigma                      ")')
            do i1 = 1_int32, n_pchrg
               write(*,'(3X,A2,2X,3F20.8)') pchrgs(i1)%pchrg_name, pchrgs(i1)%pchrg_q, pchrgs(i1)%pchrg_sigma
            enddo
            call write_empty_line(stdout,0,mpi_rank)
            write(*,'(3X,"Total number of point charges                       : ",I10)') n_pchrg
         endif
      endif ! (mpi_rank.eq.0)
      allocate( d_chrg(1:n_qdo) ) ; d_chrg = 0_dp
      allocate( d_mass(1:n_qdo) ) ; d_mass = 0_dp
      do i1 = 1_int32, n_qdo
         d_chrg(i1) = -qdos(i1)%qdo_q
         d_mass(i1) =  qdos(i1)%qdo_m
      enddo ! ( n_qdo )
      call cmp_qdo_dist()
      call cmp_pchrg_dist()
      call cmp_qdopchrg_dist()
      if (molsys_prs) then
         call cmp_qdonuc_dist()
         call cmp_nucpchrg_dist()
      endif
      allocate( drd_cnf(1:n_wlk_max) )
      do iw = 1_int32, n_wlk_max
         call drd_cnf(iw)%ini( iw )
      enddo
      call load_drudonic_confs() 
   end subroutine ini_qdo_sys
   subroutine load_drudonic_confs( )
      logical        :: file_present
      character(50)  :: svf_fle
      real(dp),       allocatable, dimension(:,:) :: r_buf
      integer(int32), allocatable, dimension(:)   :: i_buf
      integer(int32) :: i_mpi_task, iw, icom
      file_present = .false.
      if (mpi_rank.eq.0_int32) then
         svf_fle = 'qmc.save/drd_conf.sav'
         inquire(file=svf_fle, exist=file_present)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(file_present, 1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif
      if (file_present) then
         call write_line_check(stdout,0,mpi_rank,2,"Loading drudonic configurations")
         if (mpi_rank.eq.0_int32) then
            open(unit=1,file=svf_fle,action='read',form='unformatted',status='old',access='sequential')
            do iw = 1_int32, n_wlk_max
               read(1) drd_cnf(iw)%r_drd(1:3,1:n_qdo)
               read(1) drd_cnf(iw)%n_drd_tmv
               read(1) drd_cnf(iw)%i_drd_tmv(1:n_qdo)
            enddo
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         if (n_mpi_tasks.gt.1_int32) then
            if (mpi_rank.eq.0_int32) then 
               allocate( r_buf(1:3,1:n_qdo) ) ; r_buf = 0.0_dp
               allocate( i_buf(1:n_qdo) ) ; i_buf = 0_int32
            endif
            do i_mpi_task = 1_int32, n_mpi_tasks-1
               do iw = 1_int32, n_wlk_max
                  icom = n_wlk_max * (i_mpi_task-1) + iw
                  if (mpi_rank.eq.0_int32) then
                     read(1) r_buf(1:3,1:n_qdo)
                     call mpi_send (r_buf(1:3,1:n_qdo), 3*n_qdo, MPI_DOUBLE_PRECISION, i_mpi_task, icom, MPI_COMM_WORLD, mpierr)
                  else if (mpi_rank.eq.i_mpi_task) then
                     call mpi_recv (drd_cnf(iw)%r_drd(1:3,1:n_qdo), 3*n_qdo, MPI_DOUBLE_PRECISION, 0, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
                  endif
                  call mpi_barrier(MPI_COMM_WORLD,mpierr)
                  if (mpi_rank.eq.0_int32) then
                     read(1) i_buf(1)
                     call mpi_send (i_buf(1), 1, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, mpierr)
                  else if (mpi_rank.eq.i_mpi_task) then
                     call mpi_recv (drd_cnf(iw)%n_drd_tmv, 1, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
                  endif
                  call mpi_barrier(MPI_COMM_WORLD,mpierr)
                  if (mpi_rank.eq.0_int32) then
                     read(1) i_buf(1:n_qdo)
                     call mpi_send (i_buf(1:n_qdo), n_qdo, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, mpierr)
                  else if (mpi_rank.eq.i_mpi_task) then
                     call mpi_recv (drd_cnf(iw)%i_drd_tmv(1:n_qdo), n_qdo, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
                  endif
                  call mpi_barrier(MPI_COMM_WORLD,mpierr)
               enddo
            enddo
         endif
#endif
         call write_done(stdout,0,mpi_rank)
         if (mpi_rank.eq.0_int32) then
            close(1)
            deallocate(r_buf,i_buf)
         endif
      endif
   end subroutine load_drudonic_confs
   subroutine save_drudonic_confs( )
      character(50) :: svf_fle
      real(dp),       allocatable, dimension(:,:) :: r_buf
      integer(int32), allocatable, dimension(:)   :: i_buf
      integer(int32) :: icom, i_mpi_task, iw
      svf_fle = 'qmc.save/drd_conf.sav'
      if (mpi_rank.eq.0_int32) then
         open(unit=1,file=svf_fle,action='write',form='unformatted',status='unknown',access='sequential')
         do iw = 1_int32, n_wlk_max
            write(1) drd_cnf(iw)%r_drd(1:3,1:n_qdo)
            write(1) drd_cnf(iw)%n_drd_tmv
            write(1) drd_cnf(iw)%i_drd_tmv(1:n_qdo)
         enddo
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      if (n_mpi_tasks.gt.1_int32) then
         if (mpi_rank.eq.0_int32) then 
            allocate( r_buf(1:3,1:n_qdo) ) ; r_buf = 0.0_dp
            allocate( i_buf(1:n_qdo) ) ; i_buf = 0_int32
         endif
         do i_mpi_task = 1_int32, n_mpi_tasks-1
            do iw = 1_int32, n_wlk_max
               icom = n_wlk_max * (i_mpi_task-1) + iw
               if (mpi_rank.eq.i_mpi_task) then
                  call mpi_send (drd_cnf(iw)%r_drd(1:3,1:n_qdo), 3*n_qdo, MPI_DOUBLE_PRECISION, 0, icom, MPI_COMM_WORLD, mpierr)
               else if (mpi_rank.eq.0) then
                  call mpi_recv (r_buf, 3*n_qdo, MPI_DOUBLE_PRECISION, i_mpi_task, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
               endif
               call mpi_barrier(MPI_COMM_WORLD, mpierr)
               if (mpi_rank.eq.0_int32) write(1) r_buf(1:3,1:n_qdo)
               if (mpi_rank.eq.i_mpi_task) then
                  call mpi_send (drd_cnf(iw)%n_drd_tmv, 1, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, mpierr)
               else if (mpi_rank.eq.0) then
                  call mpi_recv (i_buf(1), 1, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
               endif
               call mpi_barrier(MPI_COMM_WORLD, mpierr)
               if (mpi_rank.eq.0_int32)  write(1) i_buf(1)
               if (mpi_rank.eq.i_mpi_task) then
                  call mpi_send (drd_cnf(iw)%i_drd_tmv(1:n_qdo), n_qdo, MPI_INTEGER, 0, icom, MPI_COMM_WORLD, mpierr)
               else if (mpi_rank.eq.0) then
                  call mpi_recv (i_buf, n_qdo, MPI_INTEGER, i_mpi_task, icom, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
               endif
               call mpi_barrier(MPI_COMM_WORLD, mpierr)
               if (mpi_rank.eq.0_int32)  write(1) i_buf(1:n_qdo)
            enddo
         enddo
      endif
#endif
      if (mpi_rank.eq.0_int32) then
         close(1)
         deallocate(r_buf,i_buf)
      endif
   end subroutine save_drudonic_confs
   subroutine cmp_qdo_dist
      integer(int32) :: i1, i2 ,i3
      allocate( d_qq(1:n_qdo*(n_qdo-1)/2) )
      if( n_qdo.ge.2_int32 ) then
         i3 = 0_int32
         do i2 = 1_int32, n_qdo
            do i1 = i2 + 1_int32, n_qdo
               i3 = i3 + 1_int32
               d_qq(i3)%v = r_qdo(:,i2) - r_qdo(:,i1)
               d_qq(i3)%m    = dsqrt(sum(d_qq(i3)%v(:)**2))
            enddo ; enddo 
      endif
   end subroutine cmp_qdo_dist
   subroutine cmp_qdonuc_dist()
      integer(int32) :: i1, i2
      allocate( d_qn(1:n_qdo,1:n_at) )
      do i1 = 1, n_qdo
         do i2 = 1, n_at
            d_qn(i1,i2)%v = r_qdo(:,i1) - r_at(:,i2)
            d_qn(i1,i2)%m = dsqrt( sum( d_qn(i1,i2)%v(:)**2 ) )
         enddo 
      enddo 
   end subroutine cmp_qdonuc_dist
   subroutine cmp_pchrg_dist()
      integer(int32) :: i1, i2 ,i3
      allocate( d_pchpch(1:n_pchrg*(n_pchrg-1)/2) )
      if( n_pchrg.ge.2_int32 ) then
         i3 = 0_int32
         do i2 = 1_int32, n_pchrg
            do i1 = i2 + 1_int32, n_pchrg
               i3 = i3 + 1_int32
               d_pchpch(i3)%v(:) = r_pchrg(:,i2) - r_pchrg(:,i1)
               d_pchpch(i3)%m    = dsqrt(sum(d_pchpch(i3)%v(:)**2))
            enddo ; enddo 
      endif
   end subroutine cmp_pchrg_dist
   subroutine cmp_qdopchrg_dist()
      integer(int32) :: i1, i2
      allocate( d_qpch(1:n_qdo,1:n_pchrg) )
      do i1 = 1, n_qdo
         do i2 = 1, n_pchrg
            d_qpch(i1,i2)%v = r_qdo(:,i1) - r_pchrg(:,i2)
            d_qpch(i1,i2)%m = dsqrt( sum( d_qpch(i1,i2)%v(:)**2 ) )
         enddo 
      enddo 
   end subroutine cmp_qdopchrg_dist
   subroutine cmp_nucpchrg_dist()
      integer(int32) :: i1, i2
      allocate( d_npch(1:n_at,1:n_pchrg) )
      do i1 = 1, n_at
         do i2 = 1, n_pchrg
            d_npch(i1,i2)%v = r_at(:,i1) - r_pchrg(:,i2)
            d_npch(i1,i2)%m = dsqrt( sum( d_npch(i1,i2)%v(:)**2 ) )
         enddo
      enddo 
   end subroutine cmp_nucpchrg_dist
   subroutine read_qcoord_fle( )
      integer(int32) :: ierr
      character(4)   :: units
      integer(int32) :: i1
      namelist /qdosys/ n_qdo, n_pchrg, el_sigma, ext_E_field, units, qdoham
      units    = 'bohr'
      n_qdo    = 0_int32    
      n_pchrg  = 0_int32     
      el_sigma = 0.0_dp     
      qdoham = 'cou'
      ierr = 0_int32
      if (trim(fle_coord_qdo).ne.'none') then
         call write_simple_line(stdout,0,mpi_rank,2,"l","_______________________")
         call write_variable_line(stdout,0,mpi_rank,2,"Reading QDO system file", trim(fle_coord_qdo))
         qdosys_prs = .true.
         if(mpi_rank.eq.0_int32) then
            open(unit=1,file=fle_coord_qdo,action='read',form='formatted')
            read(1,nml=qdosys,IOSTAT=ierr)
            call check_namelist_mandatory( ierr, "&qdosys" )
            read(1,*)  ! comment line
            allocate( qdos(1:n_qdo), r_qdo(1:3,1:n_qdo) )
            do i1 = 1_int32, n_qdo
               read(1,*) qdos(i1)%qdo_name, qdos(i1)%qdo_q, qdos(i1)%qdo_omega, qdos(i1)%qdo_m, qdos(i1)%qdo_sigma_c, qdos(i1)%qdo_sigma_d, r_qdo(1:3,i1)
            enddo
            if (units.eq.'angs') r_qdo(:,:) = r_qdo(:,:) * ang_to_au
            if (n_pchrg.gt.0_int32) read(1,*)
            allocate( pchrgs(1:n_pchrg), r_pchrg(1:3,1:n_pchrg) )
            do i1 = 1_int32, n_pchrg
               read(1,*) pchrgs(i1)%pchrg_name, pchrgs(i1)%pchrg_q, pchrgs(i1)%pchrg_sigma , r_pchrg(1:3,i1)
            enddo
            if (units.eq.'angs') r_pchrg(:,:) = r_pchrg(:,:) * ang_to_au
            close(1)
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         call mpi_bcast(ierr,        1, MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(n_qdo,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(n_pchrg,     1, MPI_INTEGER,          0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(el_sigma,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(ext_E_field, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         call mpi_bcast(qdoham,      3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
         if(mpi_rank.ne.0_int32) allocate( qdos(1:n_qdo),     r_qdo(1:3,1:n_qdo    ) )
         do i1 = 1_int32, n_qdo
            call mpi_bcast(qdos(i1)%qdo_name,    4, MPI_CHARACTER,        0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast(qdos(i1)%qdo_q,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast(qdos(i1)%qdo_omega,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast(qdos(i1)%qdo_m,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast(qdos(i1)%qdo_sigma_c, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast(qdos(i1)%qdo_sigma_d, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         enddo
         call mpi_bcast(r_qdo, 3 * n_qdo, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         if(mpi_rank.ne.0_int32) allocate( pchrgs(1:n_pchrg), r_pchrg(1:3,1:n_pchrg) )
         do i1 = 1_int32, n_pchrg
            call mpi_bcast(pchrgs(i1)%pchrg_name,  3, MPI_CHARACTER,        0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast(pchrgs(i1)%pchrg_q,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
            call mpi_bcast(pchrgs(i1)%pchrg_sigma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
         enddo
         call mpi_bcast(r_pchrg, 3 * n_pchrg, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
#endif
      else
         qdosys_prs = .false.
      endif
      if ( ierr .ne. 0_int32 ) then
         call fnlz_ompmpi_env()
         stop
      endif
   end subroutine read_qcoord_fle
end module qdo_system_m
