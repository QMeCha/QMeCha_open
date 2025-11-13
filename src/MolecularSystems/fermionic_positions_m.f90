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
module fermionic_positions_m
   use fortran_kinds_v,  only: int32, dp
   use openmp_mpi_m
   use molecular_system_v,  only: dist_t, n_fe, n_el, n_el_u, n_po, n_po_u, n_at
   implicit none
   integer(int32), save, public, allocatable, dimension(:,:) :: frmdst
contains
   subroutine read_frmdst( )
      logical        :: file_exists
      character(100) :: svf_fle
      integer(int32) :: i1
      svf_fle = 'wvfn.save/' // 'frmdst.sav'
      inquire(FILE=svf_fle, EXIST=file_exists)
      if (n_po.gt.0_int32) then
         allocate( frmdst(1:4,1:n_at) ) ; frmdst = 0_int32
      else
         allocate( frmdst(1:2,1:n_at) ) ; frmdst = 0_int32
      endif
      if ( file_exists ) then
         if ( mpi_rank.eq.0_int32 ) then
            open(unit=1,file=svf_fle,action='read',form='formatted',status='unknown')
            do i1 = 1_int32, n_at
               if (n_po.gt.0_int32) then
                  read(1,*) frmdst(1:4,i1)
               else
                  read(1,*) frmdst(1:2,i1)
               endif
            enddo
            close(1)
         endif
#if defined _MPI || defined _MPIh || defined _MPI08
         if (n_po.gt.0_int32) then
            call mpi_bcast(frmdst,        4*n_at, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         else
            call mpi_bcast(frmdst,        2*n_at, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
         endif
#endif
      endif
   end subroutine read_frmdst
   subroutine deallocate_frmdst( )
      deallocate( frmdst )
   end subroutine deallocate_frmdst
   subroutine save_frmdst( d_fn )
      type(dist_t), dimension(n_at,1:n_fe), intent(in)  :: d_fn
      integer(int32) :: i_at_min
      character(100) :: svf_fle
      integer(int32) :: i1
      svf_fle = 'wvfn.save/' // 'frmdst.sav'
      if ( mpi_rank.eq.0_int32 ) then
         if (n_po.gt.0_int32) then
            allocate( frmdst(1:4,1:n_at) ) ; frmdst = 0_int32
         else
            allocate( frmdst(1:2,1:n_at) ) ; frmdst = 0_int32
         endif
         do i1 = 1_int32, n_el
            call nearest_atom( d_fn(:,i1), i_at_min )
            if (i_at_min.gt.0_int32) then
               if (i1.le.n_el_u) then
                  frmdst(1,i_at_min) = frmdst(1,i_at_min) + 1_int32
               else
                  frmdst(2,i_at_min) = frmdst(2,i_at_min) + 1_int32
               endif
            endif
         enddo
         if (n_po.gt.0_int32) then
            do i1 = 1_int32, n_po
               call nearest_atom( d_fn(:,i1+n_el), i_at_min )
               if (i_at_min.gt.0_int32) then
                  if (i1.le.n_po_u) then
                     frmdst(3,i_at_min) = frmdst(3,i_at_min) + 1_int32
                  else
                     frmdst(4,i_at_min) = frmdst(4,i_at_min) + 1_int32
                  endif
               endif
            enddo
         endif
         open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
         do i1 = 1_int32, n_at
            if (n_po.gt.0_int32) then
               write(1,*) frmdst(1:4,i1)
            else
               write(1,*) frmdst(1:2,i1)
            endif
         enddo
         close(1)
         deallocate( frmdst )
      endif
   end subroutine save_frmdst
end module fermionic_positions_m
