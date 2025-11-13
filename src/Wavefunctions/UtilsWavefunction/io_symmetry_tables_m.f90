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
module io_symmetry_tables_m
   use fortran_kinds_v, only: dp, int32
   use geminal_symmetries_m, only: lmbd_ele_t, lmbd_sym_t, conv_indx
   implicit none
   public :: wrte_symGtb, wrte_symltb, read_symGtb, read_symltb
contains
   subroutine wrte_symGtb( ifile, n_nz, n_sym, nz, sym )
      integer(int32),                             intent(in) :: ifile
      integer(int32),                             intent(in) :: n_nz, n_sym
      type(lmbd_ele_t), dimension(n_nz),        intent(in) :: nz
      type(lmbd_sym_t), dimension(n_sym), target, intent(in) :: sym
      integer(int32)          :: i1, i2
      integer(int32), pointer :: o1
      write(ifile,'("# Table of non zero elements of pairing matrix")')
      write(ifile,'(I6)') n_nz
      if (n_nz.gt.0_int32) then
         do i1 = 1_int32, n_nz
            write(ifile,'(2I6)') nz(i1)%o1, nz(i1)%o2
         enddo
      endif
      write(ifile,'("# Symmetry table of pairing matrix")')
      write(ifile,'(I6)') n_sym
      if (n_sym.gt.0_int32) then
         do i1 = 1_int32, n_sym
            o1 => sym(i1)%n_c
            write(ifile,'(1000I6)') o1,( sign(1,sym(i1)%c(i2)) * nz(abs(sym(i1)%c(i2)))%o1, &
            & nz(abs(sym(i1)%c(i2)))%o2, i2 = 1_int32, o1 )
         enddo
      endif
   end subroutine wrte_symGtb
   subroutine wrte_symltb( ifile, n_nz, n_sym, nz, sym )
      integer(int32),                             intent(in) :: ifile
      integer(int32),                             intent(in) :: n_nz, n_sym
      type(lmbd_ele_t), dimension(n_nz),        intent(in) :: nz
      type(lmbd_sym_t), dimension(n_sym), target, intent(in) :: sym
      integer(int32)          :: i1, i2
      integer(int32), pointer :: o1
      write(ifile,'("# Table of non zero elements of the unpaired orbital")')
      write(ifile,'(I6)') n_nz
      if (n_nz.gt.0_int32) then
         do i1 = 1_int32, n_nz
            write(ifile,'(I6)') nz(i1)%o1
         enddo
      endif
      write(ifile,'("# Symmetry table of the unpaired orbital")')
      write(ifile,'(I6)') n_sym
      if (n_sym.gt.0_int32) then
         do i1 = 1_int32, n_sym
            o1 => sym(i1)%n_c
            write(ifile,'(1000I6)') o1,( sign(1,sym(i1)%c(i2)) * nz(abs(sym(i1)%c(i2)))%o1, i2 = 1_int32, o1)
         enddo
      endif
   end subroutine wrte_symltb
   subroutine read_symGtb( ifile, n_nz, n_sym, n_orbs, nz, sym )
      integer(int32),                               intent(in)    :: ifile
      integer(int32),                               intent(inout) :: n_nz, n_sym
      integer(int32),                               intent(in)    :: n_orbs
      type(lmbd_ele_t), allocatable, dimension(:),  intent(inout) :: nz
      type(lmbd_sym_t), allocatable, dimension(:),  intent(inout) :: sym
      integer(int32)          :: i1, i2
      integer(int32)          :: o1
      integer(int32) :: n_c
      integer(int32), allocatable, dimension(:) :: c
      n_nz  = 0_int32
      n_sym = 0_int32
      read(ifile,*)
      read(ifile,*) n_nz
      if (n_nz.gt.0_int32) then
         allocate( nz(1:n_nz) )
         do i1 = 1_int32, n_nz
            read(ifile,*) nz(i1)%o1, nz(i1)%o2
            nz(i1)%io = nz(i1)%o1 + n_orbs * ( nz(i1)%o2 -1 )
         enddo
         allocate( c(1:2*n_nz) )
      endif
      read(ifile,*) 
      read(ifile,*) n_sym
      if (n_sym.gt.0_int32) then
         allocate(sym(1:n_sym))
         do i1 = 1_int32, n_sym
            read(ifile,*) n_c, ( c(2*i2-1:2*i2), i2 = 1_int32, n_c )
            sym(i1)%is_indp = .true.
            sym(i1)%n_c     = n_c
            allocate( sym(i1)%c(1:n_c) )
            do i2 = 1_int32, n_c
               o1 = abs(c(2*i2-1)) + n_orbs * ( c(2*i2) - 1 )
               call conv_indx( n_nz, nz, o1 )
               sym(i1)%c(i2) = sign(1,c(2*i2-1)) * o1
            enddo ! i2 (n_c)
         enddo
         if (allocated(c)) deallocate(c)
      endif
   end subroutine read_symGtb
   subroutine read_symltb( n_nz, n_sym, nz, sym )
      integer(int32),                               intent(inout) :: n_nz, n_sym
      type(lmbd_ele_t), allocatable, dimension(:),  intent(inout) :: nz
      type(lmbd_sym_t), allocatable, dimension(:),  target, intent(inout) :: sym
      integer(int32)          :: i1
      integer(int32) :: n_c
      integer(int32), allocatable, dimension(:) :: c
      n_nz  = 0_int32
      n_sym = 0_int32
      read(1,*) 
      read(1,*) n_nz
      if (n_nz.gt.0_int32) then
         allocate( nz(1:n_nz) )
         do i1 = 1_int32, n_nz
            read(1,*) nz(i1)%o1
            nz(i1)%io = nz(i1)%o1
         enddo
         allocate( c(1:n_nz) )
      endif
      read(1,*) 
      read(1,*) n_sym
      if (n_sym.gt.0_int32) then
         allocate( sym(1:n_sym) )
         do i1 = 1_int32, n_sym
            read(1,*) n_c, c(1:n_c)
            sym(i1)%is_indp = .true.
            sym(i1)%n_c     = n_c
            allocate( sym(i1)%c(1:n_c) )
            sym(i1)%c(1:n_c) = c(1:n_c)
         enddo
         if (allocated(c)) deallocate(c)
      endif
   end subroutine read_symltb
end module io_symmetry_tables_m
