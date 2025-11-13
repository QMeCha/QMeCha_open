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
module orbspq_cls
   use fortran_kinds_v, only: dp, int32
   use atomic_orbital_c, only: orb_t
   use qdo_system_v, only: qdos
   use basisset_v, only: n_qbs, bssset_t
   implicit none
   type, public :: orbspq_t
      integer(int32) :: n_orbs
      integer(int32) :: l_max
      type(orb_t), allocatable, dimension(:) :: orb
   contains
      procedure :: ini => ini_orbs
   end type orbspq_t
contains
   subroutine ini_orbs( obj, i_qdo, i_orbs, bsst )
      class(orbspq_t), intent(inout) :: obj
      integer(int32),                    intent(in)    :: i_qdo
      integer(int32),                    intent(inout) :: i_orbs
      type(bssset_t), dimension(n_qbs), intent(in)    :: bsst
      integer(int32) :: i1, i2, i3, l, n
      obj%n_orbs = 0_int32
      do i1 = 1_int32, n_qbs
         if(trim(bsst(i1)%atm_name).eq.trim(qdos(i_qdo)%qdo_name)) then
            obj%n_orbs = bsst(i1)%n_orbs
            qdos(i_qdo)%i_bsQ   = i1
            allocate( obj%orb(1:obj%n_orbs) )
         endif
      enddo
      obj%l_max = 0_int32
      i3 = 1_int32
      do i1 = 1_int32, bsst(qdos(i_qdo)%i_bsQ)%n_shlls
         l = bsst(qdos(i_qdo)%i_bsQ)%shlls(i1)%l
         n = bsst(qdos(i_qdo)%i_bsQ)%shlls(i1)%n
         if( obj%l_max.lt.l ) obj%l_max = l
         do i2 = l**2, l*(l+2)
            allocate( obj%orb(i3)%c(1:n), obj%orb(i3)%z(1:n) )
            obj%orb(i3)       = bsst(qdos(i_qdo)%i_bsQ)%shlls(i1)
            obj%orb(i3)%l_z   = i2
            obj%orb(i3)%i_orb = i_orbs
            i_orbs = i_orbs + 1_int32
            i3 = i3 + 1_int32
         enddo
      enddo
   end subroutine ini_orbs
end module orbspq_cls
