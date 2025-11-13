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
module orbbssq_cls
   use fortran_kinds_v, only: int32, dp
   use orbspq_cls, only: orbspq_t
   implicit none
   type, public :: qsysorb_t
      integer(int32)                            :: n_orbs
      type(orbspq_t), allocatable, dimension(:) :: orbs
      logical, pointer                          :: oc_opt
      logical, pointer                          :: oe_opt
      integer(int32)                            :: n_par_orbs
      integer(int32)                            :: n_par_orbs_s
   contains
      procedure :: ini     => ini_sys_orbs
      procedure :: upd     => upd_sys_orbs
   end type qsysorb_t
   private :: ini_sys_orbs, cnt_orbs_par, upd_sys_orbs
contains
   subroutine ini_sys_orbs( obj, bsst, qc_opt, qe_opt )
      use basisset_v, only: n_qbs, bssset_t
      use qdo_system_v, only: n_qdo
      class(qsysorb_t), intent(inout) :: obj
      type(bssset_t), dimension(n_qbs),         intent(in) :: bsst
      logical,                             target, intent(in) :: qc_opt, qe_opt
      integer(int32) :: io, iq
      obj%oc_opt => qc_opt
      obj%oe_opt => qe_opt
      if ( obj%oe_opt ) obj%oc_opt = .true.
      allocate( obj%orbs(1:n_qdo) )
      io = 1_int32
      do iq = 1_int32, n_qdo
         call obj%orbs(iq)%ini( iq, io, bsst )
      enddo
      obj%n_orbs = sum( obj%orbs(:)%n_orbs )
      obj%n_par_orbs   = 0_int32
      obj%n_par_orbs_s = 0_int32
      if ( obj%n_orbs.eq.0_int32 ) then
         obj%oe_opt = .false.
         obj%oc_opt = .false.
      endif
      if ( obj%oc_opt ) then
         call cnt_orbs_par( obj )
      endif
   end subroutine ini_sys_orbs
   subroutine cnt_orbs_par( obj )
      use qdo_system_v, only: n_qdo
      class(qsysorb_t), intent(inout) :: obj
      integer(int32) :: i1, i2, i3
      obj%n_par_orbs = 0_int32
      obj%n_par_orbs_s = 0_int32
      if ( obj%oc_opt ) then
         i3 = 1_int32
         do i1 = 1_int32, n_qdo ; do i2 = 1_int32, obj%orbs(i1)%n_orbs
               obj%orbs(i1)%orb(i2)%i_par = i3
               if ( obj%oe_opt ) then
                  if( obj%orbs(i1)%orb(i2)%n.gt.1_int32 ) then
                     obj%orbs(i1)%orb(i2)%n_par = 2 * obj%orbs(i1)%orb(i2)%n
                  else
                     obj%orbs(i1)%orb(i2)%n_par = 1_int32
                  endif
               else
                  if ( obj%orbs(i1)%orb(i2)%n.gt.1_int32 ) then
                     obj%orbs(i1)%orb(i2)%n_par = obj%orbs(i1)%orb(i2)%n
                  else
                     obj%orbs(i1)%orb(i2)%n_par = 0_int32
                  endif
               endif
               i3 = i3 + obj%orbs(i1)%orb(i2)%n_par
               obj%n_par_orbs = obj%n_par_orbs + obj%orbs(i1)%orb(i2)%n_par
            enddo ; enddo 
      endif
      if( obj%n_par_orbs.eq.0_int32 ) then
         obj%oc_opt = .false.
         obj%oe_opt = .false.
      else
         obj%n_par_orbs_s = obj%n_par_orbs
      endif
   end subroutine cnt_orbs_par
   subroutine upd_sys_orbs( obj, da, z_cut )
      use qdo_system_v, only: n_qdo
      class(qsysorb_t), intent(inout) :: obj
      real(dp), dimension(obj%n_par_orbs), intent(in) :: da
      real(dp)                               , intent(in) :: z_cut
      real(dp)                :: z_new
      integer(int32) :: i1, i2, i3
      integer(int32) :: ip
      ip = 1_int32
      do i1 = 1_int32, n_qdo ; do i2 = 1_int32, obj%orbs(i1)%n_orbs
            if ( obj%orbs(i1)%orb(i2)%n_par .gt. 0_int32 ) then
               if ( obj%orbs(i1)%orb(i2)%n .eq. 1_int32 ) then
                  if( obj%oe_opt ) then
                     z_new = obj%orbs(i1)%orb(i2)%z(1) + da(ip)
                     if ( z_new .gt. z_cut ) then
                        obj%orbs(i1)%orb(i2)%z(1) = z_new
                     else
                        obj%orbs(i1)%orb(i2)%z(1) = z_cut
                     endif
                     ip = ip + 1_int32
                  endif
               else
                  if( obj%oe_opt ) then
                     obj%orbs(i1)%orb(i2)%c(1:obj%orbs(i1)%orb(i2)%n) &
                     & = obj%orbs(i1)%orb(i2)%c(1:obj%orbs(i1)%orb(i2)%n) &
                     & + da(ip:ip+obj%orbs(i1)%orb(i2)%n-1)
                     ip = ip + obj%orbs(i1)%orb(i2)%n
                     do i3 = 1_int32, obj%orbs(i1)%orb(i2)%n
                        z_new = obj%orbs(i1)%orb(i2)%z(i3) + da(ip)
                        if ( z_new .gt. z_cut ) then
                           obj%orbs(i1)%orb(i2)%z(i3) = z_new
                        else
                           obj%orbs(i1)%orb(i2)%z(i3) = z_cut
                        endif
                        ip = ip + 1_int32
                     enddo
                  else
                     obj%orbs(i1)%orb(i2)%c(1:obj%orbs(i1)%orb(i2)%n) &
                     & = obj%orbs(i1)%orb(i2)%c(1:obj%orbs(i1)%orb(i2)%n) &
                     & + da(ip:ip+obj%orbs(i1)%orb(i2)%n-1)
                     ip = ip + obj%orbs(i1)%orb(i2)%n
                  endif
               endif
            endif
         enddo ; enddo
   end subroutine upd_sys_orbs
end module orbbssq_cls
