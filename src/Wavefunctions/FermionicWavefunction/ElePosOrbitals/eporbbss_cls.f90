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
module eporbbss_cls
   use fortran_kinds_v, only: int32, dp
   use atomic_orbital_c, only: orb_t
   use orbbss_cls, only: sym_tab_t
   use basisset_v, only: n_fbs, bssset_t
   implicit none
   type, public :: epsysorb_t
      integer(int32)                            :: n_orbs
      integer(int32)                            :: l_max
      type(orb_t), allocatable, dimension(:)    :: orb
      logical, pointer                          :: oc_opt
      logical, pointer                          :: oe_opt
      integer(int32)                            :: n_par_orbs
      integer(int32)                            :: n_par_orbs_s
   contains
      procedure :: ini     => ini_sys_orbs
      procedure :: upd     => upd_sys_orbs
      procedure :: sym     => apl_sym_orbs_vec
   end type epsysorb_t
   private :: ini_sys_orbs, cnt_orbs_par, upd_sys_orbs
contains
   subroutine ini_sys_orbs( obj, bsst, dc_opt, de_opt )
      class(epsysorb_t), intent(inout) :: obj
      type(bssset_t), dimension(n_fbs),         intent(in) :: bsst
      logical,                           target, intent(in) :: dc_opt, de_opt
      integer(int32) :: i1, i2, i3, i_bs
      integer(int32) :: n, l
      obj%n_orbs = 0_int32 ; i_bs = 0_int32
      do i1 = 1_int32, n_fbs
         if(trim(bsst(i1)%atm_name).eq.'Ps') then
            obj%n_orbs = bsst(i1)%n_orbs
            i_bs   = i1
            exit
         endif
      enddo 
      obj%oc_opt => dc_opt
      obj%oe_opt => de_opt
      if ( obj%oe_opt ) obj%oc_opt = .true.
      if ( obj%n_orbs.gt.0_int32 ) then
         allocate( obj%orb(1:obj%n_orbs) )
         obj%l_max = 0_int32
         i3 = 1_int32
         do i1 = 1_int32, bsst(i_bs)%n_shlls
            l = bsst(i_bs)%shlls(i1)%l
            n = bsst(i_bs)%shlls(i1)%n
            if( obj%l_max.lt.l ) obj%l_max = l
            do i2 = l**2, l*(l+2)
               allocate( obj%orb(i3)%c(1:n), obj%orb(i3)%z(1:n) )
               obj%orb(i3)       = bsst(i_bs)%shlls(i1)
               obj%orb(i3)%l_z   = i2
               obj%orb(i3)%i_orb = i3
               i3 = i3 + 1_int32
            enddo 
         enddo 
         if (obj%oc_opt .or. obj%oe_opt) then
            call cnt_orbs_par( obj )
            obj%n_par_orbs_s = obj%n_par_orbs
         else
            obj%n_par_orbs   = 0_int32
            obj%n_par_orbs_s = 0_int32
         endif
      else
         obj%n_par_orbs   = 0_int32
         obj%n_par_orbs_s = 0_int32
      endif
   end subroutine ini_sys_orbs
   subroutine cnt_orbs_par( obj )
      class(epsysorb_t), intent(inout) :: obj
      integer(int32) :: i2, i3
      obj%n_par_orbs = 0_int32
      if ( obj%oc_opt ) then
         i3 = 1_int32
         do i2 = 1_int32, obj%n_orbs
            obj%orb(i2)%i_par = i3
            if ( obj%oe_opt ) then
               if( obj%orb(i2)%n.gt.1_int32 ) then
                  obj%orb(i2)%n_par = 2 * obj%orb(i2)%n
               else
                  obj%orb(i2)%n_par = 1_int32
               endif
            else
               if ( obj%orb(i2)%n.gt.1_int32 ) then
                  obj%orb(i2)%n_par = obj%orb(i2)%n
               else
                  obj%orb(i2)%n_par = 0_int32
               endif
            endif
            i3 = i3 + obj%orb(i2)%n_par
            obj%n_par_orbs = obj%n_par_orbs + obj%orb(i2)%n_par
         enddo 
      endif
   end subroutine cnt_orbs_par
   subroutine upd_sys_orbs( obj, da, z_cut )
      class(epsysorb_t), intent(inout) :: obj
      real(dp), dimension(obj%n_par_orbs), intent(in) :: da
      real(dp)    , intent(in)            :: z_cut
      real(dp)                :: z_new
      integer(int32) :: i2, i3
      integer(int32) :: ip
      ip = 1_int32
      do i2 = 1_int32, obj%n_orbs
         if ( obj%orb(i2)%n_par .gt. 0_int32 ) then
            if ( obj%orb(i2)%n .eq. 1_int32 ) then
               if( obj%oe_opt ) then
                  z_new = obj%orb(i2)%z(1) + da(ip)
                  if ( z_new .gt. z_cut ) then
                     obj%orb(i2)%z(1) = z_new
                  else
                     obj%orb(i2)%z(1) = z_cut
                  endif
                  ip = ip + 1_int32
               endif
            else
               if( obj%oe_opt ) then
                  obj%orb(i2)%c(1:obj%orb(i2)%n) &
                  & = obj%orb(i2)%c(1:obj%orb(i2)%n) &
                  & + da(ip:ip+obj%orb(i2)%n-1)
                  ip = ip + obj%orb(i2)%n
                  do i3 = 1_int32, obj%orb(i2)%n
                     z_new = obj%orb(i2)%z(i3) + da(ip)
                     if ( z_new .gt. z_cut ) then
                        obj%orb(i2)%z(i3) = z_new
                     else
                        obj%orb(i2)%z(i3) = z_cut
                     endif
                     ip = ip + 1_int32
                  enddo
               else
                  obj%orb(i2)%c(1:obj%orb(i2)%n) &
                  & = obj%orb(i2)%c(1:obj%orb(i2)%n) &
                  & + da(ip:ip+obj%orb(i2)%n-1)
                  ip = ip + obj%orb(i2)%n
               endif
            endif
         endif
      enddo
   end subroutine upd_sys_orbs
   subroutine apl_sym_orbs_vec( obj, da, da_s )
      class(epsysorb_t), intent(inout) :: obj
      real(dp), dimension(obj%n_par_orbs),   intent(in)    :: da
      real(dp), dimension(obj%n_par_orbs_s), intent(inout) :: da_s
      da_s = da
   end subroutine apl_sym_orbs_vec
end module eporbbss_cls
