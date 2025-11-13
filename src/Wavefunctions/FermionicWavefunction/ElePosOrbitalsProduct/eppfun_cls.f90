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
module eppfun_cls
   use fortran_kinds_v, only: dp, int32
   use fermionic_wavefunction_v, only: spin
   use molecular_system_v, only: n_el, n_po
   use eppcfs_mod, only: eppcfs_t
   use quantum_monte_carlo_v, only: grad_updat, n_wlk_max
   use feporb_mod, only: n_eporbs, n_par_eporbs, eporbs, eporbs_s
   implicit none
   type, public :: eppfun_t
      real(dp), allocatable, dimension(:)       :: gem, inv_gem, g
      integer(int32)                            :: n_max_dim
      real(dp), allocatable, dimension(:,:) :: vec_S
      real(dp), allocatable, dimension(:,:) :: vec_S_dif
      real(dp), allocatable, dimension(:,:) :: DD2_ln_gem
      real(dp), allocatable, dimension(:,:) :: DD2_ln_gem_new
      real(dp), allocatable, dimension(:,:) :: DD2_ln_gem_t
      real(dp), allocatable, dimension(:,:) :: DD2_ln_gem_new_t
   contains
      procedure :: ini   => ini_epo_fun
      procedure :: cmp   => cmp_epo_fun
      procedure :: ratio => ratio_epo_fun
      procedure :: upd   => upd_epo_fun
      procedure :: cmp_D => cmp_DD2_epo_fun
      procedure :: new_D => new_DD2_epo_fun
      procedure :: upd_D => upd_DD2_epo_fun
      procedure :: cmp_dl => cmp_dl_epo_fun
      procedure :: cmp_da => cmp_da_epo_fun
   end type eppfun_t
   private :: ini_epo_fun, cmp_epo_fun, ratio_epo_fun, upd_epo_fun, &
   & cmp_DD2_epo_fun, upd_DD2_epo_fun, new_DD2_epo_fun, cmp_dl_epo_fun, &
   & cmp_da_epo_fun
contains
   subroutine ini_epo_fun( obj )
      class(eppfun_t), intent(inout) :: obj
      if (n_el.gt.n_po) then
         obj%n_max_dim = n_el
      else
         obj%n_max_dim = n_po
      endif
      allocate( obj%gem(1:n_wlk_max) )     ; obj%gem = 1.0_dp
      allocate( obj%g(1:n_wlk_max) )       ; obj%g = 1.0_dp
      allocate( obj%inv_gem(1:n_wlk_max) ) ; obj%inv_gem = 1.0_dp
      allocate( obj%vec_S(1:n_po*n_el,1:n_wlk_max) ) ; obj%vec_S = 0.0_dp
      allocate( obj%vec_S_dif(1:obj%n_max_dim,1:n_wlk_max) )    ; obj%vec_S_dif = 0.0_dp
      allocate( obj%DD2_ln_gem(1:4*(n_el+n_po),1:n_wlk_max) )  ; obj%DD2_ln_gem = 0.0_dp
      allocate( obj%DD2_ln_gem_t(1:4*n_el*n_po,1:n_wlk_max ) ) ; obj%DD2_ln_gem_t = 0.0_dp
      if ( grad_updat ) then
         allocate( obj%DD2_ln_gem_new_t(1:4*obj%n_max_dim,1:n_wlk_max) ) ; obj%DD2_ln_gem_new_t = 0.0_dp
         allocate( obj%DD2_ln_gem_new(1:4,1:n_wlk_max) ) ; obj%DD2_ln_gem_new = 0.0_dp
      endif
   end subroutine ini_epo_fun
   subroutine cmp_epo_fun( obj, iw, n_el, n_po, epo_c, mat_o )
      class(eppfun_t), intent(inout) :: obj
      integer(int32),  intent(in) :: iw
      integer(int32),                              intent(in) :: n_el
      integer(int32),                              intent(in) :: n_po
      type(eppcfs_t),                              intent(in) :: epo_c
      real(dp), dimension(n_eporbs,n_el*n_po), intent(in) :: mat_o
      integer(int32)                                          :: i1
      obj%gem(iw) = 1.0_dp
      do i1 = 1, n_el*n_po
         obj%vec_S(i1,iw) = dot_product(epo_c%mat_L(:,i1) , mat_o(:,i1))
         obj%gem(iw) = obj%gem(iw) * obj%vec_S(i1,iw)
      end do
      obj%inv_gem(iw) = 1.0_dp / obj%gem(iw)
   end subroutine cmp_epo_fun
   subroutine ratio_epo_fun( obj, iw, i_fe, n_el, n_po, epo_c, vec_o_dif, g )
      class(eppfun_t), intent(inout) :: obj
      integer(int32),  intent(in) :: iw
      integer(int32),                                  intent(in)  :: i_fe
      integer(int32),                                  intent(in)  :: n_el
      integer(int32),                                  intent(in)  :: n_po
      type(eppcfs_t),                                  intent(in)  :: epo_c
      real(dp), dimension(n_eporbs,obj%n_max_dim), intent(in)  :: vec_o_dif
      real(dp),                                        intent(out) :: g
      integer(int32)                                               :: i1, i2
      obj%g(iw) = 1.0_dp
      if ( i_fe .le. n_el ) then 
         i1 = n_po*(i_fe-1)
         do i2 = 1_int32, n_po
            i1 = i1 + 1_int32
            obj%vec_S_dif(i2,iw) = dot_product(epo_c%mat_L(:,i1) , vec_o_dif(:,i2) )
            obj%g(iw) = obj%g(iw) * ( 1.0_dp + obj%vec_S_dif(i2,iw) / obj%vec_S(i1,iw) )
         end do
      else 
         i1 = i_fe-n_el
         do i2 = 1_int32, n_el
            obj%vec_S_dif(i2,iw) = dot_product(epo_c%mat_L(:,i1) , vec_o_dif(:,i2) )
            obj%g(iw) = obj%g(iw) * ( 1.0_dp + obj%vec_S_dif(i2,iw) / obj%vec_S(i1,iw) )
            i1 = i1 + n_po
         end do
      endif
      g = obj%g(iw)
   end subroutine ratio_epo_fun
   subroutine new_DD2_epo_fun( obj, iw, i_fe, n_el, n_po, epo_c, vec_DD2_o_new )
      class(eppfun_t), intent(inout) :: obj
      integer(int32),  intent(in) :: iw
      integer(int32),                                  intent(in)  :: i_fe
      integer(int32),                                  intent(in)  :: n_el
      integer(int32),                                  intent(in)  :: n_po
      type(eppcfs_t),                                  intent(in)  :: epo_c
      real(dp), dimension(n_eporbs,4*obj%n_max_dim), intent(in) :: vec_DD2_o_new
      integer(int32)                                              :: i1, i2
      obj%DD2_ln_gem_new(:,iw) = 0.0_dp
      if ( i_fe .le. n_el ) then 
         i2 = n_po*(i_fe-1)
         do i1 = 1_int32, n_po
            i2 = i2 + 1_int32
            call dgemv('T',n_eporbs, 4_int32, 1.0_dp, vec_DD2_o_new(:,4*i1-3:4*i1), n_eporbs, &
            & epo_c%mat_L(:,i2), 1_int32, 0.0_dp, obj%DD2_ln_gem_new_t(4*i1-3:4*i1,iw), 1_int32 )
            obj%DD2_ln_gem_new_t(4*i1-3:4*i1,iw) = obj%DD2_ln_gem_new_t(4*i1-3:4*i1,iw) / (obj%vec_S(i2,iw) + obj%vec_S_dif(i1,iw))
            obj%DD2_ln_gem_new_t(4*i1,iw) =  obj%DD2_ln_gem_new_t(4*i1,iw) - sum( obj%DD2_ln_gem_new_t(4*i1-3:4*i1-1,iw)**2)
            obj%DD2_ln_gem_new(1:4,iw) = obj%DD2_ln_gem_new(1:4,iw) + obj%DD2_ln_gem_new_t(4*i1-3:4*i1,iw)
         enddo
      else 
         i2 = i_fe - n_el
         do i1 = 1_int32, n_el
            call dgemv('T',n_eporbs, 4_int32, 1.0_dp, vec_DD2_o_new(:,4*i1-3:4*i1), n_eporbs, &
            & epo_c%mat_L(:,i2), 1_int32, 0.0_dp, obj%DD2_ln_gem_new_t(4*i1-3:4*i1,iw), 1_int32 )
            obj%DD2_ln_gem_new_t(4*i1-3:4*i1,iw) = obj%DD2_ln_gem_new_t(4*i1-3:4*i1,iw) / (obj%vec_S(i2,iw) + obj%vec_S_dif(i1,iw))
            obj%DD2_ln_gem_new_t(4*i1,iw) =  obj%DD2_ln_gem_new_t(4*i1,iw) - sum( obj%DD2_ln_gem_new_t(4*i1-3:4*i1-1,iw)**2)
            obj%DD2_ln_gem_new(1:3,iw) = obj%DD2_ln_gem_new(1:3,iw) - obj%DD2_ln_gem_new_t(4*i1-3:4*i1-1,iw)
            obj%DD2_ln_gem_new(4,iw) = obj%DD2_ln_gem_new(4,iw) + obj%DD2_ln_gem_new_t(4*i1,iw)
            i2 = i2 + n_po
         enddo
      endif
   end subroutine new_DD2_epo_fun
   subroutine cmp_DD2_epo_fun( obj, iw, n_el, n_po, epo_c, DD2_o )
      class(eppfun_t), intent(inout) :: obj
      integer(int32),  intent(in) :: iw
      integer(int32),                                intent(in) :: n_el
      integer(int32),                                intent(in) :: n_po
      type(eppcfs_t),                                intent(in) :: epo_c
      real(dp), dimension(n_eporbs,4*n_el*n_po), intent(in) :: DD2_o
      integer(int32)                                            :: i1, i2, i3
      obj%DD2_ln_gem(:,iw) = 0.0_dp
      i3 = 0_int32
      do i1 = 1_int32, n_el; do i2 = n_el+1_int32, n_po+n_el
            i3 = i3 + 1_int32
            call dgemv('T',n_eporbs, 4_int32, 1.0_dp/ obj%vec_S(i3,iw), DD2_o(:,4*i3-3:4*i3), n_eporbs, &
            & epo_c%mat_L(:,i3), 1_int32, 0.0_dp, obj%DD2_ln_gem_t(4*i3-3:4*i3,iw), 1_int32 )
            obj%DD2_ln_gem_t(4*i3,iw) = obj%DD2_ln_gem_t(4*i3,iw) - sum(obj%DD2_ln_gem_t(4*i3-3:4*i3-1,iw)**2)
            obj%DD2_ln_gem(4*i1-3:4*i1,iw) = obj%DD2_ln_gem(4*i1-3:4*i1,iw) + obj%DD2_ln_gem_t(4*i3-3:4*i3,iw)
            obj%DD2_ln_gem(4*i2-3:4*i2-1,iw) = obj%DD2_ln_gem(4*i2-3:4*i2-1,iw) - obj%DD2_ln_gem_t(4*i3-3:4*i3-1,iw)
            obj%DD2_ln_gem(4*i2,iw) = obj%DD2_ln_gem(4*i2,iw) + obj%DD2_ln_gem_t(4*i3,iw)
         enddo ; enddo
   end subroutine cmp_DD2_epo_fun
   subroutine upd_epo_fun( obj, iw, i_fe, n_el, n_po )
      class(eppfun_t), intent(inout) :: obj
      integer(int32),  intent(in) :: iw
      integer(int32),                  intent(in)  :: i_fe
      integer(int32),                  intent(in)  :: n_el
      integer(int32),                  intent(in)  :: n_po
      integer(int32)                               :: i1, i2
      if ( i_fe .le. n_el ) then 
         i1 = n_po*(i_fe-1)
         do i2 = 1_int32, n_po
            i1 = i1 + 1_int32
            obj%vec_S(i1,iw) = obj%vec_S(i1,iw) + obj%vec_S_dif(i2,iw)
         end do
      else 
         i1 = i_fe-n_el
         do i2 = 1_int32, n_el
            obj%vec_S(i1,iw) = obj%vec_S(i1,iw) + obj%vec_S_dif(i2,iw)
            i1 = i1 + n_po
         end do
      endif
      obj%gem(iw) = obj%gem(iw) * obj%g(iw)
      obj%inv_gem(iw) = 1.0_dp / obj%gem(iw)
   end subroutine upd_epo_fun
   subroutine upd_DD2_epo_fun( obj, iw, i_fe, n_el, n_po, epo_c, DD2_o )
      class(eppfun_t), intent(inout) :: obj
      integer(int32),  intent(in) :: iw
      integer(int32),                                intent(in) :: i_fe
      integer(int32),                                intent(in) :: n_el
      integer(int32),                                intent(in) :: n_po
      type(eppcfs_t),                                intent(in) :: epo_c
      real(dp), dimension(n_eporbs,4*n_el*n_po), intent(in) :: DD2_o
      integer(int32)                                            :: i1, i2, i3
      obj%DD2_ln_gem(:,iw) = 0.0_dp
      if ( i_fe .le. n_el ) then 
         if ( grad_updat) then
            i1 = 4*n_po*(i_fe-1)
            obj%DD2_ln_gem_t(i1+1:i1+4*n_po,iw) = obj%DD2_ln_gem_new_t(1:4*n_po,iw)
         else
            i2 = n_po*(i_fe-1)
            do i1 = 1_int32, n_po
               i2 = i2 + 1_int32
               call dgemv('T',n_eporbs, 4_int32, 1.0_dp, DD2_o(:,4*i2-3:4*i2), n_eporbs, &
               & epo_c%mat_L(:,i2), 1_int32, 0.0_dp, obj%DD2_ln_gem_t(4*i2-3:4*i2,iw), 1_int32 )
               obj%DD2_ln_gem_t(4*i2-3:4*i2,iw)  = obj%DD2_ln_gem_t(4*i2-3:4*i2,iw) / obj%vec_S(i2,iw)
               obj%DD2_ln_gem_t(4*i2,iw) = obj%DD2_ln_gem_t(4*i2,iw) - sum(obj%DD2_ln_gem_t(4*i2-3:4*i2-1,iw)**2)
            enddo
         endif
      else 
         if ( grad_updat) then
            i1 = i_fe - n_el
            do i2 = 1_int32, n_el
               obj%DD2_ln_gem_t(4*i1-3:4*i1,iw) = obj%DD2_ln_gem_new_t(4*i2-3:4*i2,iw)
               i1 = i1 + n_po
            enddo
         else
            i2 = i_fe - n_el
            do i1 = 1_int32, n_el
               call dgemv('T',n_eporbs, 4_int32, 1.0_dp, DD2_o(:,4*i2-3:4*i2), n_eporbs, &
               & epo_c%mat_L(:,i2), 1_int32, 0.0_dp, obj%DD2_ln_gem_t(4*i2-3:4*i2,iw), 1_int32 )
               obj%DD2_ln_gem_t(4*i2-3:4*i2,iw)  = obj%DD2_ln_gem_t(4*i2-3:4*i2,iw) / obj%vec_S(i2,iw)
               obj%DD2_ln_gem_t(4*i2,iw) = obj%DD2_ln_gem_t(4*i2,iw) - sum(obj%DD2_ln_gem_t(4*i2-3:4*i2-1,iw)**2)
               i2 = i2 + n_po
            enddo
         endif
      endif
      do i1 = 1_int32, n_el
         i3 = n_po*(i1-1)
         do i2 = n_el+1_int32, n_el + n_po
            i3 = i3 + 1_int32
            obj%DD2_ln_gem(4*i1-3:4*i1,iw) = obj%DD2_ln_gem(4*i1-3:4*i1,iw) + obj%DD2_ln_gem_t(4*i3-3:4*i3,iw)
            obj%DD2_ln_gem(4*i2-3:4*i2-1,iw) = obj%DD2_ln_gem(4*i2-3:4*i2-1,iw) - obj%DD2_ln_gem_t(4*i3-3:4*i3-1,iw)
            obj%DD2_ln_gem(4*i2,iw) = obj%DD2_ln_gem(4*i2,iw) + obj%DD2_ln_gem_t(4*i3,iw)
         enddo
      enddo
   end subroutine upd_DD2_epo_fun
   subroutine cmp_dl_epo_fun( obj, iw, n_el, n_po, mat_o, n_par_det, dl_ln_gem )
      class(eppfun_t), intent(inout) :: obj
      integer(int32),  intent(in) :: iw
      integer(int32),                              intent(in)    :: n_el
      integer(int32),                              intent(in)    :: n_po
      real(dp), dimension(n_eporbs,n_el*n_po), intent(in)    :: mat_o
      integer(int32),                              intent(in)    :: n_par_det
      real(dp),          dimension(n_par_det),   intent(inout) :: dl_ln_gem
      integer(int32) :: i1, ip
      ip = 0_int32
      do i1 = 1_int32, n_el * n_po
         dl_ln_gem(ip+1:ip+n_eporbs) = mat_o(:,i1) / obj%vec_S(i1,iw)
         ip = ip + n_eporbs
      enddo
   end subroutine cmp_dl_epo_fun
   subroutine cmp_da_epo_fun( obj, iw, n_el, n_po, epo_c, mat_da_o, da_ln_gem )
      class(eppfun_t), intent(inout) :: obj
      integer(int32),  intent(in) :: iw
      integer(int32),                                  intent(in)    :: n_el
      integer(int32),                                  intent(in)    :: n_po
      type(eppcfs_t),                                  intent(in)    :: epo_c
      real(dp), dimension(n_par_eporbs,n_el*n_po), intent(in)    :: mat_da_o
      real(dp), dimension(n_par_eporbs),             intent(inout) :: da_ln_gem
      integer(int32)                                                 :: pi, pf
      integer(int32)                                                 :: i1, i2
      da_ln_gem = 0.0_dp
      do i1 = 1_int32, n_el * n_po
         do i2 = 1_int32, eporbs_s%n_orbs
            if( eporbs(i2)%n_par.gt.0_int32 ) then
               pi = eporbs(i2)%i_par
               pf = pi + eporbs(i2)%n_par - 1_int32
               da_ln_gem(pi:pf) = da_ln_gem(pi:pf) + epo_c%mat_L(i2,i1) * mat_da_o(pi:pf,i1) / obj%vec_S(i1,iw)
            endif
         enddo  
      enddo
   end subroutine cmp_da_epo_fun
end module eppfun_cls
