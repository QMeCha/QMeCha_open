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
module orbmat_cls
   use fortran_kinds_v,    only: dp, int32
   use cubic_harmonics_m,  only: cubhrm_t
   use molecular_system_v,         only: dist_t, n_fe, n_at
   use fermionic_config_c, only: frm_cnf
   use orbbss_cls,         only: orbspa_t, sysorb_t
   use cubic_harmonics_m,  only: y
   use quantum_monte_carlo_v,         only: grad_updat, n_wlk_max
   implicit none
   type, public :: prim_po_t
      real(dp), allocatable, dimension(:) :: p
   end type prim_po_t
   type, public :: matorb_t
      integer(int32),                 pointer :: i_fe
      type(dist_t),  dimension(:,:),  pointer :: d_fn
      type(dist_t),  dimension(:),    pointer :: d_fn_new
      type(cubhrm_t),                 pointer :: y
      type(prim_po_t), allocatable, dimension(:) :: prim_po
      real(dp), allocatable, dimension(:,:) :: o
      real(dp), allocatable, dimension(:)   :: o_new
      real(dp), allocatable, dimension(:)   :: o_dif
      real(dp), allocatable, dimension(:,:) :: DD2_o
      real(dp), allocatable, dimension(:,:) :: DD2_o_new
      real(dp), allocatable, dimension(:,:) :: da_o
      real(dp), allocatable, dimension(:)   :: da_o_diff
   contains
      procedure :: ini   => ini_orbs_mat
      procedure :: cmp   => cmp_orbs_mat
      procedure :: new   => new_orbs_vec
      procedure :: new_D => new_orbD_vec
      procedure :: upd   => upd_orbs_mat
   end type matorb_t
contains
   subroutine ini_orbs_mat( obj, iw, orbs )
      class(matorb_t), intent(inout) :: obj
      type(sysorb_t),  intent(in)    :: orbs
      integer(int32),  intent(in)    :: iw
      integer(int32) :: i1, i2, i3
      obj%y        => y(iw)
      obj%d_fn     => frm_cnf(iw)%d_fn
      obj%d_fn_new => frm_cnf(iw)%d_fn_new
      obj%i_fe     => frm_cnf(iw)%i_fe
      allocate( obj%prim_po(1:orbs%n_orbs) )
      i3 = 1_int32
      do i1 = 1_int32, n_at ; do i2 = 1_int32, orbs%orbs(i1)%n_orbs
            allocate( obj%prim_po(i3)%p(1:orbs%orbs(i1)%orb(i2)%n) ) ; obj%prim_po(i3)%p = 0.0_dp
            i3 = i3 + 1_int32
         enddo ; enddo 
      allocate( obj%o(1:orbs%n_orbs,1:n_fe) )          ; obj%o = 0.0_dp
      allocate( obj%o_new(1:orbs%n_orbs) )             ; obj%o_new = 0.0_dp
      allocate( obj%o_dif(1:orbs%n_orbs) )             ; obj%o_dif = 0.0_dp
      allocate( obj%DD2_o(1:orbs%n_orbs,1:4*n_fe) )    ; obj%DD2_o = 0.0_dp
      allocate( obj%DD2_o_new(1:orbs%n_orbs,1:4) )     ; obj%DD2_o_new = 0.0_dp
      if ( orbs%n_par_orbs.gt.0_int32 ) then
         allocate( obj%da_o(1:orbs%n_par_orbs,1:n_fe) ) ; obj%da_o = 0.0_dp
         allocate( obj%da_o_diff(1:orbs%n_par_orbs) )   ; obj%da_o_diff = 0.0_dp
      endif
   end subroutine ini_orbs_mat
   subroutine cmp_orbs_mat( obj, orbs )
      class(matorb_t), intent(inout) :: obj
      type(sysorb_t),  intent(inout) :: orbs
      integer(int32) :: i1, i2, i3, i4
      do i1 = 1_int32, n_fe
         i4 = 1_int32
         do i2 = 1_int32, n_at
            call obj%y%cmp( obj%d_fn(i2,i1), orbs%orbs(i2)%l_max )
            call obj%y%cmp_D( obj%d_fn(i2,i1), orbs%orbs(i2)%l_max )
            do i3 = 1_int32, orbs%orbs(i2)%n_orbs
               call orbs%orbs(i2)%orb(i3)%cmp( obj%d_fn(i2,i1)%m, obj%y, obj%prim_po(i4)%p, obj%o(i4,i1) )
               if ( orbs%oc_opt ) then
                  call orbs%orbs(i2)%orb(i3)%cmp_der( obj%d_fn(i2,i1), obj%y, orbs%n_orbs, orbs%n_par_orbs, &
                  & orbs%oe_opt, obj%prim_po(i4)%p, obj%o(i4,i1), &
                  & obj%DD2_o(:,4*i1-3:4*i1), obj%da_o(:,i1) )
               else
                  call orbs%orbs(i2)%orb(i3)%cmp_DD2( obj%d_fn(i2,i1), obj%y, orbs%n_orbs, obj%prim_po(i4)%p,&
                  & obj%o(i4,i1), obj%DD2_o(:,4*i1-3:4*i1) )
               endif
               i4 = i4 + 1_int32
            enddo ; enddo ; enddo 
   end subroutine cmp_orbs_mat
   subroutine new_orbs_vec( obj, orbs )
      class(matorb_t), intent(inout) :: obj
      type(sysorb_t),  intent(inout)    :: orbs
      integer(int32) :: i1, i2, i3
      i3 = 1_int32
      do i1 = 1_int32, n_at
         call obj%y%cmp( obj%d_fn_new(i1), orbs%orbs(i1)%l_max )
         do i2 = 1_int32, orbs%orbs(i1)%n_orbs
            call orbs%orbs(i1)%orb(i2)%cmp( obj%d_fn_new(i1)%m, obj%y, obj%prim_po(i3)%p, obj%o_new(i3) )
            i3 = i3 + 1_int32
         enddo ; enddo  
      obj%o_dif(:) = obj%o_new(:) - obj%o(:,obj%i_fe)
   end subroutine new_orbs_vec
   subroutine new_orbD_vec( obj, orbs )
      class(matorb_t), intent(inout) :: obj
      type(sysorb_t),  intent(inout)    :: orbs
      integer(int32) :: i1, i2, i3
      i3 = 1_int32
      do i1 = 1_int32, n_at
         call obj%y%cmp( obj%d_fn_new(i1), orbs%orbs(i1)%l_max )
         call obj%y%cmp_D( obj%d_fn_new(i1), orbs%orbs(i1)%l_max )
         do i2 = 1_int32, orbs%orbs(i1)%n_orbs
            call orbs%orbs(i1)%orb(i2)%cmp( obj%d_fn_new(i1)%m, obj%y, obj%prim_po(i3)%p, obj%o_new(i3) )
            call orbs%orbs(i1)%orb(i2)%cmp_DD2( obj%d_fn_new(i1), obj%y, orbs%n_orbs, obj%prim_po(i3)%p,&
            & obj%o_new(i3), obj%DD2_o_new(:,1:4) )
            i3 = i3 + 1_int32
         enddo ; enddo  
      obj%o_dif(:) = obj%o_new(:) - obj%o(:,obj%i_fe)
   end subroutine new_orbD_vec
   subroutine upd_orbs_mat( obj, orbs)
      class(matorb_t), intent(inout) :: obj
      type(sysorb_t),  intent(inout) :: orbs
      integer(int32)          :: i1, i2, i3
      obj%o(:,obj%i_fe) = obj%o(:,obj%i_fe) + obj%o_dif(:)
      if ( orbs%oc_opt ) obj%da_o_diff(:) = - obj%da_o(:,obj%i_fe)
      if ( grad_updat ) then
         obj%DD2_o(:,4*obj%i_fe-3:4*obj%i_fe) = obj%DD2_o_new(:,1:4)
         if ( orbs%oc_opt ) then
            i3 = 1_int32
            do i1 = 1_int32, n_at
               do i2 = 1_int32, orbs%orbs(i1)%n_orbs
                  call orbs%orbs(i1)%orb(i2)%cmp_dc( obj%d_fn(i1,obj%i_fe), orbs%n_par_orbs, &
                  & orbs%oe_opt, obj%prim_po(i3)%p, obj%da_o(:,obj%i_fe) )
                  i3 = i3 + 1_int32
               enddo ; enddo 
         endif
      else
         i3 = 1_int32
         do i1 = 1_int32, n_at
            call obj%y%cmp( obj%d_fn(i1,obj%i_fe), orbs%orbs(i1)%l_max )
            call obj%y%cmp_D( obj%d_fn(i1,obj%i_fe), orbs%orbs(i1)%l_max )
            do i2 = 1_int32, orbs%orbs(i1)%n_orbs
               if ( orbs%oc_opt ) then
                  call orbs%orbs(i1)%orb(i2)%cmp_der( obj%d_fn(i1,obj%i_fe), obj%y, orbs%n_orbs, orbs%n_par_orbs, &
                  & orbs%oe_opt, obj%prim_po(i3)%p, obj%o(i3,obj%i_fe), &
                  & obj%DD2_o(:,4*obj%i_fe-3:4*obj%i_fe), obj%da_o(:,obj%i_fe) )
               else
                  call orbs%orbs(i1)%orb(i2)%cmp_DD2( obj%d_fn(i1,obj%i_fe), obj%y, orbs%n_orbs, obj%prim_po(i3)%p,&
                  & obj%o(i3,obj%i_fe), obj%DD2_o(:,4*obj%i_fe-3:4*obj%i_fe) )
               endif
               i3 = i3 + 1_int32
            enddo ; enddo 
         obj%DD2_o_new(:,1:4) = obj%DD2_o(:,4*obj%i_fe-3:4*obj%i_fe)
      endif
      if ( orbs%oc_opt ) obj%da_o_diff(:) = obj%da_o_diff(:) + obj%da_o(:,obj%i_fe)
   end subroutine upd_orbs_mat
end module orbmat_cls
