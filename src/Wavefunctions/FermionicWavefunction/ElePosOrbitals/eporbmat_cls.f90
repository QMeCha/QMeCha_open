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
module eporbmat_cls
   use fortran_kinds_v, only: dp, int32
   use cubic_harmonics_m, only: cubhrm_t
   use atomic_orbital_c, only: orb_t
   use molecular_system_v, only: dist_t, n_fe, n_el, n_po
   use fermionic_config_c, only: frm_cnf
   use quantum_monte_carlo_v, only: grad_updat
   use eporbbss_cls, only: epsysorb_t
   use cubic_harmonics_m, only: y
   implicit none
   type, public :: prim_po_t
      real(dp), allocatable, dimension(:) :: p
   end type prim_po_t
   type, public :: epmatorb_t
      integer(int32),                 pointer :: i_fe
      integer(int32),                 pointer :: l_max
      type(dist_t),  dimension(:,:),  pointer :: d_ep
      type(dist_t),  dimension(:),    pointer :: d_ff_new
      type(cubhrm_t),                 pointer :: y
      logical,                        pointer :: dc_opt
      logical,                        pointer :: de_opt
      integer(int32),                 pointer :: n_orbs
      integer(int32),                 pointer :: n_par_orbs
      type(orb_t), dimension(:),      pointer :: orb
      type(prim_po_t), allocatable, dimension(:) :: prim_po
      real(dp), allocatable, dimension(:,:) :: o
      real(dp), allocatable, dimension(:,:) :: o_new
      real(dp), allocatable, dimension(:,:) :: o_dif
      real(dp), allocatable, dimension(:,:) :: DD2_o
      real(dp), allocatable, dimension(:,:) :: DD2_o_new
      real(dp), allocatable, dimension(:,:) :: da_o
      real(dp), allocatable, dimension(:,:) :: da_o_diff
      real(dp), allocatable, dimension(:,:) :: da_DD2_o
   contains
      procedure :: ini => ini_orbs_mat
      procedure :: cmp => cmp_orbs_mat
      procedure :: new => new_orbs_vec
      procedure :: new_D => new_orbD_vec
      procedure :: upd => upd_orbs_mat
   end type epmatorb_t
contains
   subroutine ini_orbs_mat( obj, iw, orbs )
      class(epmatorb_t),         intent(inout) :: obj
      integer(int32),            intent(in)    :: iw
      type(epsysorb_t),  target, intent(in)    :: orbs
      integer(int32) :: i1, i2, i3
      obj%n_orbs      => orbs%n_orbs
      obj%n_par_orbs  => orbs%n_par_orbs
      obj%l_max       => orbs%l_max
      obj%orb      => orbs%orb
      obj%dc_opt   => orbs%oc_opt
      obj%de_opt   => orbs%oe_opt
      obj%y        => y(iw)
      obj%d_ep     => frm_cnf(iw)%d_ep
      obj%d_ff_new => frm_cnf(iw)%d_ff_new
      obj%i_fe     => frm_cnf(iw)%i_fe
      allocate( obj%o(1:obj%n_orbs,1:n_el*n_po) )       ; obj%o = 0.0_dp
      allocate( obj%DD2_o(1:obj%n_orbs,1:4*n_el*n_po) ) ; obj%DD2_o = 0.0_dp
      if ( obj%n_par_orbs.gt.0_int32 ) then
         allocate( obj%da_o(1:obj%n_par_orbs,1:n_el*n_po) ) ; obj%da_o = 0.0_dp
      endif
      if (n_el.gt.n_po) then
         allocate( obj%o_new(1:obj%n_orbs, 1:n_el) )  ; obj%o_new = 0.0_dp
         allocate( obj%o_dif(1:obj%n_orbs, 1:n_el) )  ; obj%o_dif = 0.0_dp
         allocate( obj%DD2_o_new(1:obj%n_orbs,1:4*n_el) ) ; obj%DD2_o_new = 0.0_dp
         allocate( obj%prim_po(1:obj%n_orbs*n_el) )
         i3 = 1_int32
         do i1 = 1_int32, n_el ; do i2 = 1_int32, obj%n_orbs
               allocate( obj%prim_po(i3)%p(1:obj%orb(i2)%n) ) ; obj%prim_po(i3)%p = 0.0_dp
               i3 = i3 + 1_int32
            enddo ; enddo 
         if ( obj%n_par_orbs.gt.0_int32 ) then
            allocate( obj%da_o_diff(1:obj%n_par_orbs,1:n_el) ) ; obj%da_o_diff = 0.0_dp
         endif
      else
         allocate( obj%o_new(1:obj%n_orbs, 1:n_po) )  ; obj%o_new = 0.0_dp
         allocate( obj%o_dif(1:obj%n_orbs, 1:n_po) )  ; obj%o_dif = 0.0_dp
         allocate( obj%DD2_o_new(1:obj%n_orbs,1:4*n_po) ) ; obj%DD2_o_new = 0.0_dp
         allocate( obj%prim_po(1:obj%n_orbs*n_po) )
         i3 = 1_int32
         do i1 = 1_int32, n_po ; do i2 = 1_int32, obj%n_orbs
               allocate( obj%prim_po(i3)%p(1:obj%orb(i2)%n) ) ; obj%prim_po(i3)%p = 0.0_dp
               i3 = i3 + 1_int32
            enddo ; enddo 
         if ( obj%n_par_orbs.gt.0_int32 ) then
            allocate( obj%da_o_diff(1:obj%n_par_orbs,1:n_po) ) ; obj%da_o_diff = 0.0_dp
         endif
      endif
   end subroutine ini_orbs_mat
   subroutine cmp_orbs_mat( obj )
      class(epmatorb_t), intent(inout) :: obj
      integer(int32) :: i1, i2, i3, o1
      i3 = 0_int32
      do i2 = 1_int32, n_el
         do i1 = 1_int32, n_po
            i3 = i3 + 1_int32
            call obj%y%cmp( obj%d_ep(i1,i2), obj%l_max )
            call obj%y%cmp_D( obj%d_ep(i1,i2), obj%l_max )
            do o1 = 1_int32, obj%n_orbs
               call obj%orb(o1)%cmp( obj%d_ep(i1,i2)%m, obj%y, obj%prim_po(o1)%p, obj%o(o1,i3) )
               if ( obj%dc_opt ) then
                  call obj%orb(o1)%cmp_der( obj%d_ep(i1,i2), obj%y, obj%n_orbs, obj%n_par_orbs, &
                  & obj%de_opt, obj%prim_po(o1)%p, obj%o(o1,i3), &
                  & obj%DD2_o(:,4*i3-3:4*i3), obj%da_o(:,i3) )
               else
                  call obj%orb(o1)%cmp_DD2( obj%d_ep(i1,i2), obj%y, obj%n_orbs, obj%prim_po(o1)%p,&
                  & obj%o(o1,i3), obj%DD2_o(:,4*i3-3:4*i3) )
               endif
            enddo 
         enddo 
      enddo 
   end subroutine cmp_orbs_mat
   subroutine new_orbs_vec( obj )
      class(epmatorb_t), intent(inout) :: obj
      type(dist_t) :: d_ff_tmp
      integer(int32) :: i1, i2, o1, o2
      o2 = 0_int32
      if ( obj%i_fe.le. n_el ) then
         do i2 = 1_int32, n_po
            call obj%y%cmp( obj%d_ff_new(i2+n_el), obj%l_max )
            do o1 = 1_int32, obj%n_orbs
               o2 = o2 + 1_int32
               call obj%orb(o1)%cmp( obj%d_ff_new(i2+n_el)%m, obj%y, obj%prim_po(o2)%p, obj%o_new(o1,i2) )
            enddo 
         enddo 
         i1 = n_po*(obj%i_fe-1)
         obj%o_dif(:,1:n_po) = obj%o_new(:,1:n_po) - obj%o(:,i1+1:i1+n_po)
      else 
         do i2 = 1_int32, n_el
            d_ff_tmp%v = -obj%d_ff_new(i2)%v
            d_ff_tmp%m =  obj%d_ff_new(i2)%m
            call obj%y%cmp( d_ff_tmp, obj%l_max )
            do o1 = 1_int32, obj%n_orbs
               o2 = o2 + 1_int32
               call obj%orb(o1)%cmp( d_ff_tmp%m, obj%y, obj%prim_po(o2)%p, obj%o_new(o1,i2) )
            enddo 
         enddo 
         i2 = obj%i_fe - n_el
         do i1 = 1_int32, n_el
            obj%o_dif(:,i1) = obj%o_new(:,i1) - obj%o(:,i2)
            i2 = i2 + n_po
         enddo
      endif
   end subroutine new_orbs_vec
   subroutine new_orbD_vec( obj )
      class(epmatorb_t), intent(inout) :: obj
      type(dist_t) :: d_ff_tmp
      integer(int32) :: i1, i2, o1, o2
      i2 = 0_int32
      o2 = 0_int32
      if ( obj%i_fe.le.n_el ) then
         do i1 = n_el+1_int32, n_fe
            i2 = i2 + 1_int32
            call obj%y%cmp(   obj%d_ff_new(i1), obj%l_max )
            call obj%y%cmp_D( obj%d_ff_new(i1), obj%l_max )
            do o1 = 1_int32, obj%n_orbs
               o2 = o2 + 1_int32
               call obj%orb(o1)%cmp( obj%d_ff_new(i1)%m, obj%y, obj%prim_po(o2)%p, obj%o_new(o1,i2) )
               call obj%orb(o1)%cmp_DD2(obj%d_ff_new(i1), obj%y, obj%n_orbs, obj%prim_po(o2)%p,&
               & obj%o_new(o1,i2), obj%DD2_o_new(:,4*i2-3:4*i2) )
            enddo 
         enddo
         i1 = n_po*(obj%i_fe-1)
         obj%o_dif(:,1:n_po) = obj%o_new(:,1:n_po) - obj%o(:,i1+1:i1+n_po)
      else
         do i1 = 1_int32, n_el
            i2 = i2 + 1_int32
            d_ff_tmp%v = -obj%d_ff_new(i1)%v
            d_ff_tmp%m =  obj%d_ff_new(i1)%m
            call obj%y%cmp(   d_ff_tmp, obj%l_max )
            call obj%y%cmp_D( d_ff_tmp, obj%l_max )
            do o1 = 1_int32, obj%n_orbs
               o2 = o2 + 1_int32
               call obj%orb(o1)%cmp( d_ff_tmp%m, obj%y, obj%prim_po(o2)%p, obj%o_new(o1,i2) )
               call obj%orb(o1)%cmp_DD2( d_ff_tmp, obj%y, obj%n_orbs, obj%prim_po(o2)%p,&
               & obj%o_new(o1,i2), obj%DD2_o_new(:,4*i2-3:4*i2) )
            enddo 
         enddo
         i2 = obj%i_fe - n_el
         do i1 = 1_int32, n_el
            obj%o_dif(:,i1) = obj%o_new(:,i1) - obj%o(:,i2)
            i2 = i2 + n_po
         enddo
      endif
   end subroutine new_orbD_vec
   subroutine upd_orbs_mat( obj )
      class(epmatorb_t), intent(inout) :: obj
      integer(int32) :: i1, i2, o1, o2
      if ( obj%i_fe .le. n_el ) then
         i1 = n_po*(obj%i_fe-1)
         obj%o(:,i1+1:i1+n_po) = obj%o(:,i1+1:i1+n_po) + obj%o_dif(:,1:n_po)
         if ( obj%dc_opt ) obj%da_o_diff(:,1:n_po) = - obj%da_o(:,i1+1:i1+n_po)
         if (grad_updat) then
            obj%DD2_o(:,4*i1+1:4*i1+4*n_po) = obj%DD2_o_new(:,1:4*n_po)
            if ( obj%dc_opt ) then
               o2 = 0_int32
               do i2 = 1_int32, n_po
                  i1 = i1 + 1_int32
                  do o1 = 1_int32, obj%n_orbs
                     o2 = o2 + 1_int32
                     call obj%orb(o1)%cmp_dc( obj%d_ep(i2,obj%i_fe), obj%n_par_orbs, &
                     & obj%de_opt, obj%prim_po(o2)%p, obj%da_o(:,i1) )
                  enddo 
               enddo 
            endif
         else
            i1 = n_po*(obj%i_fe-1)
            o2 = 0_int32
            do i2 = 1_int32, n_po
               i1 = i1 + 1_int32
               call obj%y%cmp(   obj%d_ep(i2,obj%i_fe), obj%l_max )
               call obj%y%cmp_D( obj%d_ep(i2,obj%i_fe), obj%l_max )
               do o1 = 1_int32, obj%n_orbs
                  o2 = o2 + 1_int32
                  if ( obj%dc_opt ) then
                     call obj%orb(o1)%cmp_der( obj%d_ep(i2,obj%i_fe), obj%y, obj%n_orbs, obj%n_par_orbs, &
                     & obj%de_opt, obj%prim_po(o2)%p, obj%o(o1,i1), &
                     & obj%DD2_o(:,4*i1-3:4*i1), obj%da_o(:,i1) )
                  else
                     call obj%orb(o1)%cmp_DD2( obj%d_ep(i2,obj%i_fe), obj%y, obj%n_orbs, obj%prim_po(o2)%p,&
                     & obj%o(o1,i1), obj%DD2_o(:,4*i1-3:4*i1) )
                  endif
               enddo 
               if ( obj%dc_opt ) obj%da_o_diff(:,i2) = obj%da_o_diff(:,i2) + obj%da_o(:,i1)
            enddo 
         endif
      else 
         i2 = obj%i_fe - n_el
         do i1 = 1_int32, n_el
            obj%o(:,i2) = obj%o(:,i2) + obj%o_dif(:,i1)
            if ( obj%dc_opt ) obj%da_o_diff(:,i1) = - obj%da_o(:,i2)
            i2 = i2 + n_po
         enddo 
         if (grad_updat) then
            i2 = obj%i_fe - n_el
            do i1 = 1_int32, n_el
               obj%DD2_o(:,4*i2-3:4*i2) = obj%DD2_o_new(:,4*i1-3:4*i1)
               i2 = i2 + n_po
            enddo 
            if ( obj%dc_opt ) then
               o2 = 0_int32
               i2 = obj%i_fe - n_po
               do i1 = 1_int32, n_el
                  i2 = i2 + n_po
                  do o1 = 1_int32, obj%n_orbs
                     o2 = o2 + 1_int32
                     call obj%orb(o1)%cmp_dc( obj%d_ep(obj%i_fe-n_el,i1), obj%n_par_orbs, &
                     & obj%de_opt, obj%prim_po(o2)%p, obj%da_o(:,i2) )
                  enddo 
               enddo 
            endif 
         else 
            o2 = 0_int32
            i2 = obj%i_fe - n_el
            do i1 = 1_int32, n_el
               call obj%y%cmp(   obj%d_ep(obj%i_fe-n_el,i1), obj%l_max )
               call obj%y%cmp_D( obj%d_ep(obj%i_fe-n_el,i1), obj%l_max )
               do o1 = 1_int32, obj%n_orbs
                  o2 = o2 + 1_int32
                  if ( obj%dc_opt ) then
                     call obj%orb(o1)%cmp_der( obj%d_ep(obj%i_fe-n_el,i1), obj%y, obj%n_orbs, obj%n_par_orbs, &
                     & obj%de_opt, obj%prim_po(o2)%p, obj%o(o1,i2), &
                     & obj%DD2_o(:,4*i2-3:4*i2), obj%da_o(:,i2) )
                  else
                     call obj%orb(o1)%cmp_DD2( obj%d_ep(obj%i_fe-n_el,i1), obj%y, obj%n_orbs, obj%prim_po(o2)%p,&
                     & obj%o(o1,i2), obj%DD2_o(:,4*i2-3:4*i2) )
                  endif
               enddo 
               if ( obj%dc_opt ) obj%da_o_diff(:,i1) = obj%da_o_diff(:,i1) + obj%da_o(:,i2)
               i2 = i2 + n_po
            enddo 
         endif
      endif 
   end subroutine upd_orbs_mat
end module eporbmat_cls
