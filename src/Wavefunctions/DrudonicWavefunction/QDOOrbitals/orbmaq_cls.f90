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
module orbmaq_cls
   use fortran_kinds_v, only: dp, int32
   use cubic_harmonics_m, only: cubhrm_t
   use orbspq_cls, only: orbspq_t
   use molecular_system_v, only: dist_t
   use orbmat_cls, only: prim_po_t
   use qdo_system_v, only: n_qdo
   implicit none
   type, public :: matorbq_t
      integer(int32),                 pointer :: i_drd
      type(dist_t),  dimension(:,:),  pointer :: d_dq
      type(dist_t),  dimension(:),    pointer :: d_dq_new
      type(cubhrm_t),                 pointer :: y
      logical,                        pointer :: dc_opt
      logical,                        pointer :: de_opt
      integer(int32),                 pointer :: n_orbs
      integer(int32),                 pointer :: n_par_orbs
      type(orbspq_t), dimension(:), pointer :: orbs
      type(prim_po_t), allocatable, dimension(:) :: prim_po
      real(dp), allocatable, dimension(:,:) :: o
      real(dp), allocatable, dimension(:)   :: o_new
      real(dp), allocatable, dimension(:,:) :: DD2_o
      real(dp), allocatable, dimension(:,:) :: DD2_o_new
      real(dp), allocatable, dimension(:,:) :: da_o
   contains
      procedure :: ini => ini_orbs_mat
      procedure :: cmp => cmp_orbs_mat
      procedure :: new => new_orbs_vec
      procedure :: new_D => new_orbD_vec
      procedure :: upd => upd_orbs_mat
   end type matorbq_t
contains
   subroutine ini_orbs_mat( obj, iw, orbs )
      use drudonic_config_c, only: drd_cnf
      use orbbssq_cls, only: qsysorb_t
      use cubic_harmonics_m, only: y
      class(matorbq_t),      intent(inout) :: obj
      integer(int32),          intent(in)    :: iw
      type(qsysorb_t),  target, intent(in) :: orbs
      integer(int32) :: i1, i2, i3
      obj%n_orbs      => orbs%n_orbs
      obj%n_par_orbs  => orbs%n_par_orbs
      obj%orbs      => orbs%orbs
      obj%dc_opt    => orbs%oc_opt
      obj%de_opt    => orbs%oe_opt
      obj%y         => y(iw)
      obj%d_dq      => drd_cnf(iw)%d_dq
      obj%d_dq_new  => drd_cnf(iw)%d_dq_new
      obj%i_drd     => drd_cnf(iw)%i_drd
      allocate( obj%prim_po(1:obj%n_orbs) )
      i3 = 1_int32
      do i1 = 1_int32, n_qdo ; do i2 = 1_int32, obj%orbs(i1)%n_orbs
            allocate( obj%prim_po(i3)%p(1:obj%orbs(i1)%orb(i2)%n) ) ; obj%prim_po(i3)%p = 0.0_dp
            i3 = i3 + 1_int32
         enddo ; enddo
      allocate( obj%o(1:obj%n_orbs,1:n_qdo) )           ; obj%o         = 0.0_dp
      allocate( obj%o_new(1:obj%n_orbs) )               ; obj%o_new     = 0.0_dp
      allocate( obj%DD2_o(1:obj%n_orbs,1:4*n_qdo) )     ; obj%DD2_o     = 0.0_dp
      allocate( obj%DD2_o_new(1:obj%n_orbs,1:4) )       ; obj%DD2_o_new = 0.0_dp
      if ( obj%n_par_orbs.gt.0_int32 ) then
         allocate( obj%da_o(1:obj%n_par_orbs,1:n_qdo) ) ; obj%da_o      = 0.0_dp
      endif
   end subroutine ini_orbs_mat
   subroutine cmp_orbs_mat( obj )
      class(matorbq_t), intent(inout) :: obj
      integer(int32) :: i1, i2, i3, i4
      do i1 = 1_int32, n_qdo
         i4 = 1_int32
         do i2 = 1_int32, n_qdo
            call obj%y%cmp( obj%d_dq(i2,i1), obj%orbs(i2)%l_max )
            call obj%y%cmp_D( obj%d_dq(i2,i1), obj%orbs(i2)%l_max )
            do i3 = 1_int32, obj%orbs(i2)%n_orbs
               call obj%orbs(i2)%orb(i3)%cmp( obj%d_dq(i2,i1)%m, obj%y, obj%prim_po(i4)%p, obj%o(i4,i1) )
               if ( obj%dc_opt ) then
                  call obj%orbs(i2)%orb(i3)%cmp_der( obj%d_dq(i2,i1), obj%y, obj%n_orbs, obj%n_par_orbs, &
                  & obj%de_opt, obj%prim_po(i4)%p, obj%o(i4,i1), &
                  & obj%DD2_o(:,4*i1-3:4*i1), obj%da_o(:,i1) )
               else
                  call obj%orbs(i2)%orb(i3)%cmp_DD2( obj%d_dq(i2,i1), obj%y, obj%n_orbs, obj%prim_po(i4)%p,&
                  & obj%o(i4,i1), obj%DD2_o(:,4*i1-3:4*i1) )
               endif
               i4 = i4 + 1_int32
            enddo ; enddo ; enddo
   end subroutine cmp_orbs_mat
   subroutine new_orbs_vec( obj )
      class(matorbq_t), intent(inout) :: obj
      integer(int32) :: i1, i2, i3
      i3 = 1_int32
      do i1 = 1_int32, n_qdo
         call obj%y%cmp( obj%d_dq_new(i1), obj%orbs(i1)%l_max )
         do i2 = 1_int32, obj%orbs(i1)%n_orbs
            call obj%orbs(i1)%orb(i2)%cmp( obj%d_dq_new(i1)%m, obj%y, obj%prim_po(i3)%p, obj%o_new(i3) )
            i3 = i3 + 1_int32
         enddo ; enddo 
   end subroutine new_orbs_vec
   subroutine new_orbD_vec( obj )
      class(matorbq_t), intent(inout) :: obj
      integer(int32) :: i1, i2, i3
      i3 = 1_int32
      do i1 = 1_int32, n_qdo
         call obj%y%cmp( obj%d_dq_new(i1), obj%orbs(i1)%l_max )
         call obj%y%cmp_D( obj%d_dq_new(i1), obj%orbs(i1)%l_max )
         do i2 = 1_int32, obj%orbs(i1)%n_orbs
            call obj%orbs(i1)%orb(i2)%cmp( obj%d_dq_new(i1)%m, obj%y, obj%prim_po(i3)%p, obj%o_new(i3) )
            call obj%orbs(i1)%orb(i2)%cmp_DD2( obj%d_dq_new(i1), obj%y, obj%n_orbs, obj%prim_po(i3)%p,&
            & obj%o_new(i3), obj%DD2_o_new(:,1:4) )
            i3 = i3 + 1_int32
         enddo ; enddo 
   end subroutine new_orbD_vec
   subroutine upd_orbs_mat( obj )
      use quantum_monte_carlo_v, only: grad_updat
      class(matorbq_t), intent(inout) :: obj
      integer(int32), pointer :: idrd
      integer(int32) :: i1, i2, i3
      idrd => obj%i_drd
      obj%o(:,idrd) = obj%o_new
      if ( grad_updat ) then
         obj%DD2_o(:,4*idrd-3:4*idrd) = obj%DD2_o_new(:,1:4)
         if ( obj%dc_opt ) then
            i3 = 1_int32
            do i1 = 1_int32, n_qdo
               do i2 = 1_int32, obj%orbs(i1)%n_orbs
                  call obj%orbs(i1)%orb(i2)%cmp_dc( obj%d_dq(i1,idrd), obj%n_par_orbs, &
                  & obj%de_opt, obj%prim_po(i3)%p, obj%da_o(:,idrd) )
                  i3 = i3 + 1_int32
               enddo ; enddo 
         endif
      else
         i3 = 1_int32
         do i1 = 1_int32, n_qdo
            call obj%y%cmp( obj%d_dq(i1,idrd), obj%orbs(i1)%l_max )
            call obj%y%cmp_D( obj%d_dq(i1,idrd), obj%orbs(i1)%l_max )
            do i2 = 1_int32, obj%orbs(i1)%n_orbs
               if ( obj%dc_opt ) then
                  call obj%orbs(i1)%orb(i2)%cmp_der( obj%d_dq(i1,idrd), obj%y, obj%n_orbs, obj%n_par_orbs, &
                  & obj%de_opt, obj%prim_po(i3)%p, obj%o(i3,idrd), &
                  & obj%DD2_o(:,4*idrd-3:4*idrd), obj%da_o(:,idrd) )
               else
                  call obj%orbs(i1)%orb(i2)%cmp_DD2( obj%d_dq(i1,idrd), obj%y, obj%n_orbs, obj%prim_po(i3)%p,&
                  & obj%o(i3,idrd), obj%DD2_o(:,4*idrd-3:4*idrd) )
               endif
               i3 = i3 + 1_int32
            enddo ; enddo 
         obj%DD2_o_new(:,1:4) = obj%DD2_o(:,4*idrd-3:4*idrd)
      endif
   end subroutine upd_orbs_mat
end module orbmaq_cls
