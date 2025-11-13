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
module qedip_cls
   use fortran_kinds_v, only: dp, int32
   use jstqepar_var
   use qedip_mod,  only: mat_Aqe, mu_nuc
   use qdo_system_v, only: n_qdo
   use drudonic_config_c, only: drd_cnf
   use molecular_system_v, only: n_fe
   use fermionic_config_c, only: frm_cnf
   implicit none
   type, public :: qedip_t
      real(dp), allocatable, dimension(:)    :: x
      real(dp), allocatable, dimension(:)    :: dx
      real(dp), allocatable, dimension(:)    :: mu
      real(dp), allocatable, dimension(:)    :: dmu
      real(dp), allocatable, dimension(:)    :: w_e
      real(dp), allocatable, dimension(:)    :: dw_e
      real(dp), allocatable, dimension(:)    :: w_d
      real(dp), allocatable, dimension(:)    :: dw_d
      integer(int32)                         :: i_b
      integer(int32)                         :: i_e
   contains
      procedure :: ini     =>  ini_qedip
      procedure :: cmp     =>  cmp_qedip
      procedure :: ratio   =>  ratio_qedip
      procedure :: upd     =>  upd_qedip
      procedure :: cmp_D   =>  cmp_DD2_qedip
      procedure :: new_D   =>  new_DD2_qedip
      procedure :: upd_D   =>  upd_DD2_qedip
      procedure :: cmp_dp  =>  cmp_dp_qedip
      procedure :: upd_dp  =>  upd_dp_qedip
   end type qedip_t
   private  :: ini_qedip, cmp_qedip, ratio_qedip, upd_qedip, &
      cmp_DD2_qedip, new_DD2_qedip, upd_DD2_qedip, &
      cmp_dp_qedip, upd_dp_qedip
contains
   subroutine ini_qedip( obj )
      class(qedip_t), intent(inout) :: obj
      obj%i_b = 0_int32
      obj%i_e = 0_int32
      allocate( obj%x(1:3_int32*n_qdo) )       ;  obj%x     = 0.0_dp
      allocate( obj%dx(1:3_int32) )            ;  obj%dx    = 0.0_dp
      allocate( obj%mu(1:3_int32) )            ;  obj%mu    = 0.0_dp
      allocate( obj%dmu(1:3_int32) )           ;  obj%dmu   = 0.0_dp
      allocate( obj%w_e(1:3_int32*n_qdo) )     ;  obj%w_e   = 0.0_dp
      allocate( obj%dw_e(1:3_int32*n_qdo) )    ;  obj%dw_e  = 0.0_dp
      allocate( obj%w_d(1:3_int32) )           ;  obj%w_d   = 0.0_dp
      allocate( obj%dw_d(1:3_int32) )          ;  obj%dw_d  = 0.0_dp
   end subroutine ini_qedip
   subroutine cmp_qedip( obj, iw, jstqe)
      class(qedip_t),    intent(inout) :: obj
      integer(int32),    intent(in)    :: iw
      real(dp),          intent(inout) :: jstqe
      integer(int32) :: i1, i2, i3
      i3 = 0_int32
      do i1 = 1_int32, n_qdo
         do i2 = 1_int32, 3_int32
            i3 = i3 + 1_int32
            obj%x(i3) = drd_cnf(iw)%d_dq(i1,i1)%v(i2)
         enddo 
      enddo 
      obj%mu(:) = mu_nuc(:)
      do i1 = 1_int32, n_fe
         obj%mu(:) = obj%mu(:) + frm_cnf(iw)%r_fe(:,i1)
      enddo 
      do i1 = 1_int32, 3_int32*n_qdo
         obj%w_e(i1) = dot_product( obj%mu(:), mat_Aqe(:,i1) )
      enddo 
      do i1 = 1_int32, 3_int32
         obj%w_d(i1) = dot_product( mat_Aqe(i1,:), obj%x(:) )
      enddo 
      jstqe = jstqe + dot_product( obj%mu(:), obj%w_d(:) )
   end subroutine cmp_qedip
   subroutine ratio_qedip( obj, iw, g_jstqe )
      class(qedip_t),    intent(inout)  :: obj
      integer(int32),    intent(in)     :: iw
      real(dp),          intent(inout)  :: g_jstqe
      select case ( drd_cnf(iw)%i_type ) 
       case(0) 
         obj%dmu(:)  = frm_cnf(iw)%r_fe_new(:) - frm_cnf(iw)%r_fe(:, frm_cnf(iw)%i_fe )
         g_jstqe = g_jstqe * exp( dot_product( obj%dmu(:), obj%w_d(:) ) )
       case(1) 
         obj%i_e = 3_int32*drd_cnf(iw)%i_drd
         obj%i_b = obj%i_e - 2_int32
         obj%dx(:) = drd_cnf(iw)%d_dq_new( drd_cnf(iw)%i_drd )%v(1:3) - obj%x(obj%i_b : obj%i_e)
         g_jstqe = g_jstqe * exp( dot_product( obj%dx(1:3), obj%w_e( obj%i_b:obj%i_e )  ) )
      end select
   end subroutine ratio_qedip
   subroutine upd_qedip( obj, iw )
      class(qedip_t),    intent(inout)  :: obj
      integer(int32),    intent(in)     :: iw
      integer(int32) :: i1
      select case ( drd_cnf(iw)%i_type ) 
       case(0) 
         obj%mu(:)  = obj%mu(:) + obj%dmu(:)
         do i1 = 1_int32, 3_int32*n_qdo
            obj%dw_e(i1) = dot_product( obj%dmu(:), mat_Aqe(:,i1) )
         enddo 
         obj%w_e(:) = obj%w_e(:) + obj%dw_e(:)
       case(1) 
         obj%x( obj%i_b:obj%i_e ) = obj%x( obj%i_b:obj%i_e ) + obj%dx(:)
         do i1 = 1_int32, 3_int32
            obj%dw_d(i1) = dot_product( mat_Aqe(i1,obj%i_b:obj%i_e), obj%dx(:) )
         enddo 
         obj%w_d(:) = obj%w_d(:) + obj%dw_d(:)
      end select
   end subroutine upd_qedip
   subroutine cmp_DD2_qedip( obj, DD2d_jstqe, DD2e_jstqe )
      class(qedip_t),    intent(inout) :: obj
      real(dp), dimension(4*n_qdo), intent(inout)     :: DD2d_jstqe
      real(dp), dimension(4*n_fe ), intent(inout)     :: DD2e_jstqe
      integer(int32) :: i1
      do i1 = 1_int32, n_qdo
         DD2d_jstqe( (4*i1-3):(4*i1-1) ) = obj%w_e( (3*i1-2):(3*i1) )
         DD2d_jstqe( 4*i1 ) = 0.0
      enddo 
      do i1 = 1_int32, n_fe
         DD2e_jstqe((4*i1-3):(4*i1-1)) = obj%w_d(:)
         DD2e_jstqe(4*i1) = 0.0
      enddo 
   end subroutine cmp_DD2_qedip
   subroutine new_DD2_qedip( obj, iw, DD2d_jstqe, DD2e_jstqe, new_DD2_jstqe )
      class(qedip_t),    intent(inout) :: obj
      integer(int32),    intent(in)     :: iw
      real(dp), dimension(4*n_qdo), intent(in)     :: DD2d_jstqe
      real(dp), dimension(4*n_fe ), intent(in)     :: DD2e_jstqe
      real(dp), dimension(4),       intent(inout)  :: new_DD2_jstqe
      integer(int32) :: im
      select case ( drd_cnf(iw)%i_type ) 
       case(0) 
         im = frm_cnf(iw)%i_fe
         new_DD2_jstqe(1:3) = new_DD2_jstqe(1:3) +  DD2e_jstqe( (4*im-3):(4*im-1) )
         new_DD2_jstqe(4)   = new_DD2_jstqe(4) + 0.0_dp
       case(1) 
         im = drd_cnf(iw)%i_drd
         new_DD2_jstqe(1:3) = new_DD2_jstqe(1:3) + DD2d_jstqe( (4*im-3):(4*im-1) )
         new_DD2_jstqe(4)   = new_DD2_jstqe(4) + 0.0_dp
      end select
   end subroutine new_DD2_qedip
   subroutine upd_DD2_qedip( obj, iw, DD2d_jstqe, DD2e_jstqe )
      class(qedip_t),    intent(inout)  :: obj
      integer(int32),    intent(in)     :: iw
      real(dp), dimension(4*n_qdo), intent(inout)     :: DD2d_jstqe
      real(dp), dimension(4*n_fe ), intent(inout)     :: DD2e_jstqe
      integer(int32) :: i1
      select case ( drd_cnf(iw)%i_type ) ! type of moved particle
       case(0)
         do i1 = 1_int32, n_qdo
            DD2d_jstqe( (4*i1-3):(4*i1-1) ) = DD2d_jstqe( (4*i1-3):(4*i1-1) ) + obj%dw_e( (3*i1-2):(3*i1) )
         enddo 
       case(1) 
         do i1 = 1_int32, n_fe
            DD2e_jstqe( (4*i1-3):(4*i1-1) ) = DD2e_jstqe( (4*i1-3):(4*i1-1) ) + obj%dw_d(:)
         enddo 
      end select
   end subroutine upd_DD2_qedip
   subroutine cmp_dp_qedip( obj, dp_jstqe)
      class(qedip_t),                     intent(inout) :: obj
      real(dp), dimension(n_par_jstqe), intent(inout)   :: dp_jstqe
      integer(int32) :: i1, i2, i3
      i3 = 0_int32
      if(jqe_opt) then
         do i1 = 1_int32, 3_int32
            do i2 = 1_int32, 3_int32*n_qdo
               i3 = i3 + 1_int32
               dp_jstqe(i3) = obj%mu(i1)*obj%x(i2)
            enddo 
         enddo 
      endif
   end subroutine cmp_dp_qedip
   subroutine upd_dp_qedip( obj, dp_jstqe )
      class(qedip_t),    intent(inout)  :: obj
      real(dp), dimension(n_par_jstqe), intent(inout)   :: dp_jstqe
      integer(int32) :: i1, i2, i3
      i3 = 0_int32
      if(jqe_opt) then
         do i1 = 1_int32, 3_int32
            do i2 = 1_int32, 3_int32*n_qdo
               i3 = i3 + 1_int32
               dp_jstqe(i3) = obj%mu(i1)*obj%x(i2)
            enddo 
         enddo 
      endif
   end subroutine upd_dp_qedip
end module qedip_cls
