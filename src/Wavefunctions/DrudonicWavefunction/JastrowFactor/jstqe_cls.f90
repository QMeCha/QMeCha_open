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
module jstqe_cls
   use fortran_kinds_v, only: dp, int32
   use jstqepar_var
   use qedip_cls,  only: qedip_t
   implicit none
   type, public :: jstqefct_t
      real(dp)                            :: jstqe
      real(dp)                            :: g_jstqe
      real(dp), allocatable, dimension(:) :: DD2d_jstqe
      real(dp), allocatable, dimension(:) :: DD2e_jstqe
      real(dp), allocatable, dimension(:) :: new_DD2_jstqe
      real(dp), allocatable, dimension(:) :: dp_jstqe
      type(qedip_t)                       :: qedip
   contains
      procedure :: ini    => ini_jstqe_fct
      procedure :: cmp    => cmp_jstqe_fct
      procedure :: ratio  => ratio_jstqe_fct
      procedure :: upd    => upd_jstqe_fct
      procedure :: cmp_D  => cmp_DD2_jstqe_fct
      procedure :: new_D  => new_DD2_jstqe_fct
      procedure :: upd_D  => upd_DD2_jstqe_fct
      procedure :: cmp_dp => cmp_dp_jstqe_fct
      procedure :: upd_dp => upd_dp_jstqe_fct
   end type jstqefct_t
   private :: ini_jstqe_fct, cmp_jstqe_fct, ratio_jstqe_fct, upd_jstqe_fct, &
      cmp_DD2_jstqe_fct, new_DD2_jstqe_fct, upd_DD2_jstqe_fct, &
      cmp_dp_jstqe_fct, upd_dp_jstqe_fct
contains
   subroutine ini_jstqe_fct( obj )
      use molecular_system_v, only: n_el
      use qdo_system_v, only: n_qdo
      use quantum_monte_carlo_v, only: grad_updat
      class(jstqefct_t), intent(inout) :: obj
      obj%jstqe = 0.0_dp
      allocate( obj%DD2d_jstqe( 1:4*n_qdo ) );        obj%DD2d_jstqe     = 0.0_dp
      allocate( obj%DD2e_jstqe( 1:4*n_el  ) );        obj%DD2e_jstqe     = 0.0_dp
      if ( grad_updat) then
         allocate( obj%new_DD2_jstqe( 1:4 ) );         obj%new_DD2_jstqe  = 0.0_dp
      endif
      if ( jqe_opt ) then
         allocate( obj%dp_jstqe( 1:n_par_jstqe ) );      obj%dp_jstqe       = 0.0_dp
      endif
      select case (jstqe_type)
       case (1) 
         call obj%qedip%ini()
      end select
   end subroutine ini_jstqe_fct
   subroutine cmp_jstqe_fct( obj, iw )
      class(jstqefct_t), intent(inout) :: obj
      integer(int32),    intent(in)    :: iw
      obj%jstqe = 0.0_dp
      select case (jstqe_type)
       case (1) 
         call obj%qedip%cmp( iw, obj%jstqe )
      end select
      obj%jstqe = exp(obj%jstqe)
   end subroutine cmp_jstqe_fct
   subroutine ratio_jstqe_fct( obj, iw, g )
      class(jstqefct_t), intent(inout) :: obj
      integer(int32),    intent(in)    :: iw
      real(dp),          intent(out)   :: g
      obj%g_jstqe = 1.0_dp
      select case (jstqe_type)
       case(1) 
         call obj%qedip%ratio( iw, obj%g_jstqe )
      end select
      g = obj%g_jstqe
   end subroutine ratio_jstqe_fct
   subroutine upd_jstqe_fct( obj, iw )
      class(jstqefct_t), intent(inout) :: obj
      integer(int32),    intent(in)    :: iw
      select case (jstqe_type)
       case(1)
         call obj%qedip%upd( iw )
      end select
      obj%jstqe = obj%jstqe * obj%g_jstqe
   end subroutine upd_jstqe_fct
   subroutine cmp_DD2_jstqe_fct( obj )
      class(jstqefct_t), intent(inout) :: obj
      select case (jstqe_type)
       case (1) 
         call obj%qedip%cmp_D( obj%DD2d_jstqe, obj%DD2e_jstqe )
      end select
   end subroutine cmp_DD2_jstqe_fct
   subroutine new_DD2_jstqe_fct( obj , iw )
      class(jstqefct_t), intent(inout) :: obj
      integer(int32),    intent(in)    :: iw
      obj%new_DD2_jstqe = 0.0_dp
      select case (jstqe_type)
       case (1) 
         call obj%qedip%new_D( iw, obj%DD2d_jstqe, obj%DD2e_jstqe, obj%new_DD2_jstqe )
      end select
   end subroutine new_DD2_jstqe_fct
   subroutine upd_DD2_jstqe_fct( obj , iw )
      class(jstqefct_t), intent(inout) :: obj
      integer(int32),    intent(in)    :: iw
      select case (jstqe_type)
       case (1) 
         call obj%qedip%upd_D( iw, obj%DD2d_jstqe, obj%DD2e_jstqe )
      end select
   end subroutine upd_DD2_jstqe_fct
   subroutine cmp_dp_jstqe_fct( obj, dp_ln_jstqe )
      class(jstqefct_t), intent(inout) :: obj
      real(dp), dimension(n_par_jstqe), intent(inout) :: dp_ln_jstqe
      select case (jstqe_type)
       case(1) 
         call obj%qedip%cmp_dp( obj%dp_jstqe )
      end select
      dp_ln_jstqe(:) = obj%dp_jstqe(:)
   end subroutine cmp_dp_jstqe_fct
   subroutine upd_dp_jstqe_fct( obj, dp_ln_jstqe )
      class(jstqefct_t), intent(inout) :: obj
      real(dp), dimension(n_par_jstqe), intent(inout) :: dp_ln_jstqe
      select case (jstqe_type)
       case(1) 
         call obj%qedip%upd_dp( obj%dp_jstqe )
      end select
      dp_ln_jstqe(:) = obj%dp_jstqe(:)
   end subroutine upd_dp_jstqe_fct
end module jstqe_cls
