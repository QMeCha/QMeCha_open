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
module dipqfun_cls
   use fortran_kinds_v,   only: dp, int32
   use qdo_wavefunction_v,        only: ql_opt
   use dipqcfs_mod,       only: dipqcfs_t
   use drudonic_config_c, only: drdcnf_t
   use openmp_mpi_m
   implicit none
   type, public :: dipqfun_t
      real(dp)                               :: psi_dip
      real(dp)                               :: r
      real(dp), allocatable, dimension(:)    :: x
      real(dp), allocatable, dimension(:)    :: dx
      real(dp), allocatable, dimension(:)    :: v
      real(dp), allocatable, dimension(:)    :: dv
      integer(int32)                         :: i_b
      integer(int32)                         :: i_e
      real(dp), allocatable, dimension(:)    :: DD2_ln_dipq
      real(dp), allocatable, dimension(:)    :: DD2_ln_dipq_new
   contains
      procedure :: ini    => ini_dipq_fun
      procedure :: cmp    => cmp_dipq_fun
      procedure :: ratio  => ratio_dipq_fun
      procedure :: upd    => upd_dipq_fun
      procedure :: cmp_D  => cmp_DD2_dipq_fun
      procedure :: new_D  => new_DD2_dipq_fun
      procedure :: upd_D  => upd_DD2_dipq_fun
      procedure :: cmp_dl => cmp_dl_dipq_fun
   end type dipqfun_t
   private :: ini_dipq_fun, cmp_dipq_fun, ratio_dipq_fun, upd_dipq_fun, &
   & cmp_DD2_dipq_fun, new_DD2_dipq_fun, upd_DD2_dipq_fun, &
   & cmp_dl_dipq_fun
contains
   subroutine ini_dipq_fun( obj, n_qdo )
      use quantum_monte_carlo_v, only: grad_updat
      class(dipqfun_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: n_qdo
      obj%psi_dip = 0.0_dp
      obj%i_b     = 0_int32
      allocate( obj%x(1:3_int32*n_qdo) )               ; obj%x             = 0.0_dp
      allocate( obj%dx(1:3_int32) )                    ; obj%dx            = 0.0_dp
      allocate( obj%v(1:3_int32*n_qdo) )               ; obj%v             = 0.0_dp
      allocate( obj%dv(1:3_int32*n_qdo) )              ; obj%dv            = 0.0_dp
      allocate( obj%DD2_ln_dipq(1:4*n_qdo) )         ; obj%DD2_ln_dipq   = 0.0_dp
      if (grad_updat) then
         allocate( obj%DD2_ln_dipq_new(1:4) )           ; obj%DD2_ln_dipq_new = 0.0_dp
      endif
   end subroutine ini_dipq_fun
   subroutine cmp_dipq_fun( obj, n_qdo, dipq_c, drd_cnf)
      class(dipqfun_t), intent(inout) :: obj
      integer(int32),                             intent(in) :: n_qdo
      type(dipqcfs_t),                            intent(in) :: dipq_c
      type(drdcnf_t),                             intent(in) :: drd_cnf
      integer(int32) :: i1
      do i1 = 1_int32, n_qdo
         obj%x(3*i1-2:3*i1) = drd_cnf%d_dq(i1,i1)%v(:) - dipq_c%vec_xi(3*i1-2:3*i1)
      enddo 
      call DGEMM('n','n',3_int32*n_qdo,1_int32,3_int32*n_qdo,1.0_dp, &
      &dipq_c%mat_A,3_int32*n_qdo,obj%x,3_int32*n_qdo,0.0_dp,obj%v,3_int32*n_qdo)
      obj%psi_dip = dot_product(obj%x, obj%v)
      obj%psi_dip = exp(-obj%psi_dip)
   end subroutine cmp_dipq_fun
   subroutine ratio_dipq_fun( obj, n_qdo, dipq_c, drd_cnf, g)
      class(dipqfun_t), intent(inout) :: obj
      integer(int32),                            intent(in)  :: n_qdo
      type(dipqcfs_t),                           intent(in)  :: dipq_c
      type(drdcnf_t),                            intent(in)  :: drd_cnf
      real(dp),                                  intent(out) :: g
      obj%i_b = (drd_cnf%i_drd - 1_int32)*3_int32 + 1_int32
      obj%i_e = obj%i_b + 2_int32
      obj%dx(:) = drd_cnf%d_dq_new( drd_cnf%i_drd )%v(:) - dipq_c%vec_xi(obj%i_b : obj%i_e) - obj%x(obj%i_b : obj%i_e)
      call DGEMM('n','n',3_int32*n_qdo,1_int32,3_int32,1.0_dp, \
      dipq_c%mat_A(:,obj%i_b:obj%i_e),3_int32*n_qdo,obj%dx,3_int32,0.0_dp,obj%dv,3_int32*n_qdo)
      g = exp(-dot_product( obj%dx(:), 2*obj%v(obj%i_b:obj%i_e) + obj%dv(obj%i_b:obj%i_e) ))
      obj%r = g
   end subroutine ratio_dipq_fun
   subroutine upd_dipq_fun( obj )
      class(dipqfun_t), intent(inout) :: obj
      obj%x( obj%i_b:obj%i_e ) = obj%x( obj%i_b:obj%i_e ) + obj%dx(:)
      obj%v(:) = obj%v(:) + obj%dv(:)
      obj%psi_dip = obj%r * obj%psi_dip
   end subroutine upd_dipq_fun
   subroutine cmp_DD2_dipq_fun( obj, n_qdo, dipq_c)
      class(dipqfun_t), intent(inout) :: obj
      integer(int32),                            intent(in)  :: n_qdo
      type(dipqcfs_t),                           intent(in)  :: dipq_c
      integer(int32) :: i1
      do i1 = 1_int32, n_qdo
         obj%DD2_ln_dipq( (4*i1-3):(4*i1-1) ) = -2.0_dp*obj%v( (3*i1-2):(3*i1)  )
         obj%DD2_ln_dipq( 4*i1 ) = - 2.0_dp*( dipq_c%mat_A( (3*i1-2), (3*i1-2) ) + &
            dipq_c%mat_A( (3*i1-1), (3*i1-1) ) + &
            dipq_c%mat_A( (3*i1-0), (3*i1-0) ) )
      enddo 
   end subroutine cmp_DD2_dipq_fun
   subroutine new_DD2_dipq_fun( obj, i_drd, dipq_c)
      class(dipqfun_t), intent(inout) :: obj
      integer(int32),                                 intent(in)  :: i_drd
      type(dipqcfs_t),                                intent(in)  :: dipq_c
      obj%DD2_ln_dipq_new(:) = 0.0_dp
      obj%DD2_ln_dipq_new(1:3) = obj%DD2_ln_dipq( (4*i_drd-3):(4*i_drd-1) ) - &
         2.0_dp*obj%dv( obj%i_b:obj%i_e )
      obj%DD2_ln_dipq_new(4)   = - 2.0_dp*( dipq_c%mat_A( (obj%i_b + 0), (obj%i_b + 0) ) + &
         dipq_c%mat_A( (obj%i_b + 1), (obj%i_b + 1) ) + &
         dipq_c%mat_A( (obj%i_b + 2), (obj%i_b + 2) ) )
   end subroutine new_DD2_dipq_fun
   subroutine upd_DD2_dipq_fun( obj, n_qdo )
      !use quantum_monte_carlo_v, only: smp_type
      class(dipqfun_t), intent(inout) :: obj
      integer(int32),                                 intent(in)    :: n_qdo
      integer(int32) :: i1
      do i1 = 1_int32, n_qdo
         obj%DD2_ln_dipq( (4*i1-3):(4*i1-1) ) = obj%DD2_ln_dipq( (4*i1-3):(4*i1-1) ) &
            -2.0_dp*obj%dv( (3*i1-2):(3*i1)  )
      enddo ! (n_qdo)
   end subroutine upd_DD2_dipq_fun
   subroutine cmp_dl_dipq_fun( obj, n_qdo, n_par_drd, dl_ln_prd )
      class(dipqfun_t), intent(inout) :: obj
      integer(int32),                           intent(in)    :: n_qdo
      integer(int32),                           intent(in)    :: n_par_drd
      real(dp),        dimension(n_par_drd), intent(inout)  :: dl_ln_prd
      integer(int32) :: i1, i2, i3
      i3 = 0_int32
      do i1 = 1_int32, 3_int32*n_qdo
         do i2 = i1, 3_int32*n_qdo
            i3 = i3 + 1_int32
            dl_ln_prd(i3) = -obj%x(i1)*obj%x(i2)
         enddo
      enddo
      dl_ln_prd(i3+1:i3+3*n_qdo) = 2.0_dp*obj%v(:)
   end subroutine cmp_dl_dipq_fun
end module dipqfun_cls
