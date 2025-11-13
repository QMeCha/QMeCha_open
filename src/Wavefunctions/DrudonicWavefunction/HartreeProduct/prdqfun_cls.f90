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
module prdqfun_cls
use fortran_kinds_v,  only: dp, int32 
use qdo_wavefunction_v,  only: ql_opt
use drudonic_orbitals_m,  only: n_dorbs, n_par_dorbs, n_par_dorbs_s, dorbs
use prdqcfs_mod, only: prdqcfs_t 
use quantum_monte_carlo_v, only: grad_updat
implicit none
type, public :: prdqfun_t
  real(dp)                               :: psi_prd
  real(dp)                               :: g
  real(dp), allocatable, dimension(:)    :: qso
  real(dp)                               :: qso_new
  real(dp), allocatable, dimension(:)    :: DD2_ln_prdq
  real(dp), allocatable, dimension(:)    :: DD2_ln_prdq_new
  contains
  procedure :: ini    => ini_prdq_fun
  procedure :: cmp    => cmp_prdq_fun
  procedure :: ratio  => ratio_prdq_fun
  procedure :: upd    => upd_prdq_fun
  procedure :: cmp_D  => cmp_DD2_prdq_fun
  procedure :: new_D  => new_DD2_prdq_fun
  procedure :: upd_D  => upd_DD2_prdq_fun
  procedure :: cmp_dl => cmp_dl_prdq_fun
  procedure :: cmp_da => cmp_da_prdq_fun
end type prdqfun_t
private :: ini_prdq_fun, cmp_prdq_fun, ratio_prdq_fun, upd_prdq_fun, &
         & cmp_DD2_prdq_fun, new_DD2_prdq_fun, upd_DD2_prdq_fun, &
         & cmp_dl_prdq_fun, cmp_da_prdq_fun
contains
subroutine ini_prdq_fun( obj, n_qdo )
  class(prdqfun_t), intent(inout) :: obj
  integer(int32),  intent(in)    :: n_qdo
  if (n_qdo.gt.0_int32) then 
    allocate( obj%qso(1:n_qdo) )                 ; obj%qso             = 0.0_dp
    allocate( obj%DD2_ln_prdq(1:4*n_qdo) )       ; obj%DD2_ln_prdq     = 0.0_dp
    if (grad_updat) then
      allocate( obj%DD2_ln_prdq_new(1:4) )       ; obj%DD2_ln_prdq_new = 0.0_dp
    endif 
  endif 
end subroutine ini_prdq_fun
subroutine cmp_prdq_fun( obj, n_qdo, prdq_c, o )
  class(prdqfun_t), intent(inout) :: obj
  integer(int32),                             intent(in) :: n_qdo
  type(prdqcfs_t),                            intent(in) :: prdq_c
  real(dp),  dimension(n_dorbs,n_qdo),    intent(in) :: o
  integer(int32) :: i1
  do i1 = 1_int32, n_qdo
    obj%qso(i1) = dot_product( prdq_c%mat_L(1:n_dorbs,i1), o(1:n_dorbs,i1) )
  enddo ! i1 (n_qdo)
  obj%psi_prd = product(obj%qso(:))
end subroutine cmp_prdq_fun
subroutine ratio_prdq_fun( obj, i_drd, prdq_c, o_new, g )
  class(prdqfun_t), intent(inout) :: obj
  integer(int32),                             intent(in)  :: i_drd
  type(prdqcfs_t),                            intent(in)  :: prdq_c
  real(dp),      dimension(n_dorbs),        intent(in)  :: o_new    
  real(dp),                                   intent(out) :: g   
  obj%qso_new = dot_product( prdq_c%mat_L(1:n_dorbs,i_drd), o_new )
  obj%g = obj%qso_new / obj%qso(i_drd)
  g = obj%g 
end subroutine ratio_prdq_fun
subroutine upd_prdq_fun( obj, i_drd )
  class(prdqfun_t), intent(inout) :: obj
  integer(int32),                             intent(in)  :: i_drd
  obj%psi_prd = obj%psi_prd * obj%g
  obj%qso(i_drd) = obj%qso_new
end subroutine upd_prdq_fun
subroutine cmp_DD2_prdq_fun( obj, n_qdo, prdq_c, DD2_o )
  class(prdqfun_t), intent(inout) :: obj
  integer(int32),                             intent(in) :: n_qdo
  type(prdqcfs_t),                            intent(in) :: prdq_c
  real(dp),  dimension(n_dorbs,4*n_qdo),  intent(in) :: DD2_o
  integer(int32) :: i1
  do i1 = 1_int32, n_qdo
    call dgemv('T', n_dorbs, 4, 1.0_dp/obj%qso(i1), DD2_o(:,4*i1-3:4*i1), n_dorbs, &
              & prdq_c%mat_L(:,i1), 1_int32, 0.0_dp, obj%DD2_ln_prdq(4*i1-3:4*i1), 1_int32)
    obj%DD2_ln_prdq(4*i1) = obj%DD2_ln_prdq(4*i1) - sum( obj%DD2_ln_prdq(4*i1-3:4*i1-1)**2 )
  enddo ! i1 (n_qdo)
end subroutine cmp_DD2_prdq_fun
subroutine new_DD2_prdq_fun( obj, i_drd, prdq_c, DD2_o_new )
  class(prdqfun_t), intent(inout) :: obj
  integer(int32),                           intent(in)  :: i_drd
  type(prdqcfs_t),                          intent(in)  :: prdq_c
  real(dp),      dimension(n_dorbs,4),  intent(in)  :: DD2_o_new  
  call dgemv('T', n_dorbs, 4, 1.0_dp/obj%qso_new, DD2_o_new(:,4), n_dorbs, &
            & prdq_c%mat_L(:,i_drd), 1_int32, 0.0_dp, obj%DD2_ln_prdq_new(:) , 1_int32)
  obj%DD2_ln_prdq_new(4) = obj%DD2_ln_prdq_new(4)  - sum( obj%DD2_ln_prdq_new(1:3)**2 )
end subroutine new_DD2_prdq_fun
subroutine upd_DD2_prdq_fun( obj, n_qdo, i_drd, prdq_c, DD2_o )
  class(prdqfun_t), intent(inout) :: obj
  !? we might compute DD2_o_new also for grad_updat, send it here
  !? and evaluate DD2_o_new and then add it into DD2_o
  integer(int32),                                 intent(in) :: n_qdo
  integer(int32),                                 intent(in) :: i_drd
  type(prdqcfs_t),                                intent(in) :: prdq_c
  real(dp),      dimension(n_dorbs,4*n_qdo),  intent(in) :: DD2_o
  integer(int32) :: i1
  do i1 = 1_int32, n_qdo
    if (i1.eq.i_drd) then
      obj%DD2_ln_prdq(4*i1-3:4*i1) = obj%DD2_ln_prdq_new(1:4)
    else
      call dgemv('T', n_dorbs, 4, 1.0_dp/obj%qso(i1), DD2_o(:,4*i1-3:4*i1), n_dorbs, &
                & prdq_c%mat_L(:,i1), 1_int32, 0.0_dp, obj%DD2_ln_prdq(4*i1-3:4*i1), 1_int32)
      obj%DD2_ln_prdq(4*i1) = obj%DD2_ln_prdq(4*i1) - sum( obj%DD2_ln_prdq(4*i1-3:4*i1-1)**2 )
    endif
  enddo ! i1 (n_qdo)
end subroutine upd_DD2_prdq_fun
subroutine cmp_dl_prdq_fun( obj, n_qdo, o, n_par_drd, dl_ln_prd )
  class(prdqfun_t), intent(inout) :: obj
  integer(int32),                                intent(in)    :: n_qdo
  real(dp),        dimension(n_dorbs,n_qdo), intent(in)    :: o
  integer(int32),                                intent(in)    :: n_par_drd
  real(dp),        dimension(n_par_drd),       intent(inout) :: dl_ln_prd
  integer(int32) :: i1, ip
  ip = 0.0_dp
  do i1 = 1_int32, n_qdo
    dl_ln_prd(ip+1:ip+n_dorbs) = o(1:n_dorbs,i1) / obj%qso(i1)
    ip = ip + n_dorbs
  enddo ! i1 (n_qdo)
end subroutine cmp_dl_prdq_fun
subroutine cmp_da_prdq_fun( obj, n_qdo, prdq_c, da_o, da_ln_prd )
  class(prdqfun_t), intent(inout) :: obj
  integer(int32),                             intent(in)    :: n_qdo
  type(prdqcfs_t),                            intent(in)    :: prdq_c
  real(dp), dimension(n_par_dorbs,n_qdo), intent(in)    :: da_o
  real(dp), dimension(n_par_dorbs),         intent(inout) :: da_ln_prd
  integer(int32) :: i1, i2, i3, io, pi, pf
  da_ln_prd = 0.0_dp
  do i1 = 1_int32, n_qdo
    do i2 = 1_int32, n_qdo ; do i3 = 1_int32, dorbs(i2)%n_orbs
      io = dorbs(i2)%orb(i3)%i_orb
      if( dorbs(i2)%orb(i3)%n_par.gt.0_int32 ) then
        pi = dorbs(i2)%orb(i3)%i_par
        pf = pi + dorbs(i2)%orb(i3)%n_par - 1_int32
        da_ln_prd(pi:pf) = da_ln_prd(pi:pf) + prdq_c%mat_L(io,i1) * da_o(pi:pf,i1) / obj%qso(i1) 
      endif
    enddo ; enddo  !  i2 (n_at) i3 (orbs(i2)%n_orbs)
  enddo ! i1 (n_el_u)
end subroutine cmp_da_prdq_fun
end module prdqfun_cls
