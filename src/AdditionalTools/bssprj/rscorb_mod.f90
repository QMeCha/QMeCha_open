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
module rscorb_mod
use fortran_kinds_v, only: dp, int32
use atomic_orbital_c, only: orb_t
implicit none
public :: cmp_orb, cmp_dpar_orb
contains
subroutine cmp_orb( orb_par, r, P_n, orb )
  use primorbs_operations_m, only: cmp_prim_orb
  type(orb_t),                      intent(inout) :: orb_par
  real(dp),                         intent(in)    :: r
  real(dp), dimension(orb_par%n), intent(inout) :: P_n
  real(dp),                         intent(inout) :: orb
  integer(int32) :: i1
  do i1 = 1_int32, orb_par%n
    P_n(i1) = cmp_prim_orb( r, orb_par%z(i1), orb_par%l, orb_par%prm_n(i1), orb_par%prm_t(i1) )
  enddo
  orb = dot_product( orb_par%c(:), P_n(:) )
  orb = orb * r ** orb_par%l
end subroutine cmp_orb
subroutine cmp_dpar_orb( orb_par, r, n_par_orbs, de_opt, P_n, da_orb )
  use primorbs_operations_m, only: cmp_Z_fact
  type(orb_t),                       intent(inout) :: orb_par
  real(dp),                          intent(in)    :: r
  integer(int32),                    intent(in)    :: n_par_orbs
  logical,                           intent(in)    :: de_opt
  real(dp), dimension(orb_par%n),  intent(in)    :: P_n
  real(dp), dimension(n_par_orbs), intent(inout) :: da_orb
  real(dp) :: Z_n
  integer(int32) :: i1
  da_orb(1:orb_par%n) = P_n(1:orb_par%n) 
  if( de_opt ) then
    do i1 = 1_int32, orb_par%n 
      call cmp_Z_fact( r, orb_par%prm_n(i1), orb_par%l, orb_par%z(i1), orb_par%prm_t(i1), Z_n )
      da_orb(orb_par%n+i1) = orb_par%c(i1) * Z_n * P_n(i1)
    enddo   
  endif 
  da_orb = da_orb * r ** orb_par%l
end subroutine cmp_dpar_orb
end module rscorb_mod
