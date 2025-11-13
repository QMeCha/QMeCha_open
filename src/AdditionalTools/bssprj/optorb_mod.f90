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
module optorb_mod
use fortran_kinds_v, only: dp, int32
use intgrd_mod, only: n_int_points, trgt_grid
use trgtgo_mod, only: trgt_orb_vals, n_tot_shlls
use atomic_orbital_c, only: orb_t
use cusps_functions_params_m, only: b_en, jcen_prs, n_par_jcen
implicit none
logical, public, save :: de_opt
integer(int32), public, save :: n_par_orbs
integer(int32), save, public :: Z_a
type(orb_t),    save, public, pointer :: shll
integer(int32), save, public          :: i_shll
real(dp),       save, public, allocatable, dimension(:,:) :: P_n
real(dp),       save, public, allocatable, dimension(:)   :: vec_orb_par
public :: ufparm, run_opt_orbs, orbs_res, orbs_jcb
contains
subroutine ini_opt_orbs
  use fermions_nuclei_cusps_m, only: init_fermions_nuclei_cusps
  use datpio_mod, only: oe_opt, b_par, l_par
  if ( allocated(b_en%b) ) deallocate(b_en%b)
  call init_fermions_nuclei_cusps( b_en, jcen_prs, n_par_jcen, .false., .false. )
  if(jcen_prs) b_en%b(1,:) = b_par 
  de_opt = oe_opt
end subroutine ini_opt_orbs
subroutine run_opt_orbs
  use basisset_v, only: n_fbs, frm_bsst
  use rscorb_mod, only: cmp_orb
  real(dp):: orb
  integer(int32), allocatable, dimension(:) :: iv
  real(dp),       allocatable, dimension(:) :: v
  integer(int32)                :: uiparm(1)
  real(dp)                      :: urparm(1)
  integer(int32), external :: comp_atomic_charge
  integer(int32) :: i1, i2, i3
  i_shll = 1_int32
  do i1 = 1_int32, n_fbs 
    if ( frm_bsst(i1)%atm_name(1:1) .eq. '*') then
      Z_a = comp_atomic_charge( frm_bsst(i1)%atm_name(2:3) )
    else
      Z_a = comp_atomic_charge( frm_bsst(i1)%atm_name(1:2) )
    endif
    do i2 = 1_int32, frm_bsst(i1)%n_shlls
    shll => frm_bsst(i1)%shlls(i2)
    allocate( P_n(1:shll%n,1:n_int_points) ) ; P_n = 0.0_dp
    if ( de_opt ) then
      n_par_orbs = 2_int32 * shll%n
    else
      n_par_orbs = shll%n
    endif
    allocate(  vec_orb_par(1:n_par_orbs) )
    vec_orb_par ( 1:shll%n ) = shll%c
    if (de_opt) then
      vec_orb_par ( shll%n+1:n_par_orbs ) = shll%z
    endif
    allocate( iv (1:60+n_par_orbs) ) ; iv = 0_int32
    allocate( v (1:93+n_int_points*n_par_orbs+3*n_int_points+n_par_orbs*(3*n_par_orbs+33)/2) ) ; v = 0.0_dp
    call nl2sol ( n_int_points, n_par_orbs, vec_orb_par, orbs_res, orbs_jcb, iv, v, uiparm, urparm, ufparm )
    do i3 = 1_int32, n_int_points
      call cmp_orb( shll, trgt_grid(i3,i_shll), P_n, orb)
      if ( jcen_prs ) then
        orb = orb * exp( f_fncfun( 1, b_en, trgt_grid(i3,i_shll) ) )
      endif  
    enddo
    deallocate(vec_orb_par,P_n, iv, v)
    i_shll = i_shll + 1_int32
  enddo ; enddo ! i1 (n_fbs) i2 (frm_bsst(i1)%n_shlls)
end subroutine run_opt_orbs
subroutine orbs_res ( n_grid_p, n_par, vec_par, nf, vec_res, uiparm, urparm, ufdummy )
  use rscorb_mod, only: cmp_orb
  integer(int32),                  intent(inout) :: n_grid_p
  integer(int32),                  intent(inout) :: n_par
  real(dp), dimension(n_par)   , intent(inout) :: vec_par
  integer(int32),                  intent(inout) :: nf
  real(dp), dimension(n_grid_p), intent(inout) :: vec_res
  external                      :: ufdummy
  integer(int32)                :: uiparm(*)
  real(dp)                      :: urparm(*)
  integer(int32) :: i1
  vec_res = 0.0_dp
  shll%c = vec_par ( 1:shll%n )
  if (de_opt) then
    do i1 = 1_int32, shll%n 
      if (vec_par(shll%n+i1).gt.0_int32) then
        shll%z(i1) = vec_par(shll%n+i1)
      else
        vec_par(shll%n+i1) = shll%z(i1)
      endif
    enddo
  endif
  do i1 = 1_int32, n_grid_p
    call cmp_orb( shll, trgt_grid(i1,i_shll), P_n(:,i1), vec_res(i1) )
    if ( jcen_prs ) then
      vec_res(i1) = vec_res(i1) * exp( f_fncfun( 1, b_en, trgt_grid(i1,i_shll) ) )
    endif
  enddo
  vec_res = vec_res - trgt_orb_vals(:,i_shll)
end subroutine orbs_res
subroutine orbs_jcb ( n_grid_p, n_par, vec_par, nf, mat_jcb, uiparm, urparm, ufdummy )
  use rscorb_mod, only: cmp_dpar_orb
  integer(int32),                  intent(inout) :: n_grid_p
  integer(int32),                  intent(inout) :: n_par
  real(dp), dimension(n_par)   , intent(inout) :: vec_par
  integer(int32),                  intent(inout) :: nf
  real(dp), dimension(n_grid_p, 1:n_par), intent(inout) :: mat_jcb
  external                      :: ufdummy
  integer(int32)                :: uiparm(*)
  real(dp)                      :: urparm(*)
  integer(int32) :: i1
  mat_jcb = 0.0_dp
  shll%c = vec_par ( 1:shll%n )
  if (de_opt) then
    do i1 = 1_int32, shll%n 
      if (vec_par(shll%n+i1).gt.0_int32) then
        shll%z(i1) = vec_par(shll%n+i1)
      else
        vec_par(shll%n+i1) = shll%z(i1)
      endif
    enddo
  endif
  do i1 = 1_int32, n_grid_p
    call cmp_dpar_orb( shll, trgt_grid(i1,i_shll), n_par, de_opt, P_n(:,i1), mat_jcb(i1,:) )
    if ( jcen_prs ) then
      mat_jcb(i1,:) = mat_jcb(i1,:) * exp( f_fncfun( 1, b_en, trgt_grid(i1,i_shll) ) )
    endif
  enddo
end subroutine orbs_jcb
subroutine ufparm (  )
end subroutine ufparm
subroutine prnt_opt_orbs
  use basisset_v, only: n_fbs, frm_bsst
  use rscorb_mod, only: cmp_orb
  real(dp) :: orb
  real(dp) :: w(1:n_tot_shlls)
  integer(int32), external :: comp_atomic_charge
  integer(int32) :: i1, i2, i3
  i_shll = 1_int32
  w = 0.0_dp
  do i1 = 1_int32, n_fbs 
    if ( frm_bsst(i1)%atm_name(1:1) .eq. '*') then
      Z_a = comp_atomic_charge( frm_bsst(i1)%atm_name(2:3) )
    else
      Z_a = comp_atomic_charge( frm_bsst(i1)%atm_name(1:2) )
    endif
    do i2 = 1_int32, frm_bsst(i1)%n_shlls
      shll => frm_bsst(i1)%shlls(i2)
      allocate( P_n(1:shll%n,1:n_int_points) ) ; P_n = 0.0_dp
      do i3 = 1_int32, n_int_points
        call cmp_orb( shll, trgt_grid(i3,i_shll), P_n, orb)
        if ( jcen_prs ) then
          orb = orb * exp(f_fncfun( 1, b_en, trgt_grid(i3,i_shll) ))
        endif  
        write(10+i_shll,'(4E16.7)') trgt_grid(i3,i_shll), trgt_orb_vals(i3,i_shll), orb, orb - trgt_orb_vals(i3,i_shll)
        w(i_shll) = w(i_shll) + (orb - trgt_orb_vals(i3,i_shll))**2
      enddo
      deallocate(P_n)
      i_shll = i_shll + 1_int32
  enddo ; enddo ! i1 (n_fbs) i2 (frm_bsst(i1)%n_shlls)
  write(*,'("W :",200E12.4)') sqrt(w(1:n_tot_shlls)), sqrt(sum(w(1:n_tot_shlls)))
end subroutine prnt_opt_orbs
end module optorb_mod
