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
module trgtgo_mod
use fortran_kinds_v, only: dp, int32
use intgrd_mod, only: n_int_points, orb_cut_off, dorb, trgt_grid, nuc_cut_off
use cusps_functions_params_m, only: b_en, jcen_prs, n_par_jcen
implicit none
integer(int32), save, public                              :: n_tot_shlls
real(dp)      , save, public, allocatable, dimension(:,:) :: trgt_orb_vals
public  :: ini_trgt_orbs
!private ::
contains
subroutine ini_trgt_orbs( )
  use basisset_v, only: n_fbs, frm_bsst
  use orbmat_cls, only: prim_po_t
  use rscorb_mod, only: cmp_orb
  use fermions_nuclei_cusps_m, only: init_fermions_nuclei_cusps
  type(prim_po_t), allocatable, dimension(:)   :: prim_po
  real(dp)                                     :: orb
  integer(int32)          :: Z_a
  integer(int32)          :: i1, i2, i3, i4
  integer(int32), pointer :: l
  integer(int32), external :: comp_atomic_charge
  call init_fermions_nuclei_cusps( b_en, jcen_prs, n_par_jcen, .false., .false. )
  n_tot_shlls = sum(frm_bsst(1:n_fbs)%n_shlls)
  allocate( prim_po(1:n_tot_shlls) ) 
  i3 = 1_int32
  do i1 = 1_int32, n_fbs ; do i2 = 1_int32, frm_bsst(i1)%n_shlls
    allocate( prim_po(i3)%p(1:frm_bsst(i1)%shlls(i2)%n) ) ; prim_po(i3)%p = 0.0_dp 
    i3 = i3 + 1_int32
  enddo ; enddo ! i1 (n_fbs)
  allocate(trgt_grid(1:n_int_points,1:n_tot_shlls)) ; trgt_grid(1:n_int_points,1:n_tot_shlls) = 0.0_dp
  allocate(trgt_orb_vals(1:n_int_points,1:n_tot_shlls)) ; trgt_orb_vals(1:n_int_points,1:n_tot_shlls) = 0.0_dp
  i3 = 1_int32
  do i1 = 1_int32, n_fbs 
    if ( frm_bsst(i1)%atm_name(1:1) .eq. '*') then
      Z_a = comp_atomic_charge( frm_bsst(i1)%atm_name(2:3) )
    else
      Z_a = comp_atomic_charge( frm_bsst(i1)%atm_name(1:2) )
    endif
    do i2 = 1_int32, frm_bsst(i1)%n_shlls
    l => frm_bsst(i1)%shlls(i2)%l
    dorb = 0.0
    orb = 1.0_dp
    do while (abs(orb).gt.orb_cut_off)
      dorb = dorb + 1.0_dp / dble(n_int_points)
      !if ( l.eq.0_int32 ) then
      !  trgt_grid(1,i3) = 0.01_dp/ dble(Z_a) !01d0
      !else
        trgt_grid(1,i3) = nuc_cut_off
      !endif
      do i4 = 2_int32, n_int_points
        trgt_grid(i4,i3) = trgt_grid(i4-1,i3) + dorb !* (1.0-exp(-trgt_grid(i4-1,i3)))
      enddo
      call cmp_orb( frm_bsst(i1)%shlls(i2), trgt_grid(n_int_points,i3), prim_po(i3)%p, orb )
      if ( jcen_prs ) then
        orb = orb * exp(f_fncfun( 1, b_en, trgt_grid(n_int_points,i3) ))
      endif
    enddo
    !write(10+i2,*)  
    do i4 = 1_int32, n_int_points
      call cmp_orb( frm_bsst(i1)%shlls(i2), trgt_grid(i4,i3), prim_po(i3)%p, trgt_orb_vals(i4,i3) )
      if ( jcen_prs ) then
        trgt_orb_vals(i4,i3) = trgt_orb_vals(i4,i3) * exp(f_fncfun( 1, b_en, trgt_grid(i4,i3) ) )
      endif
      !write(10+i2,*) trgt_grid(i4,i3), trgt_orb_vals(i4,i3)
    enddo
    i3 = i3 + 1_int32
  enddo ; enddo ! i1 (n_at) i2 (obj%orbs(i2)%n_orbs)
end subroutine ini_trgt_orbs
end module trgtgo_mod
