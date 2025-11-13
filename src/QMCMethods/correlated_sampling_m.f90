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
module correlated_sampling_m
   use fortran_kinds_v,          only: dp, int32, stdout
   use openmp_mpi_m,             only: mpi_rank
   use write_lines_m
   use molecular_system_v,       only: n_fe, molsys_prs
   use qdo_system_v,             only: n_qdo, qdosys_prs
   use fermionic_config_c,       only: frm_cnf
   use drudonic_config_c,        only: drd_cnf
   use fermionic_wavefunction_m, only: wvfn
   use qdo_wavefunction_m,       only: wvfnq
   implicit none
   type, public :: crrsmp_t
      real(dp), allocatable, dimension(:)     :: ratio2
      real(dp), allocatable, dimension(:)     :: spsi
      real(dp), allocatable, dimension(:)     :: logpsi
      real(dp), allocatable, dimension(:)     :: e_b
      real(dp), allocatable, dimension(:,:,:) :: r_fe_sav
      real(dp), allocatable, dimension(:,:,:) :: r_drd_sav
   end type crrsmp_t
   real(dp), public, save :: W_avg, W_err, W2_avg, Ovrlap_avg, logshift
   type(crrsmp_t), public, save, allocatable, dimension(:) :: crrsmp
   public :: Init_CorrelatedSampling, Save_CorrelatedSampling, Fnlz_CorrelatedSampling
contains
   subroutine Init_CorrelatedSampling( bin_l, n_bin, n_wlk )
      integer(int32), intent(in) :: bin_l, n_bin, n_wlk
      integer(int32) :: iw
      allocate( crrsmp(1:n_wlk) )
      do iw = 1_int32, n_wlk
         allocate( crrsmp(iw)%logpsi(1:bin_l*n_bin) )
         allocate( crrsmp(iw)%spsi(1:bin_l*n_bin) )
         allocate( crrsmp(iw)%ratio2(1:bin_l*n_bin) )
         if (molsys_prs) allocate( crrsmp(iw)%r_fe_sav(1:3,1:n_fe,1:bin_l*n_bin) )
         if (qdosys_prs) allocate( crrsmp(iw)%r_drd_sav(1:3,1:n_qdo,1:bin_l*n_bin) )
      enddo
      W_avg = 0.0_dp
      Ovrlap_avg = 1.0_dp
      W2_avg = 0.0_dp
      W_err = 0.0_dp
   end subroutine Init_CorrelatedSampling
   subroutine Save_CorrelatedSampling( iw, is )
      integer(int32), intent(in) :: iw, is
      crrsmp(iw)%spsi(is)   = 1.0_dp
      crrsmp(iw)%logpsi(is) = 0.0_dp
      if (molsys_prs) then
         crrsmp(iw)%spsi(is)   = crrsmp(iw)%spsi(is)*wvfn%swvfn(iw)
         crrsmp(iw)%logpsi(is) = crrsmp(iw)%logpsi(is) + wvfn%logwvfn(iw)
         crrsmp(iw)%r_fe_sav(1:3,1:n_fe,is) = frm_cnf(iw)%r_fe(1:3,1:n_fe)
      endif
      if (qdosys_prs) then
         crrsmp(iw)%spsi(is)   = crrsmp(iw)%spsi(is)*wvfnq(iw)%swvfnq
         crrsmp(iw)%logpsi(is) = crrsmp(iw)%logpsi(is) + wvfnq(iw)%logwvfnq
         crrsmp(iw)%r_drd_sav(1:3,1:n_qdo,is) = drd_cnf(iw)%r_drd(1:3,1:n_qdo)
      endif
   end subroutine Save_CorrelatedSampling
   subroutine Fnlz_CorrelatedSampling( n_wlk )
      integer(int32), intent(in) :: n_wlk
      integer(int32) :: iw
      do iw = 1_int32, n_wlk
         deallocate( crrsmp(iw)%logpsi )
         deallocate( crrsmp(iw)%spsi )
         deallocate( crrsmp(iw)%ratio2 )
         if (molsys_prs) deallocate( crrsmp(iw)%r_fe_sav )
         if (qdosys_prs) deallocate( crrsmp(iw)%r_drd_sav )
      enddo
      deallocate( crrsmp )
   end subroutine Fnlz_CorrelatedSampling
end module correlated_sampling_m
