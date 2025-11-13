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
module feporb_mod
   use fortran_kinds_v, only: int32, dp, stdout
   use openmp_mpi_m, only: mpi_rank
   use fermionic_wavefunction_v, only: n_par_frm, n_par_frm_s
   use write_lines_m
   use eporbbss_cls, only: epsysorb_t
   use eporbmat_cls, only: epmatorb_t
   use orbbss_cls, only: orbspa_t
   use atomic_orbital_c, only: orb_t
   implicit none
   integer(int32), save, public,               pointer :: n_eporbs
   integer(int32), save, public,               pointer :: n_par_eporbs
   integer(int32), save, public,               pointer :: n_par_eporbs_s
   type(orb_t),    save, public, dimension(:), pointer :: eporbs !?
   type(epsysorb_t), save, public,                            target :: eporbs_s
   type(epmatorb_t), save, public, allocatable, dimension(:), target :: eporbs_mat
   logical, save, public :: poc_opt, poe_opt
   public  :: ini_eporbs_mat
contains
   subroutine ini_eporbs_mat( )
      use quantum_monte_carlo_v, only: n_wlk_max
      use basisset_v, only: frm_bsst
      integer(int32)  :: i1
      !? Allocating orbitals per atom
      call eporbs_s%ini( frm_bsst, poc_opt, poe_opt )
      eporbs          => eporbs_s%orb
      n_eporbs        => eporbs_s%n_orbs
      n_par_eporbs   => eporbs_s%n_par_orbs
      n_par_eporbs_s => eporbs_s%n_par_orbs_s
      if (n_par_eporbs_s.eq.0_int32) then
         poc_opt = .false.
         poe_opt = .false.
      endif
      n_par_frm   = n_par_frm + n_par_eporbs
      n_par_frm_s = n_par_frm_s + n_par_eporbs_s
      if ( n_eporbs.gt.0_int32 ) then
         allocate( eporbs_mat(1:n_wlk_max) )
         do i1 = 1_int32, n_wlk_max
            call eporbs_mat(i1)%ini( i1, eporbs_s )
         enddo 
      endif
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of positronic orbitals for fermionic part", n_eporbs)
   end subroutine ini_eporbs_mat
end module feporb_mod
