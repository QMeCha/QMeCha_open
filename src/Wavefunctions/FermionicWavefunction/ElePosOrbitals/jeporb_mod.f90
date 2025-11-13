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
module jeporb_mod
   use fortran_kinds_v,   only: int32, dp, stdout
   use openmp_mpi_m, only: mpi_rank
   use jastrow_factors_v,   only: n_par_jst, n_par_jst_s
   use write_lines_m
   use eporbbss_cls, only: epsysorb_t
   use eporbmat_cls, only: epmatorb_t
   use orbbss_cls,   only: orbspa_t
   use atomic_orbital_c,   only: orb_t
   implicit none
   integer(int32), save, public,               pointer :: n_jepos
   integer(int32), save, public,               pointer :: n_par_jepos
   integer(int32), save, public,               pointer :: n_par_jepos_s
   type(orb_t),    save, public, dimension(:), pointer :: jepos 
   type(epsysorb_t), save, public,                            target :: jepos_s
   type(epmatorb_t), save, public, allocatable, dimension(:), target :: jepos_mat
   logical, save, public :: jpoc_opt, jpoe_opt
   public  :: ini_jepos_mat
contains
   subroutine ini_jepos_mat( )
      use quantum_monte_carlo_v, only: n_wlk_max
      use basisset_v, only: jst_bsst
      integer(int32)  :: i1
      call jepos_s%ini( jst_bsst, jpoc_opt, jpoc_opt )
      jepos          => jepos_s%orb
      n_jepos        => jepos_s%n_orbs
      n_par_jepos   => jepos_s%n_par_orbs
      n_par_jepos_s => jepos_s%n_par_orbs_s
      if (n_par_jepos_s.eq.0_int32) then
         jpoc_opt = .false.
         jpoc_opt = .false.
      endif
      if ( n_jepos.gt.0_int32 ) then
         allocate( jepos_mat(1:n_wlk_max) )
         do i1 = 1_int32, n_wlk_max
            call jepos_mat(i1)%ini( i1, jepos_s )
         enddo
      endif
      call write_variable_line(stdout,0,mpi_rank,2,"Total number of positronic orbitals for Jastrow part", n_jepos)
   end subroutine ini_jepos_mat
end module jeporb_mod
