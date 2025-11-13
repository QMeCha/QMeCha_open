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
module kinetic_energy_m
   use fortran_kinds_v,          only: int32, dp
   use molecular_system_v,       only: n_fe
   use fermionic_wavefunction_m, only: wvfn
   use qdo_system_v,             only: n_qdo
   use qdo_wavefunction_m,       only: wvfnq
   implicit none
   public :: cmp_kinene, upd_kinene, upd_kinene_emb
contains
   subroutine cmp_kinene( iw )
      integer(int32), intent(in) :: iw
      if ( n_fe.gt.0_int32 ) then
         call wvfn%cmp_D( iw )
      endif
      if ( n_qdo.gt.0_int32 ) then
         call wvfnq(iw)%cmp_D( iw )
      endif
   end subroutine cmp_kinene
   subroutine upd_kinene( iw )
      integer(int32), intent(in) :: iw
      if ( n_fe.gt.0_int32 ) then
         call wvfn%upd_D( iw )
      endif
      if ( n_qdo.gt.0_int32 ) then
         call wvfnq(iw)%upd_D( iw )
      endif
   end subroutine upd_kinene
   subroutine upd_kinene_emb( iw )
      integer(int32), intent(in) :: iw
      if ( n_qdo.gt.0_int32 ) then
         call wvfnq(iw)%upd_D( iw )
      endif
   end subroutine upd_kinene_emb
end module kinetic_energy_m
