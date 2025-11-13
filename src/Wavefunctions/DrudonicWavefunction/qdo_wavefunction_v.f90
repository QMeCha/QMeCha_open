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
module qdo_wavefunction_v
   use fortran_kinds_v, only: int8, int32
   implicit none
   integer(int32),  public, save   :: wf_type_d
   logical, save, public          :: ql_opt
   integer(int32), public, save   :: n_par_drd
   integer(int32), public, save   :: n_par_drd_s
   integer(int32), public, save   :: n_par_drd_wvfn
   integer(int32), public, save   :: n_par_drd_wvfn_s
end module qdo_wavefunction_v
