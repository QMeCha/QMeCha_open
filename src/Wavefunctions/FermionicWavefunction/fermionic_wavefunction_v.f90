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
module fermionic_wavefunction_v
   use fortran_kinds_v, only: int32
   implicit none
   character(1),   public, save :: spin
   integer(int32), public, save :: wf_type_e, wf_type_p, wf_type_ep
   logical,        public, save :: dl_opt, edl_opt, pdl_opt
   integer(int32), public, save :: n_par_det_e, n_par_det_e_s
   integer(int32), public, save :: n_par_det_p, n_par_det_p_s
   integer(int32), public, save :: n_par_det_ep, n_par_det_ep_s
   integer(int32), public, save :: n_par_det,   n_par_det_s
   integer(int32), public, save :: n_par_frm,   n_par_frm_s
   integer(int32), public, save :: n_par_frm_wvfn
   integer(int32), public, save :: n_par_frm_wvfn_s
end module fermionic_wavefunction_v
