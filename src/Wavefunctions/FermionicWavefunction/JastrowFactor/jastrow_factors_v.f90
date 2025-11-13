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
module jastrow_factors_v
   use fortran_kinds_v, only: int32, dp
   implicit none
   logical, save, public          :: jst_prs
   logical, save, public          :: jc2_prs
   logical, save, public          :: jd1_prs, ejd1_prs, pjd1_prs
   logical, save, public          :: jd2_prs, ejd2_prs, pjd2_prs
   logical, save, public          :: jc1_opt, ejc1_opt, pjc1_opt
   logical, save, public          :: jc2_opt, ejc2_opt, pjc2_opt
   logical, save, public          :: jce_opt
   logical, save, public          :: jdc_opt
   logical, save, public          :: jd_opt,  ejd_opt,  pjd_opt
   logical, save, public          :: jd1_opt, ejd1_opt, pjd1_opt
   logical, save, public          :: jd2_opt, ejd2_opt, pjd2_opt
   integer(int32), save, public   :: n_par_jstc
   integer(int32), save, public   :: n_par_jstd, n_par_jstd_s
   integer(int32), save, public   :: n_par_jepd, n_par_jepd_s
   integer(int32), save, public   :: n_par_jepa, n_par_jepa_s
   integer(int32), save, public   :: n_par_jst,  n_par_jst_s
   logical, save, public          :: jdep_opt
   logical, save, public          :: jdep_prs, jdte_prs, jdtp_prs, jqep_prs
   logical, save, public :: jaep_opt, jaep_prs
   logical, save, public :: eTepa_prs, pTepa_prs, eQepa_prs, pQepa_prs
end module jastrow_factors_v
