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
module jstqepar_var
   use fortran_kinds_v, only: int32
   implicit none
   logical, save, public          :: jstqe_prs
   integer(int32), save, public   :: jstqe_type
   logical, save, public          :: jqe_opt
   integer(int32), save, public   :: n_par_jstqe
end module jstqepar_var
