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
module monte_carlo_averages_v
   use fortran_kinds_v, only: int32, dp
   implicit none
   integer(int32), public, save, target :: n_obs, n_obs2, n_pp_obs, n_tot_obs
   real(dp),       public, save, allocatable, dimension(:,:), target :: obs_b
   integer(int32), public, save         :: min_blck_size, max_blck_size
   integer(int32), public, save         :: n_blck_sizes
   integer(int32), public, save, allocatable, dimension(:) :: blck_sizes
   integer(int32), public, save, allocatable, dimension(:) :: blck_l
   real(dp), public, save, allocatable, dimension(:), target :: obs_avg
   real(dp), public, save, allocatable, dimension(:), target :: obs_var
   real(dp), public, save, allocatable, dimension(:,:,:), target :: obs_pp_b
   real(dp), public, save, allocatable, dimension(:,:),   target :: obs_pp_avg
   real(dp), public, save, allocatable, dimension(:),     target :: obs_pp_var
   real(dp), public, save, allocatable, dimension(:,:)         :: blck_avg
   real(dp), public, save, allocatable, dimension(:,:)         :: blck_var
   real(dp), public, save, allocatable, dimension(:,:), target :: blck_err
   real(dp), pointer               :: w_avg, w2_avg
   real(dp), pointer               :: e_avg, e_var, k_avg, v_avg, e_blck_err, k_blck_err, v_blck_err
   real(dp), dimension(:), pointer :: e_pp_avg, e_pp_var, e_pp_blck_err
   real(dp), dimension(:), pointer :: dip_avg, dip_blck_err
   real(dp), dimension(:), pointer :: qud_avg, qud_blck_err
   real(dp), pointer               :: e_block_err, k_block_err, v_block_err
   real(dp), pointer               :: v_ee_avg, v_ee_blck_err    
   real(dp), dimension(:), pointer :: v_fn_avg, v_fn_blck_err
   real(dp), pointer               :: v_dd_avg, v_dd_blck_err    
   real(dp), dimension(:), pointer :: v_dq_avg, v_dq_blck_err
   real(dp), dimension(:), pointer :: r_fe_avg, r_fe_blck_err
   real(dp), dimension(:), pointer :: r_drd_avg, r_drd_blck_err
   real(dp), pointer :: d_ee_avg, d_ee_blck_err  
   real(dp), pointer :: d_pp_avg, d_pp_blck_err 
   real(dp), pointer :: d_ep_avg, d_ep_blck_err 
   real(dp), pointer :: d_dd_avg, d_dd_blck_err 
end module monte_carlo_averages_v
