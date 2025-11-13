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
module triplet_geminal_c
   use fortran_kinds_v,           only: dp, int32
   use triplet_geminal_params_m,  only: tgmcfs_t
   use fermionic_wavefunction_v,  only: spin, n_par_frm, n_par_frm_s
   use fermionic_orbitals_m,      only: n_forbs, n_par_forbs, aoc_opt
   use quantum_monte_carlo_v,     only: grad_updat, n_wlk_max
   use triplet_geminal_operats_m, only: cmp_Z_mat, cmp_T_mat, cmp_T_vec, upd_Z_mat,&
   & upd_T_inv, cmp_D_mat, cmp_DD2_vec, upd_D_mat, upd_DD2_vec, cmp_H_mat, cmp_da_vec
   implicit none
   type, public :: tgmfun_t
      integer(int32)                            :: d_Tu, d_Td
      real(dp),       allocatable, dimension(:) :: pff_T_u, pff_T_d, g
      integer(int32), allocatable, dimension(:) :: n_mat_chngs_u, n_mat_chngs_d
      real(dp),       allocatable, dimension(:,:,:) :: mat_T_u, mat_T_d
      real(dp),       allocatable, dimension(:,:,:) :: inv_T_u, inv_T_d
      real(dp),       allocatable, dimension(:,:,:) :: mat_Z
      real(dp),       allocatable, dimension(:,:,:) :: mat_D
      real(dp),       allocatable, dimension(:,:,:) :: mat_OTi
      real(dp),       allocatable, dimension(:,:,:) :: mat_H
      real(dp),       allocatable, dimension(:,:)   :: da_ln_pff_T_u, da_ln_pff_T_d
      real(dp),       allocatable, dimension(:,:)   :: vec_T
      real(dp),       allocatable, dimension(:,:)   :: vec_Tit
      real(dp),       allocatable, dimension(:,:)   :: vec_Ti
      real(dp),       allocatable, dimension(:,:)   :: vec_D
      real(dp),       allocatable, dimension(:,:)   :: vec_Dt
      real(dp),       allocatable, dimension(:,:)   :: DD2_ln_pff_T
      real(dp),       allocatable, dimension(:,:)   :: DD2_ln_pff_T_new
      real(dp),       allocatable, dimension(:,:)   :: work
      integer(int32), allocatable, dimension(:,:)   :: ipiv
   contains
      procedure :: ini    => ini_tgm_fun
      procedure :: cmp    => cmp_tgm_fun
      procedure :: ratio  => ratio_tgm_fun
      procedure :: upd    => upd_tgm_fun
      procedure :: cmp_D  => cmp_DD2_tgm_fun
      procedure :: new_D  => new_DD2_tgm_fun
      procedure :: upd_D  => upd_DD2_tgm_fun
      procedure :: cmp_dl => cmp_dl_tgm_fun
      procedure :: cmp_da => cmp_da_tgm_fun
      procedure :: cpy_dl => cpy_dl_tgm_fun
   end type tgmfun_t
   private :: ini_tgm_fun, cmp_tgm_fun, ratio_tgm_fun, upd_tgm_fun, cmp_DD2_tgm_fun,&
   & upd_DD2_tgm_fun, cmp_dl_tgm_fun, new_DD2_tgm_fun, cmp_da_tgm_fun, cpy_dl_tgm_fun
contains
   subroutine ini_tgm_fun( tgmfun_t_in, n_fe, n_fe_u, n_fe_d, dl_opt, tgm_c )
      class(tgmfun_t), intent(inout) :: tgmfun_t_in
      integer(int32),  intent(in)    :: n_fe, n_fe_u, n_fe_d
      logical,         intent(in)    :: dl_opt
      type(tgmcfs_t),  intent(in)    :: tgm_c
      if( mod(n_fe_u,2_int32).gt.0_int32 ) then
         tgmfun_t_in%d_Tu = n_fe_u + 1_int32
      else
         tgmfun_t_in%d_Tu = n_fe_u
      endif
      if(n_fe_d.gt.0_int32) then
         if (n_fe_d.eq.1_int32) then
            tgmfun_t_in%d_Td = 1_int32
         else
            if( mod(n_fe_d,2_int32).gt.0_int32 ) then
               tgmfun_t_in%d_Td = n_fe_d + 1_int32
            else
               tgmfun_t_in%d_Td = n_fe_d
            endif
         endif
      else
         tgmfun_t_in%d_Td = 0_int32
      endif
      allocate( tgmfun_t_in%pff_T_u(1:n_wlk_max) ) ; tgmfun_t_in%pff_T_u = 1.0_dp
      allocate( tgmfun_t_in%pff_T_d(1:n_wlk_max) ) ; tgmfun_t_in%pff_T_d = 1.0_dp
      if (n_fe_u.gt.1_int32) then
         allocate( tgmfun_t_in%n_mat_chngs_u(1:n_wlk_max) ) ; tgmfun_t_in%n_mat_chngs_u = 0_int32
      endif
      if (n_fe_d.gt.1_int32) then
         allocate( tgmfun_t_in%n_mat_chngs_d(1:n_wlk_max) ) ; tgmfun_t_in%n_mat_chngs_d = 0_int32
      endif
      allocate( tgmfun_t_in%g(1:n_wlk_max) ) ; tgmfun_t_in%g = 1.0_dp
      allocate( tgmfun_t_in%mat_T_u(1:tgmfun_t_in%d_Tu,1:tgmfun_t_in%d_Tu,1:n_wlk_max) )    ; tgmfun_t_in%mat_T_u = 0.0_dp
      allocate( tgmfun_t_in%inv_T_u(1:tgmfun_t_in%d_Tu,1:tgmfun_t_in%d_Tu,1:n_wlk_max) )    ; tgmfun_t_in%inv_T_u = 0.0_dp
      if (tgmfun_t_in%d_Td.gt.0_int32) then
         allocate( tgmfun_t_in%mat_T_d(1:tgmfun_t_in%d_Td,1:tgmfun_t_in%d_Td,1:n_wlk_max) )    ; tgmfun_t_in%mat_T_d = 0.0_dp
         allocate( tgmfun_t_in%inv_T_d(1:tgmfun_t_in%d_Td,1:tgmfun_t_in%d_Td,1:n_wlk_max) )    ; tgmfun_t_in%inv_T_d = 0.0_dp
      endif
      allocate( tgmfun_t_in%mat_Z(1:n_forbs,1:tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,1:n_wlk_max) ) ; tgmfun_t_in%mat_Z = 0.0_dp
      allocate( tgmfun_t_in%mat_D(1:n_forbs,1:tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,1:n_wlk_max) ) ; tgmfun_t_in%mat_D = 0.0_dp
      allocate( tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Tu,1:n_wlk_max) )       ; tgmfun_t_in%vec_T = 0.0_dp
      allocate( tgmfun_t_in%vec_Ti(1:tgmfun_t_in%d_Tu,1:n_wlk_max) )      ; tgmfun_t_in%vec_Ti = 0.0_dp
      allocate( tgmfun_t_in%vec_Tit(1:tgmfun_t_in%d_Tu,1:n_wlk_max) )     ; tgmfun_t_in%vec_Tit = 0.0_dp
      allocate( tgmfun_t_in%vec_D(1:n_forbs,1:n_wlk_max) )     ; tgmfun_t_in%vec_D = 0.0_dp
      allocate( tgmfun_t_in%vec_Dt(1:n_forbs,1:n_wlk_max) )    ; tgmfun_t_in%vec_Dt = 0.0_dp
      allocate( tgmfun_t_in%DD2_ln_pff_T(1:4*n_fe,1:n_wlk_max) ) ; tgmfun_t_in%DD2_ln_pff_T = 0.0_dp
      if ( grad_updat) then
         allocate( tgmfun_t_in%DD2_ln_pff_T_new(1:4,1:n_wlk_max) ); tgmfun_t_in%DD2_ln_pff_T_new = 0.0_dp
      endif
      if ( dl_opt ) then
         if ( tgmfun_t_in%d_Tu.gt.tgmfun_t_in%d_Td) then
            allocate( tgmfun_t_in%mat_OTi(1:n_forbs,1:tgmfun_t_in%d_Tu,1:n_wlk_max) )     ; tgmfun_t_in%mat_OTi  = 0.0_dp
         else
            allocate( tgmfun_t_in%mat_OTi(1:n_forbs,1:tgmfun_t_in%d_Td,1:n_wlk_max) )     ; tgmfun_t_in%mat_OTi  = 0.0_dp
         endif
         allocate( tgmfun_t_in%mat_H(1:n_forbs,1:tgm_c%d_Zu+tgm_c%d_Zd,1:n_wlk_max) ) ; tgmfun_t_in%mat_H  = 0.0_dp
      endif
      if ( aoc_opt ) then
         allocate( tgmfun_t_in%da_ln_pff_T_u(1:n_par_forbs,1:n_wlk_max) ) ; tgmfun_t_in%da_ln_pff_T_u = 0.0_dp
         if( tgmfun_t_in%d_Td.gt.0_int32 ) then
            allocate( tgmfun_t_in%da_ln_pff_T_d(1:n_par_forbs,1:n_wlk_max) ) ; tgmfun_t_in%da_ln_pff_T_d = 0.0_dp
         endif
      endif
      if ( tgmfun_t_in%d_Tu.gt.tgmfun_t_in%d_Td) then
         allocate( tgmfun_t_in%work(1:2*tgmfun_t_in%d_Tu,1:n_wlk_max) )
         allocate( tgmfun_t_in%ipiv(1:tgmfun_t_in%d_Tu,1:n_wlk_max) )
      else
         allocate( tgmfun_t_in%work(1:2*tgmfun_t_in%d_Td,1:n_wlk_max) )
         allocate( tgmfun_t_in%ipiv(1:tgmfun_t_in%d_Td,1:n_wlk_max) )
      endif
   end subroutine ini_tgm_fun
   subroutine cmp_tgm_fun( tgmfun_t_in, iw, n_fe, n_fe_u, n_fe_d, tgm_c, mat_o )
      class(tgmfun_t), intent(inout) :: tgmfun_t_in
      integer(int32),  intent(in) :: iw
      integer(int32),  intent(in) :: n_fe
      integer(int32),  intent(in) :: n_fe_u
      integer(int32),  intent(in) :: n_fe_d
      type(tgmcfs_t),  intent(in) :: tgm_c
      real(dp), dimension(n_forbs,1:n_fe), intent(in) :: mat_o
      call cmp_Z_mat( tgmfun_t_in%d_Tu, tgm_c%d_Zu, n_fe_u, tgm_c%mat_Zu, mat_o(:,1:n_fe_u), &
      & tgmfun_t_in%mat_Z(:,1:tgmfun_t_in%d_Tu,iw) )
      if(n_fe_d.gt.0_int32) then
         if(tgmfun_t_in%d_Td.eq.1_int32) then
            tgmfun_t_in%mat_Z(:,tgmfun_t_in%d_Tu+1,iw) = tgm_c%mat_Zd(:,1)
         else
            call cmp_Z_mat( tgmfun_t_in%d_Td, tgm_c%d_Zd, n_fe_d, tgm_c%mat_Zd, mat_o(:,n_fe_u+1:n_fe), &
            & tgmfun_t_in%mat_Z(:,tgmfun_t_in%d_Tu+1:tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw) )
         endif
      endif
      call cmp_T_mat( tgmfun_t_in%d_Tu, n_fe_u, mat_o(:,1:n_fe_u), tgmfun_t_in%mat_Z(:,1:tgmfun_t_in%d_Tu,iw),&
      & tgmfun_t_in%mat_T_u(:,:,iw) )
      call hh_trn_pff_inv( tgmfun_t_in%d_Tu, tgmfun_t_in%work(:,iw), tgmfun_t_in%ipiv(:,iw), &
      & tgmfun_t_in%mat_T_u(:,:,iw), tgmfun_t_in%pff_T_u(iw), tgmfun_t_in%inv_T_u(:,:,iw) )
      if ( n_fe_d.gt.0_int32 ) then
         if( tgmfun_t_in%d_Td.eq.1_int32 ) then
            tgmfun_t_in%mat_T_d(1,1,iw)= dot_product(tgmfun_t_in%mat_Z(:,tgmfun_t_in%d_Tu+1,iw), mat_o(:,n_fe))
            tgmfun_t_in%pff_T_d(iw) = tgmfun_t_in%mat_T_d(1,1,iw)
            tgmfun_t_in%inv_T_d(1,1,iw) = 1.0_dp / tgmfun_t_in%mat_T_d(1,1,iw)
         else
            call cmp_T_mat ( tgmfun_t_in%d_Td, n_fe_d, mat_o(:,n_fe_u+1:n_fe),&
            & tgmfun_t_in%mat_Z(:,tgmfun_t_in%d_Tu+1:tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw), tgmfun_t_in%mat_T_d(:,:,iw) )
            call hh_trn_pff_inv( tgmfun_t_in%d_Td, tgmfun_t_in%work(1:2*tgmfun_t_in%d_Td,iw), tgmfun_t_in%ipiv(1:tgmfun_t_in%d_Td,iw), &
            & tgmfun_t_in%mat_T_d(:,:,iw), tgmfun_t_in%pff_T_d(iw), tgmfun_t_in%inv_T_d(:,:,iw) )
         endif
      endif
      if (n_fe_u.gt.0_int32) tgmfun_t_in%n_mat_chngs_u(iw) = 0_int32
      if (n_fe_d.gt.0_int32) tgmfun_t_in%n_mat_chngs_d(iw) = 0_int32
   end subroutine cmp_tgm_fun
   subroutine ratio_tgm_fun( tgmfun_t_in, iw, i_fe, n_fe_u, n_fe_d, vec_o_new, g )
      class(tgmfun_t), intent(inout) :: tgmfun_t_in
      integer(int32),               intent(in)    :: iw
      integer(int32),               intent(in)    :: i_fe
      integer(int32),               intent(in)    :: n_fe_u
      integer(int32),               intent(in)    :: n_fe_d
      real(dp), dimension(n_forbs), intent(in)    :: vec_o_new
      real(dp),                     intent(inout) :: g
      if ( i_fe.le.n_fe_u ) then
         call cmp_T_vec ( i_fe, tgmfun_t_in%d_Tu, vec_o_new, tgmfun_t_in%mat_Z(:,1:tgmfun_t_in%d_Tu,iw),&
         & tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Tu,iw) )
         tgmfun_t_in%vec_Ti(1:tgmfun_t_in%d_Tu,iw) = - tgmfun_t_in%inv_T_u(:,i_fe,iw)
         tgmfun_t_in%g(iw) = dot_product( tgmfun_t_in%vec_Ti(1:tgmfun_t_in%d_Tu,iw), tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Tu,iw) )
      else
         if ( n_fe_d.eq.1_int32 ) then
            tgmfun_t_in%vec_T(1,iw) = dot_product( tgmfun_t_in%mat_Z(:,tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw), vec_o_new )
            tgmfun_t_in%g(iw) = tgmfun_t_in%vec_T(1,iw) * tgmfun_t_in%inv_T_d(1,1,iw)
         else
            call cmp_T_vec ( i_fe-n_fe_u, tgmfun_t_in%d_Td, vec_o_new, tgmfun_t_in%mat_Z(:,tgmfun_t_in%d_Tu+1:tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw), &
            & tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Td,iw) )
            tgmfun_t_in%vec_Ti(1:tgmfun_t_in%d_Td,iw) = - tgmfun_t_in%inv_T_d(:,i_fe-n_fe_u,iw)
            tgmfun_t_in%g(iw) = dot_product( tgmfun_t_in%vec_Ti(1:tgmfun_t_in%d_Td,iw), tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Td,iw) )
         endif
      endif
      g = tgmfun_t_in%g(iw)
   end subroutine ratio_tgm_fun
   subroutine upd_tgm_fun( tgmfun_t_in, iw, i_fe, n_fe, n_fe_u, n_fe_d, tgm_c, mat_o )
      class(tgmfun_t), intent(inout) :: tgmfun_t_in
      integer(int32),                    intent(in) :: iw
      integer(int32),                    intent(in) :: i_fe
      integer(int32),                    intent(in) :: n_fe
      integer(int32),                    intent(in) :: n_fe_u
      integer(int32),                    intent(in) :: n_fe_d
      type(tgmcfs_t),                    intent(in) :: tgm_c
      real(dp), dimension(n_forbs,n_fe), intent(in) :: mat_o
      call upd_Z_mat ( i_fe, n_fe, n_fe_u, n_fe_d, tgmfun_t_in%d_Tu, tgmfun_t_in%d_Td, tgm_c, mat_o, tgmfun_t_in%mat_Z(:,:,iw) )
      if(i_fe.le.n_fe_u) then
         tgmfun_t_in%pff_T_u(iw) = tgmfun_t_in%pff_T_u(iw) * tgmfun_t_in%g(iw)
         tgmfun_t_in%mat_T_u(:,i_fe,iw) = tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Tu,iw)
         tgmfun_t_in%mat_T_u(i_fe,:,iw) = - tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Tu,iw)
         if( tgmfun_t_in%n_mat_chngs_u(iw).eq.n_fe_u .or. abs(tgmfun_t_in%g(iw)).lt.10d-14 ) then
            call hh_trn_pff_inv( tgmfun_t_in%d_Tu, tgmfun_t_in%work(:,iw), tgmfun_t_in%ipiv(:,iw), &
            & tgmfun_t_in%mat_T_u(:,:,iw), tgmfun_t_in%pff_T_u(iw), tgmfun_t_in%inv_T_u(:,:,iw) )
            tgmfun_t_in%n_mat_chngs_u(iw) = 0_int32
         else
            call upd_T_inv( i_fe, tgmfun_t_in%d_Tu, tgmfun_t_in%g(iw), tgmfun_t_in%vec_T(:,iw), tgmfun_t_in%vec_Ti(:,iw),&
            & tgmfun_t_in%vec_Tit(:,iw), tgmfun_t_in%inv_T_u(:,:,iw) )
            tgmfun_t_in%n_mat_chngs_u(iw) = tgmfun_t_in%n_mat_chngs_u(iw) + 1_int32
         endif
      else
         tgmfun_t_in%pff_T_d(iw) = tgmfun_t_in%pff_T_d(iw) * tgmfun_t_in%g(iw)
         if(n_fe_d.eq.1_int32) then
            tgmfun_t_in%mat_T_d(1,1,iw) = tgmfun_t_in%vec_T(1,iw)
            tgmfun_t_in%inv_T_d(1,1,iw) = 1.0_dp / tgmfun_t_in%vec_T(1,iw)
         else
            tgmfun_t_in%mat_T_d(:,i_fe-n_fe_u,iw) = tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Td,iw)
            tgmfun_t_in%mat_T_d(i_fe-n_fe_u,:,iw) = - tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Td,iw)
            if( tgmfun_t_in%n_mat_chngs_d(iw).eq.n_fe_d .or. abs(tgmfun_t_in%g(iw)).lt.10d-14) then
               call hh_trn_pff_inv( tgmfun_t_in%d_Td, tgmfun_t_in%work(1:2*tgmfun_t_in%d_Td,iw), tgmfun_t_in%ipiv(1:tgmfun_t_in%d_Td,iw), &
               & tgmfun_t_in%mat_T_d(:,:,iw), tgmfun_t_in%pff_T_d(iw), tgmfun_t_in%inv_T_d(:,:,iw) )
               tgmfun_t_in%n_mat_chngs_d(iw) = 0_int32
            else
               call upd_T_inv( i_fe-n_fe_u, tgmfun_t_in%d_Td, tgmfun_t_in%g(iw), tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Td,iw),&
               & tgmfun_t_in%vec_Ti(1:tgmfun_t_in%d_Td,iw), tgmfun_t_in%vec_Tit(1:tgmfun_t_in%d_Td,iw), tgmfun_t_in%inv_T_d(:,:,iw) )
               tgmfun_t_in%n_mat_chngs_d(iw) = tgmfun_t_in%n_mat_chngs_d(iw) + 1_int32
            endif
         endif
      endif
   end subroutine upd_tgm_fun
   subroutine cmp_DD2_tgm_fun( tgmfun_t_in, iw, n_fe, n_fe_u, n_fe_d, DD2_o )
      class(tgmfun_t), intent(inout) :: tgmfun_t_in
      integer(int32),                      intent(in) :: iw
      integer(int32),                      intent(in) :: n_fe
      integer(int32),                      intent(in) :: n_fe_u
      integer(int32),                      intent(in) :: n_fe_d
      real(dp), dimension(n_forbs,4*n_fe), intent(in) :: DD2_o
      call cmp_D_mat( tgmfun_t_in%d_Tu, tgmfun_t_in%mat_Z(:,1:tgmfun_t_in%d_Tu,iw), &
      & tgmfun_t_in%inv_T_u(:,:,iw), tgmfun_t_in%mat_D(:,1:tgmfun_t_in%d_Tu,iw) )
      if ( n_fe_d .gt. 0_int32 ) then
         if ( n_fe_d .eq. 1_int32 ) then
            tgmfun_t_in%mat_D(:,tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw) = tgmfun_t_in%inv_T_d(1,1,iw) * tgmfun_t_in%mat_Z(:,tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw)
         else
            call cmp_D_mat( tgmfun_t_in%d_Td, tgmfun_t_in%mat_Z(:,tgmfun_t_in%d_Tu+1:tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw), tgmfun_t_in%inv_T_d(:,:,iw), &
            & tgmfun_t_in%mat_D(:,tgmfun_t_in%d_Tu+1:tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw) )
         endif
      endif
      call cmp_DD2_vec( tgmfun_t_in%d_Tu, n_fe, n_fe_u, n_fe_d, tgmfun_t_in%mat_D(:,1:tgmfun_t_in%d_Tu+n_fe_d,iw),&
      & DD2_o, tgmfun_t_in%DD2_ln_pff_T(:,iw) )
   end subroutine cmp_DD2_tgm_fun
   subroutine new_DD2_tgm_fun( tgmfun_t_in, iw, i_fe, n_fe_u, DD2_o_new )
      class(tgmfun_t), intent(inout) :: tgmfun_t_in
      integer(int32),                 intent(in) :: iw
      integer(int32),                 intent(in) :: i_fe
      integer(int32),                 intent(in) :: n_fe_u
      real(dp), dimension(n_forbs,4), intent(in) :: DD2_o_new
      if (i_fe.le.n_fe_u) then
         call dgemv('T', n_forbs, 4, 1.0_dp/tgmfun_t_in%g(iw), DD2_o_new(:,1:4), n_forbs, &
         & tgmfun_t_in%mat_D(:,i_fe,iw), 1_int32, 0.0_dp, tgmfun_t_in%DD2_ln_pff_T_new(:,iw), 1_int32 )
      else
         call dgemv('T', n_forbs, 4, 1.0_dp/tgmfun_t_in%g(iw), DD2_o_new(:,1:4), n_forbs, &
         & tgmfun_t_in%mat_D(:,tgmfun_t_in%d_Tu+i_fe-n_fe_u,iw), 1_int32, 0.0_dp, tgmfun_t_in%DD2_ln_pff_T_new(:,iw), 1_int32 )
      endif
      tgmfun_t_in%DD2_ln_pff_T_new(4,iw) = tgmfun_t_in%DD2_ln_pff_T_new(4,iw) - sum(tgmfun_t_in%DD2_ln_pff_T_new(1:3,iw)**2)
   end subroutine new_DD2_tgm_fun
   subroutine upd_DD2_tgm_fun( tgmfun_t_in, iw, i_fe, n_fe, n_fe_u, n_fe_d, DD2_o )
      class(tgmfun_t), intent(inout) :: tgmfun_t_in
      integer(int32),                      intent(in) :: iw
      integer(int32),                      intent(in) :: i_fe
      integer(int32),                      intent(in) :: n_fe
      integer(int32),                      intent(in) :: n_fe_u
      integer(int32),                      intent(in) :: n_fe_d
      real(dp), dimension(n_forbs,4*n_fe), intent(in) :: DD2_o
      if( i_fe.le.n_fe_u ) then
         if ( tgmfun_t_in%n_mat_chngs_u(iw).eq.0_int32 ) then
            call cmp_D_mat( tgmfun_t_in%d_Tu, tgmfun_t_in%mat_Z(:,1:tgmfun_t_in%d_Tu,iw), tgmfun_t_in%inv_T_u(:,:,iw), tgmfun_t_in%mat_D(:,1:tgmfun_t_in%d_Tu,iw) )
         else
            call upd_D_mat( i_fe, tgmfun_t_in%d_Tu, tgmfun_t_in%mat_Z(:,i_fe,iw), tgmfun_t_in%vec_T(:,iw),&
            & tgmfun_t_in%vec_Tit(:,iw), tgmfun_t_in%vec_Ti(:,iw), tgmfun_t_in%g(iw), &
            & tgmfun_t_in%vec_D(:,iw), tgmfun_t_in%vec_Dt(:,iw), tgmfun_t_in%mat_D(:,1:tgmfun_t_in%d_Tu,iw) )
         endif
         call upd_DD2_vec( n_fe_u, tgmfun_t_in%mat_D(:,1:n_fe_u,iw), DD2_o(:,1:4*n_fe_u), tgmfun_t_in%DD2_ln_pff_T(1:4*n_fe_u,iw) )
      else
         if (tgmfun_t_in%d_Td.eq.1_int32) then
            tgmfun_t_in%mat_D(:,tgmfun_t_in%d_Tu+n_fe_d,iw) = tgmfun_t_in%inv_T_d(1,1,iw) * tgmfun_t_in%mat_Z(:,tgmfun_t_in%d_Tu+n_fe_d,iw)
         else
            if ( tgmfun_t_in%n_mat_chngs_d(iw).eq.0_int32 ) then
               call cmp_D_mat( tgmfun_t_in%d_Td, tgmfun_t_in%mat_Z(:,tgmfun_t_in%d_Tu+1:tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw), tgmfun_t_in%inv_T_d(:,:,iw), &
               & tgmfun_t_in%mat_D(:,tgmfun_t_in%d_Tu+1:tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw) )
            else
               call upd_D_mat( i_fe-n_fe_u, tgmfun_t_in%d_Td, tgmfun_t_in%mat_Z(:,i_fe-n_fe_u+tgmfun_t_in%d_Tu,iw), &
               & tgmfun_t_in%vec_T(1:tgmfun_t_in%d_Td,iw), tgmfun_t_in%vec_Tit(1:tgmfun_t_in%d_Td,iw), tgmfun_t_in%vec_Ti(1:tgmfun_t_in%d_Td,iw), &
               & tgmfun_t_in%g(iw), tgmfun_t_in%vec_D(:,iw), tgmfun_t_in%vec_Dt(:,iw), tgmfun_t_in%mat_D(:,tgmfun_t_in%d_Tu+1:tgmfun_t_in%d_Tu+tgmfun_t_in%d_Td,iw) )
            endif
         endif
         call upd_DD2_vec( n_fe_d, tgmfun_t_in%mat_D(:,tgmfun_t_in%d_Tu+1:tgmfun_t_in%d_Tu+n_fe_d,iw),&
         & DD2_o(:,4*n_fe_u+1:4*n_fe), tgmfun_t_in%DD2_ln_pff_T(4*n_fe_u+1:4*n_fe,iw) )
      endif
   end subroutine upd_DD2_tgm_fun
   subroutine cmp_dl_tgm_fun( tgmfun_t_in, iw, n_fe, n_fe_u, n_fe_d, tgm_c, mat_o, n_par_det, dp_ln_pff_T )
      class(tgmfun_t), intent(inout) :: tgmfun_t_in
      integer(int32),                    intent(in)    :: iw
      integer(int32),                    intent(in)    :: n_fe
      integer(int32),                    intent(in)    :: n_fe_u
      integer(int32),                    intent(in)    :: n_fe_d
      type(tgmcfs_t),                    intent(in)    :: tgm_c
      real(dp), dimension(n_forbs,n_fe), intent(in)    :: mat_o
      integer(int32),                    intent(in)    :: n_par_det
      real(dp), dimension(n_par_det),    intent(inout) :: dp_ln_pff_T
      call cmp_H_mat( n_fe_u, tgmfun_t_in%d_Tu, tgm_c%d_Zu, mat_o(:,1:n_fe_u), tgmfun_t_in%inv_T_u(:,:,iw),&
      & tgmfun_t_in%mat_OTi(:,1:tgmfun_t_in%d_Tu,iw), tgmfun_t_in%mat_H(:,1:tgm_c%d_Zu,iw))
      if ( n_fe_d.gt.0_int32 ) then
         if ( tgmfun_t_in%d_Td.eq.1_int32 ) then
            tgmfun_t_in%mat_H(:,tgm_c%d_Zu+1,iw) = tgmfun_t_in%inv_T_d(1,1,iw) * mat_o(:,n_fe)
         else
            call cmp_H_mat( n_fe_d, tgmfun_t_in%d_Td, tgm_c%d_Zd, mat_o(:,n_fe_u+1:n_fe), tgmfun_t_in%inv_T_d(:,:,iw), &
            & tgmfun_t_in%mat_OTi(:,1:tgmfun_t_in%d_Td,iw), tgmfun_t_in%mat_H(:,tgm_c%d_Zu+1:tgm_c%d_Zu+tgm_c%d_Zd,iw) )
         endif
      endif
      call tgmfun_t_in%cpy_dl( iw, n_fe_d, tgm_c, n_par_det, dp_ln_pff_T )
   end subroutine cmp_dl_tgm_fun
   subroutine cmp_da_tgm_fun( tgmfun_t_in, iw, n_fe, n_fe_u, n_fe_d, mat_da_o, dp_ln_pff_T )
      class(tgmfun_t), intent(inout) :: tgmfun_t_in
      integer(int32),                        intent(in)    :: iw
      integer(int32),                        intent(in)    :: n_fe
      integer(int32),                        intent(in)    :: n_fe_u
      integer(int32),                        intent(in)    :: n_fe_d
      real(dp), dimension(n_par_forbs,n_fe), intent(in)    :: mat_da_o
      real(dp), dimension(n_par_forbs),      intent(inout) :: dp_ln_pff_T
      call cmp_da_vec( n_fe_u, mat_da_o(:,1:n_fe_u), &
      & tgmfun_t_in%mat_D(:,1:n_fe_u,iw), tgmfun_t_in%da_ln_pff_T_u(:,iw) )
      dp_ln_pff_T = dp_ln_pff_T + tgmfun_t_in%da_ln_pff_T_u(:,iw)
      if ( n_fe_d.gt.0_int32 ) then
         call cmp_da_vec( n_fe_d, mat_da_o(:,n_fe_u+1:n_fe), &
         & tgmfun_t_in%mat_D(:,tgmfun_t_in%d_Tu+1:tgmfun_t_in%d_Tu+n_fe_d,iw), tgmfun_t_in%da_ln_pff_T_d(:,iw) )
         dp_ln_pff_T = dp_ln_pff_T + tgmfun_t_in%da_ln_pff_T_d(:,iw)
      endif
   end subroutine cmp_da_tgm_fun
   subroutine cpy_dl_tgm_fun( tgmfun_t_in, iw, n_fe_d, tgm_c, n_par_det, dp_ln_pff_T )
      class(tgmfun_t), intent(inout) :: tgmfun_t_in
      integer(int32),                 intent(in) :: iw
      integer(int32),                 intent(in)    :: n_fe_d
      type(tgmcfs_t),  target,        intent(in)    :: tgm_c
      integer(int32),                 intent(in)    :: n_par_det
      real(dp), dimension(n_par_det), intent(inout) :: dp_ln_pff_T
      integer(int32), pointer :: o1, o2
      integer(int32) :: ip
      integer(int32) :: i1
      ip = 1_int32
      if (spin.eq.'R') then
         do i1 = 1_int32, tgm_c%n_nz_Z
            o1 => tgm_c%nz_Z(i1)%o1
            o2 => tgm_c%nz_Z(i1)%o2
            dp_ln_pff_T(ip) = tgmfun_t_in%mat_H(o1,o2,iw)
            if ( n_fe_d.gt.1_int32 ) then
               dp_ln_pff_T(ip) = dp_ln_pff_T(ip) + tgmfun_t_in%mat_H(o1,o2+tgm_c%d_Zu,iw)
            endif
            ip = ip + 1_int32
         enddo
         if ( tgm_c%n_nz_Zu_u.gt.0_int32 ) then
            do i1 = 1_int32, tgm_c%n_nz_Zu_u
               o1 => tgm_c%nz_Zu_u(i1)%o1
               o2 => tgm_c%nz_Zu_u(i1)%o2
               dp_ln_pff_T(ip) = tgmfun_t_in%mat_H(o1,o2,iw)
               if ( tgm_c%n_nz_Zd_u.gt.0_int32  ) then
                  dp_ln_pff_T(ip) = dp_ln_pff_T(ip) + tgmfun_t_in%mat_H(o1,tgm_c%d_Zd+tgm_c%d_Zu,iw)
               endif
               ip = ip + 1_int32
            enddo
         else
            if ( tgm_c%n_nz_Zd_u.gt.0_int32  ) then
               do i1 = 1_int32, tgm_c%n_nz_Zd_u
                  o1 => tgm_c%nz_Zd_u(i1)%o1
                  dp_ln_pff_T(ip) = tgmfun_t_in%mat_H(o1,tgm_c%d_Zu+tgm_c%d_Zd,iw)
                  ip = ip + 1_int32
               enddo
            endif
         endif
      else ! if ( spin.eq.'U')
         do i1 = 1_int32, tgm_c%n_nz_Z
            o1 => tgm_c%nz_Z(i1)%o1
            o2 => tgm_c%nz_Z(i1)%o2
            dp_ln_pff_T(ip) = tgmfun_t_in%mat_H(o1,o2,iw)
            ip = ip + 1_int32
         enddo
         if ( tgm_c%n_nz_Zu_u.gt.0_int32 ) then
            do i1 = 1_int32, tgm_c%n_nz_Zu_u
               o1 => tgm_c%nz_Zu_u(i1)%o1
               o2 => tgm_c%nz_Zu_u(i1)%o2
               dp_ln_pff_T(ip) = tgmfun_t_in%mat_H(o1,o2,iw)
               ip = ip + 1_int32
            enddo
         endif
         if ( tgm_c%d_Zd.ge.n_forbs ) then
            do i1 = 1_int32, tgm_c%n_nz_Z
               o1 => tgm_c%nz_Z(i1)%o1
               o2 => tgm_c%nz_Z(i1)%o2
               dp_ln_pff_T(ip) = tgmfun_t_in%mat_H(o1,o2+tgm_c%d_Zu,iw)
               ip = ip + 1_int32
            enddo
         endif
         if ( tgm_c%n_nz_Zd_u.gt.0_int32) then
            do i1 = 1_int32, tgm_c%n_nz_Zd_u
               o1 => tgm_c%nz_Zd_u(i1)%o1
               dp_ln_pff_T(ip) = tgmfun_t_in%mat_H(o1,tgm_c%d_Zd+tgm_c%d_Zu,iw)
               ip = ip + 1_int32
            enddo
         endif
      endif
   end subroutine cpy_dl_tgm_fun
end module triplet_geminal_c
