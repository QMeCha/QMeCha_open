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
module slater_determinant_c
   use fortran_kinds_v,              only: dp, int32
   use orbmat_cls,                   only: matorb_t
   use fermionic_wavefunction_v,     only: spin
   use slater_determinant_params_m,  only: sldcfs_t
   use fermionic_orbitals_m,         only: n_forbs, n_par_forbs, aoc_opt
   use quantum_monte_carlo_v,        only: grad_updat, n_wlk_max
   use slater_determinant_operats_m, only: cmp_S_mat, cmp_S_vec, upd_S_inv_det,&
   & cmp_D_mat, upd_D_mat, cmp_K_mat, cmp_da_sld
   implicit none
   type, public :: sldfun_t
      real(dp), allocatable, dimension(:) :: det_S_u, det_S_d, g
      integer(int32), allocatable, dimension(:) :: n_mat_chngs_u,  n_mat_chngs_d
      real(dp), allocatable, dimension(:,:,:) :: mat_S_u, inv_S_u
      real(dp), allocatable, dimension(:,:,:) :: mat_S_d, inv_S_d
      real(dp), allocatable, dimension(:,:,:) :: mat_D
      real(dp), allocatable, dimension(:,:,:) :: mat_K_u, mat_K_d
      real(dp), allocatable, dimension(:,:)   :: da_ln_det_S_d, da_ln_det_S_u
      real(dp), allocatable, dimension(:,:)   :: vec_S_new
      real(dp), allocatable, dimension(:,:)   :: vec_v
      real(dp), allocatable, dimension(:,:)   :: vec_v_s
      real(dp), allocatable, dimension(:,:)   :: vec_v_d
      real(dp), allocatable, dimension(:,:)   :: DD2_ln_det_S
      real(dp), allocatable, dimension(:,:)   :: DD2_ln_det_S_new
   contains
      procedure :: ini   => ini_sld_fun
      procedure :: cmp   => cmp_sld_fun
      procedure :: ratio => ratio_sld_fun
      procedure :: upd   => upd_sld_fun
      procedure :: cmp_D => cmp_DD2_sld_fun
      procedure :: upd_D => upd_DD2_sld_fun
      procedure :: new_D => new_DD2_sld_fun
      procedure :: cmp_dl => cmp_dl_sld_fun
      procedure :: cmp_da => cmp_da_sld_fun
   end type sldfun_t
   private :: ini_sld_fun, cmp_sld_fun, ratio_sld_fun, upd_sld_fun, cmp_DD2_sld_fun, &
   & upd_DD2_sld_fun, new_DD2_sld_fun, cmp_dl_sld_fun, cmp_da_sld_fun
contains
   subroutine ini_sld_fun( sldfun_t_in, n_fe, n_fe_u, n_fe_d, dl_opt  )
      class(sldfun_t), intent(inout) :: sldfun_t_in
      integer(int32),  intent(in)    :: n_fe, n_fe_u, n_fe_d
      logical,         intent(in)    :: dl_opt
      allocate( sldfun_t_in%g(1:n_wlk_max) ) ; sldfun_t_in%g(1:n_wlk_max) = 1.0_dp
      allocate( sldfun_t_in%det_S_u(1:n_wlk_max) ) ; sldfun_t_in%det_S_u(1:n_wlk_max) = 1.0_dp
      if (n_fe_u.gt.0_int32) then
         allocate( sldfun_t_in%n_mat_chngs_u(1:n_wlk_max) ) ; sldfun_t_in%n_mat_chngs_u(1:n_wlk_max) = 0_int32
         allocate( sldfun_t_in%mat_S_u(1:n_fe_u,1:n_fe_u,1:n_wlk_max) ) ; sldfun_t_in%mat_S_u = 0.0_dp
         allocate( sldfun_t_in%inv_S_u(1:n_fe_u,1:n_fe_u,1:n_wlk_max) ) ; sldfun_t_in%inv_S_u = 0.0_dp
      endif
      allocate( sldfun_t_in%det_S_d(1:n_wlk_max) ) ; sldfun_t_in%det_S_d(1:n_wlk_max) = 1.0_dp
      if (n_fe_d.gt.0_int32) then
         allocate( sldfun_t_in%n_mat_chngs_d(1:n_wlk_max) ) ; sldfun_t_in%n_mat_chngs_d(1:n_wlk_max) = 0_int32
         allocate( sldfun_t_in%mat_S_d(1:n_fe_d,1:n_fe_d,1:n_wlk_max) ) ; sldfun_t_in%mat_S_d = 0.0_dp
         allocate( sldfun_t_in%inv_S_d(1:n_fe_d,1:n_fe_d,1:n_wlk_max) ) ; sldfun_t_in%inv_S_d = 0.0_dp
      endif
      allocate( sldfun_t_in%mat_D(1:n_forbs,1:n_fe,1:n_wlk_max) )    ; sldfun_t_in%mat_D = 0.0_dp
      if (n_fe_u.ge.n_fe_d) then
         allocate( sldfun_t_in%vec_S_new(1:n_fe_u,1:n_wlk_max) )    ; sldfun_t_in%vec_S_new = 0.0_dp
         allocate( sldfun_t_in%vec_v(1:n_fe_u,1:n_wlk_max) )        ; sldfun_t_in%vec_v = 0.0_dp
         allocate( sldfun_t_in%vec_v_s(1:n_fe_u,1:n_wlk_max) )      ; sldfun_t_in%vec_v_s = 0.0_dp
      else
         allocate( sldfun_t_in%vec_S_new(1:n_fe_d,1:n_wlk_max) )    ; sldfun_t_in%vec_S_new = 0.0_dp
         allocate( sldfun_t_in%vec_v(1:n_fe_d,1:n_wlk_max) )        ; sldfun_t_in%vec_v = 0.0_dp
         allocate( sldfun_t_in%vec_v_s(1:n_fe_d,1:n_wlk_max) )      ; sldfun_t_in%vec_v_s = 0.0_dp
      endif
      allocate( sldfun_t_in%vec_v_d(1:n_forbs,1:n_wlk_max) )     ; sldfun_t_in%vec_v_d = 0.0_dp
      allocate( sldfun_t_in%DD2_ln_det_S(1:4*n_fe,1:n_wlk_max) ) ; sldfun_t_in%DD2_ln_det_S = 0.0_dp
      if ( grad_updat ) then
         allocate( sldfun_t_in%DD2_ln_det_S_new(1:4,1:n_wlk_max) ) ; sldfun_t_in%DD2_ln_det_S_new = 0.0_dp
      endif
      if ( dl_opt ) then
         if( n_fe_u.gt.0_int32 ) then
            allocate( sldfun_t_in%mat_K_u(1:n_forbs,1:n_fe_u,1:n_wlk_max) ) ; sldfun_t_in%mat_K_u = 0.0_dp
         endif
         if( n_fe_d.gt.0_int32 ) then
            allocate( sldfun_t_in%mat_K_d(1:n_forbs,1:n_fe_d,1:n_wlk_max) ) ; sldfun_t_in%mat_K_d = 0.0_dp
         endif
      endif
      if ( aoc_opt ) then
         if( n_fe_u.gt.0_int32 ) then
            allocate( sldfun_t_in%da_ln_det_S_u(1:n_par_forbs,1:n_wlk_max) ) ; sldfun_t_in%da_ln_det_S_u = 0.0_dp
         endif
         if( n_fe_d.gt.0_int32 ) then
            allocate( sldfun_t_in%da_ln_det_S_d(1:n_par_forbs,1:n_wlk_max) ) ; sldfun_t_in%da_ln_det_S_d = 0.0_dp
         endif
      endif
   end subroutine ini_sld_fun
   subroutine cmp_sld_fun( sldfun_t_in, iw, n_fe, n_fe_u, n_fe_d, sld_c, mat_o )
      class(sldfun_t), intent(inout) :: sldfun_t_in
      integer(int32),                     intent(in) :: iw
      integer(int32),                     intent(in) :: n_fe, n_fe_u, n_fe_d
      type(sldcfs_t),                     intent(in) :: sld_c
      real(dp),  dimension(n_forbs,n_fe), intent(in) :: mat_o
      if (n_fe_u.gt.0_int32) then
         call cmp_S_mat( n_forbs, n_fe_u, mat_o(:,1:n_fe_u), sld_c%mat_L_u, &
         & sldfun_t_in%mat_S_u(:,:,iw) )
         call lu_dec_det_inv( n_fe_u, sldfun_t_in%mat_S_u(:,:,iw), sldfun_t_in%det_S_u(iw), sldfun_t_in%inv_S_u(:,:,iw) )
         sldfun_t_in%n_mat_chngs_u(iw) = 0_int32
      endif
      if (n_fe_d.gt.0_int32) then
         call cmp_S_mat( n_forbs, n_fe_d, mat_o(:,n_fe_u+1:n_fe), sld_c%mat_L_d, &
         & sldfun_t_in%mat_S_d(:,:,iw) )
         call lu_dec_det_inv( n_fe_d, sldfun_t_in%mat_S_d(:,:,iw), sldfun_t_in%det_S_d(iw), sldfun_t_in%inv_S_d(:,:,iw) )
         sldfun_t_in%n_mat_chngs_d(iw) = 0_int32
      endif
   end subroutine cmp_sld_fun
   subroutine ratio_sld_fun( sldfun_t_in, iw, i_fe, n_fe_u, n_fe_d, sld_c, vec_o_new, g )
      class(sldfun_t), intent(inout) :: sldfun_t_in
      integer(int32),                intent(in) :: iw
      integer(int32),                intent(in)  :: i_fe
      integer(int32),                intent(in)  :: n_fe_u, n_fe_d
      type(sldcfs_t),                intent(in)  :: sld_c
      real(dp),  dimension(n_forbs), intent(in)  :: vec_o_new
      real(dp),                      intent(out) :: g
      sldfun_t_in%vec_S_new(:,iw) = 0.0_dp
      if( i_fe.le.n_fe_u ) then
         call cmp_S_vec( n_forbs, n_fe_u, vec_o_new, sld_c%mat_L_u, sldfun_t_in%vec_S_new(1:n_fe_u,iw) )
         sldfun_t_in%vec_v_s(1:n_fe_u,iw) = sldfun_t_in%inv_S_u(i_fe,1:n_fe_u,iw)
         sldfun_t_in%g(iw) = dot_product( sldfun_t_in%vec_v_s(1:n_fe_u,iw), sldfun_t_in%vec_S_new(1:n_fe_u,iw) )
      else
         call cmp_S_vec( n_forbs, n_fe_d, vec_o_new, sld_c%mat_L_d, sldfun_t_in%vec_S_new(1:n_fe_d,iw) )
         sldfun_t_in%vec_v_s(1:n_fe_d,iw) = sldfun_t_in%inv_S_d(i_fe-n_fe_u,1:n_fe_d,iw)
         sldfun_t_in%g(iw) = dot_product( sldfun_t_in%vec_v_s(1:n_fe_d,iw), sldfun_t_in%vec_S_new(1:n_fe_d,iw) )
      endif
      g = sldfun_t_in%g(iw)
   end subroutine ratio_sld_fun
   subroutine new_DD2_sld_fun( sldfun_t_in, iw, i_fe, vec_DD2_o_new )
      class(sldfun_t), intent(inout) :: sldfun_t_in
      integer(int32),                 intent(in) :: iw
      integer(int32),                 intent(in) :: i_fe
      real(dp), dimension(n_forbs,4), intent(in) :: vec_DD2_o_new
      call dgemv('T', n_forbs, 4, 1.0_dp / sldfun_t_in%g(iw), vec_DD2_o_new, n_forbs, &
      & sldfun_t_in%mat_D(:,i_fe,iw), 1_int32, 0.0_dp, sldfun_t_in%DD2_ln_det_S_new(:,iw), 1_int32 )
      sldfun_t_in%DD2_ln_det_S_new(4,iw) = sldfun_t_in%DD2_ln_det_S_new(4,iw) - sum(sldfun_t_in%DD2_ln_det_S_new(1:3,iw)**2)
   end subroutine new_DD2_sld_fun
   subroutine upd_sld_fun( sldfun_t_in, iw, i_fe, n_fe_u, n_fe_d )
      class(sldfun_t), intent(inout) :: sldfun_t_in
      integer(int32),  intent(in)    :: iw
      integer(int32),  intent(in)    :: i_fe
      integer(int32),  intent(in)    :: n_fe_u, n_fe_d
      if( i_fe.le.n_fe_u ) then
         call upd_S_inv_det( i_fe, n_fe_u, sldfun_t_in%n_mat_chngs_u(iw), sldfun_t_in%mat_S_u(:,:,iw), sldfun_t_in%inv_S_u(:,:,iw), sldfun_t_in%det_S_u(iw), &
         & sldfun_t_in%vec_S_new(1:n_fe_u,iw), sldfun_t_in%vec_v_s(1:n_fe_u,iw), sldfun_t_in%vec_v(1:n_fe_u,iw), sldfun_t_in%g(iw) )
      else
         call upd_S_inv_det( i_fe-n_fe_u, n_fe_d, sldfun_t_in%n_mat_chngs_d(iw), sldfun_t_in%mat_S_d(:,:,iw), sldfun_t_in%inv_S_d(:,:,iw), sldfun_t_in%det_S_d(iw), &
         & sldfun_t_in%vec_S_new(1:n_fe_d,iw), sldfun_t_in%vec_v_s(1:n_fe_d,iw), sldfun_t_in%vec_v(1:n_fe_d,iw), sldfun_t_in%g(iw) )
      endif
   end subroutine upd_sld_fun
   subroutine cmp_DD2_sld_fun( sldfun_t_in, iw, n_fe, n_fe_u, n_fe_d, sld_c, DD2_o )
      class(sldfun_t), intent(inout) :: sldfun_t_in
      integer(int32),                       intent(in) :: iw
      integer(int32),                       intent(in) :: n_fe, n_fe_u, n_fe_d
      type(sldcfs_t),                       intent(in) :: sld_c
      real(dp),  dimension(n_forbs,4*n_fe), intent(in) :: DD2_o
      integer(int32) :: i1
      if(n_fe_u.gt.0_int32) &
      & call cmp_D_mat( n_fe_u, sld_c%mat_L_u, sldfun_t_in%inv_S_u(:,:,iw), sldfun_t_in%mat_D(:,1:n_fe_u,iw) )
      if(n_fe_d.gt.0_int32) &
      & call cmp_D_mat( n_fe_d, sld_c%mat_L_d, sldfun_t_in%inv_S_d(:,:,iw), sldfun_t_in%mat_D(:,n_fe_u+1:n_fe,iw) )
      do i1 = 1_int32, n_fe
         call dgemv('T', n_forbs, 4, 1.0_dp, DD2_o(:,4*i1-3:4*i1), n_forbs, &
         & sldfun_t_in%mat_D(:,i1,iw), 1_int32, 0.0_dp, sldfun_t_in%DD2_ln_det_S(4*i1-3:4*i1,iw), 1_int32)
         sldfun_t_in%DD2_ln_det_S(4*i1,iw) = sldfun_t_in%DD2_ln_det_S(4*i1,iw) - sum(sldfun_t_in%DD2_ln_det_S(4*i1-3:4*i1-1,iw)**2)
      enddo
   end subroutine cmp_DD2_sld_fun
   subroutine upd_DD2_sld_fun( sldfun_t_in, iw, i_fe, n_fe, n_fe_u, n_fe_d, sld_c, DD2_o )
      class(sldfun_t), intent(inout) :: sldfun_t_in
      integer(int32),                      intent(in) :: iw
      integer(int32),                      intent(in) :: i_fe
      integer(int32),                      intent(in) :: n_fe
      integer(int32),                      intent(in) :: n_fe_u
      integer(int32),                      intent(in) :: n_fe_d
      type(sldcfs_t),                      intent(in) :: sld_c
      real(dp), dimension(n_forbs,4*n_fe), intent(in) :: DD2_o
      integer(int32) :: i1
      if( i_fe.le.n_fe_u ) then
         call upd_D_mat( n_fe_u, i_fe, sldfun_t_in%n_mat_chngs_u(iw), sld_c%mat_L_u, sldfun_t_in%inv_S_u(:,:,iw),&
         & sldfun_t_in%mat_D(:,1:n_fe_u,iw ), sldfun_t_in%vec_v_d(:,iw), sldfun_t_in%vec_v(1:n_fe_u,iw), sldfun_t_in%g(iw) )
         do i1 = 1_int32, n_fe_u
            call dgemv('T', n_forbs, 4, 1.0_dp, DD2_o(:,4*i1-3:4*i1), n_forbs, &
            & sldfun_t_in%mat_D(:,i1,iw), 1_int32, 0.0_dp, sldfun_t_in%DD2_ln_det_S(4*i1-3:4*i1,iw), 1_int32)
            sldfun_t_in%DD2_ln_det_S(4*i1,iw) = sldfun_t_in%DD2_ln_det_S(4*i1,iw) - sum(sldfun_t_in%DD2_ln_det_S(4*i1-3:4*i1-1,iw)**2)
         enddo
      else
         call upd_D_mat( n_fe_d, i_fe-n_fe_u, sldfun_t_in%n_mat_chngs_d(iw), sld_c%mat_L_d, sldfun_t_in%inv_S_d(:,:,iw),&
         & sldfun_t_in%mat_D(:,n_fe_u+1:n_fe,iw), sldfun_t_in%vec_v_d(:,iw), sldfun_t_in%vec_v(1:n_fe_d,iw), sldfun_t_in%g(iw))
         do i1 = n_fe_u+1, n_fe
            call dgemv('T', n_forbs, 4, 1.0_dp, DD2_o(:,4*i1-3:4*i1), n_forbs, &
            & sldfun_t_in%mat_D(:,i1,iw), 1_int32, 0.0_dp, sldfun_t_in%DD2_ln_det_S(4*i1-3:4*i1,iw), 1_int32)
            sldfun_t_in%DD2_ln_det_S(4*i1,iw) = sldfun_t_in%DD2_ln_det_S(4*i1,iw) - sum(sldfun_t_in%DD2_ln_det_S(4*i1-3:4*i1-1,iw)**2)
         enddo
      endif
   end subroutine upd_DD2_sld_fun
   subroutine cmp_dl_sld_fun( sldfun_t_in, iw, n_fe, n_fe_u, n_fe_d, sld_c, mat_o, &
   & n_par_det, dp_ln_det_S )
      class(sldfun_t), intent(inout) :: sldfun_t_in
      integer(int32),                    intent(in) :: iw
      integer(int32),                    intent(in)    :: n_fe
      integer(int32),                    intent(in)    :: n_fe_u
      integer(int32),                    intent(in)    :: n_fe_d
      type(sldcfs_t),                    intent(in)    :: sld_c
      real(dp), dimension(n_forbs,n_fe), intent(in)    :: mat_o
      integer(int32),                    intent(in)    :: n_par_det
      real(dp), dimension(n_par_det),    intent(inout) :: dp_ln_det_S
      integer(int32) :: i1, ip
      ip = 1_int32
      call cmp_K_mat( n_fe_u, mat_o(:,1:n_fe_u), sldfun_t_in%inv_S_u(:,:,iw), sldfun_t_in%mat_K_u(:,:,iw) )
      do i1 = 1_int32, sld_c%n_nz_Lu
         dp_ln_det_S(ip) = sldfun_t_in%mat_K_u(sld_c%nz_Lu(i1)%o1, sld_c%nz_Lu(i1)%o2,iw)
         ip = ip + 1_int32
      enddo
      if( n_fe_d.gt.0_int32 ) then
         call cmp_K_mat( n_fe_d, mat_o(:,n_fe_u+1:n_fe), sldfun_t_in%inv_S_d(:,:,iw), sldfun_t_in%mat_K_d(:,:,iw) )
         if (spin.eq.'R') then
            ip = 1_int32
            do i1 = 1_int32, sld_c%n_nz_Ld
               dp_ln_det_S(ip) = dp_ln_det_S(ip) + sldfun_t_in%mat_K_d(sld_c%nz_Ld(i1)%o1, sld_c%nz_Ld(i1)%o2,iw)
               ip = ip + 1_int32
            enddo
            ip = sld_c%n_nz_Lu + 1_int32
         else
            do i1 = 1_int32, sld_c%n_nz_Ld
               dp_ln_det_S(ip) = sldfun_t_in%mat_K_d(sld_c%nz_Ld(i1)%o1, sld_c%nz_Ld(i1)%o2,iw)
               ip = ip + 1_int32
            enddo
         endif
      endif
   end subroutine cmp_dl_sld_fun
   subroutine cmp_da_sld_fun( sldfun_t_in, iw, n_fe, n_fe_u, n_fe_d, mat_da_o, da_ln_det_S )
      class(sldfun_t), intent(inout) :: sldfun_t_in
      integer(int32),                        intent(in) :: iw
      integer(int32),                        intent(in)    :: n_fe
      integer(int32),                        intent(in)    :: n_fe_u
      integer(int32),                        intent(in)    :: n_fe_d
      real(dp), dimension(n_par_forbs,n_fe), intent(in)    :: mat_da_o
      real(dp), dimension(n_par_forbs),      intent(inout) :: da_ln_det_S
      call cmp_da_sld( n_fe_u, mat_da_o(:,1:n_fe_u), sldfun_t_in%mat_D(:,1:n_fe_u,iw), sldfun_t_in%da_ln_det_S_u(:,iw) )
      da_ln_det_S = da_ln_det_S + sldfun_t_in%da_ln_det_S_u(:,iw)
      if( n_fe_d.gt.0_int32 ) then
         call cmp_da_sld( n_fe_d, mat_da_o(:,n_fe_u+1:n_fe), sldfun_t_in%mat_D(:,n_fe_u+1:n_fe,iw), sldfun_t_in%da_ln_det_S_d(:,iw) )
         da_ln_det_S  = da_ln_det_S + sldfun_t_in%da_ln_det_S_d(:,iw)
      endif
   end subroutine cmp_da_sld_fun
end module slater_determinant_c
