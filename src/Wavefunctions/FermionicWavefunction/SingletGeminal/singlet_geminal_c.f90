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
module singlet_geminal_c
   use fortran_kinds_v,           only: dp, int32
   use molecular_system_v,        only: n_at
   use singlet_geminal_params_m,  only: sgmcfs_t
   use fermionic_orbitals_m,      only: n_forbs, n_par_forbs, aoc_opt, forbs
   use quantum_monte_carlo_v,     only: grad_updat, n_wlk_max
   use singlet_geminal_operats_m, only: cmp_L_mat, cmp_G_mat, cmp_G_vec, cmp_D_mat,&
   & upd_D_mat, cmp_DD2_vec, upd_L_mat, upd_G_inv, cmp_KH_mat
   implicit none
   type, public :: sgmfun_t
      real(dp),       allocatable, dimension(:) :: det_G, g
      integer(int32), allocatable, dimension(:) :: n_mat_chngs
      real(dp), allocatable, dimension(:,:,:) :: mat_G
      real(dp), allocatable, dimension(:,:,:) :: inv_G
      real(dp), allocatable, dimension(:,:,:) :: mat_L
      real(dp), allocatable, dimension(:,:,:) :: mat_D
      real(dp), allocatable, dimension(:,:,:) :: mat_K
      real(dp), allocatable, dimension(:,:,:) :: mat_H
      real(dp), allocatable, dimension(:,:)   :: vec_G
      real(dp), allocatable, dimension(:,:)   :: vec_Gig
      real(dp), allocatable, dimension(:,:)   :: vec_Gi
      real(dp), allocatable, dimension(:,:)   :: vec_D
      real(dp), allocatable, dimension(:,:)   :: vec_Dg
      real(dp), allocatable, dimension(:,:)   :: da_ln_det_G
      real(dp), allocatable, dimension(:,:)   :: DD2_ln_det_G
      real(dp), allocatable, dimension(:,:)   :: DD2_ln_det_G_new
   contains
      procedure :: ini    => ini_sgm_fun
      procedure :: cmp    => cmp_sgm_fun
      procedure :: ratio  => ratio_sgm_fun
      procedure :: upd    => upd_sgm_fun
      procedure :: cmp_D  => cmp_DD2_sgm_fun
      procedure :: new_D  => new_DD2_sgm_fun
      procedure :: upd_D  => upd_DD2_sgm_fun
      procedure :: cmp_dl => cmp_dl_sgm_fun
      procedure :: cmp_da => cmp_da_sgm_fun
   end type sgmfun_t
   private :: ini_sgm_fun, cmp_sgm_fun, ratio_sgm_fun, upd_sgm_fun, cmp_DD2_sgm_fun, &
   & new_DD2_sgm_fun, upd_DD2_sgm_fun, cmp_dl_sgm_fun, cmp_da_sgm_fun
contains
   subroutine ini_sgm_fun( sgmfun_t_in, n_fe, n_fe_u, dl_opt )
      class(sgmfun_t), intent(inout) :: sgmfun_t_in
      integer(int32),  intent(in)    :: n_fe, n_fe_u
      logical,         intent(in)    :: dl_opt
      allocate( sgmfun_t_in%det_G(1:n_wlk_max) ) ; sgmfun_t_in%det_G = 1.0_dp
      allocate( sgmfun_t_in%g(1:n_wlk_max) ) ; sgmfun_t_in%g = 1.0_dp
      allocate( sgmfun_t_in%n_mat_chngs(1:n_wlk_max) ) ; sgmfun_t_in%n_mat_chngs = 0_int32
      allocate( sgmfun_t_in%mat_G(1:n_fe_u,1:n_fe_u,1:n_wlk_max) )     ; sgmfun_t_in%mat_G = 0.0_dp
      allocate( sgmfun_t_in%inv_G(1:n_fe_u,1:n_fe_u,1:n_wlk_max) )     ; sgmfun_t_in%inv_G = 0.0_dp
      allocate( sgmfun_t_in%mat_L(1:n_forbs,1:2*n_fe_u,1:n_wlk_max) )  ; sgmfun_t_in%mat_L = 0.0_dp
      allocate( sgmfun_t_in%mat_D(1:n_forbs,1:2*n_fe_u,1:n_wlk_max) )  ; sgmfun_t_in%mat_D = 0.0_dp
      allocate( sgmfun_t_in%vec_G(1:n_fe_u,1:n_wlk_max) )              ; sgmfun_t_in%vec_G = 0.0_dp
      allocate( sgmfun_t_in%vec_Gig(1:n_fe_u,1:n_wlk_max) )            ; sgmfun_t_in%vec_Gig = 0.0_dp
      allocate( sgmfun_t_in%vec_Gi(1:n_fe_u,1:n_wlk_max) )             ; sgmfun_t_in%vec_Gi = 0.0_dp
      allocate( sgmfun_t_in%vec_D(1:n_forbs,1:n_wlk_max) )             ; sgmfun_t_in%vec_D = 0.0_dp
      allocate( sgmfun_t_in%vec_Dg(1:n_forbs,1:n_wlk_max) )            ; sgmfun_t_in%vec_Dg = 0.0_dp
      if ( dl_opt ) then
         allocate( sgmfun_t_in%mat_K(1:n_forbs,1:n_fe_u,1:n_wlk_max) )  ; sgmfun_t_in%mat_K = 0.0_dp
         allocate( sgmfun_t_in%mat_H(1:n_forbs,1:n_forbs,1:n_wlk_max) ) ; sgmfun_t_in%mat_H = 0.0_dp
      endif
      if ( aoc_opt ) then
         allocate( sgmfun_t_in%da_ln_det_G(1:n_par_forbs,1:n_wlk_max) )  ; sgmfun_t_in%da_ln_det_G = 0.0_dp
      endif
      allocate( sgmfun_t_in%DD2_ln_det_G(1:4*n_fe,1:n_wlk_max) )      ; sgmfun_t_in%DD2_ln_det_G = 0.0_dp
      if ( grad_updat ) then
         allocate( sgmfun_t_in%DD2_ln_det_G_new(1:4,1:n_wlk_max) )       ; sgmfun_t_in%DD2_ln_det_G_new = 0.0_dp
      endif
   end subroutine ini_sgm_fun
   subroutine cmp_sgm_fun( sgmfun_t_in, iw, n_fe, n_fe_u, n_fe_d, n_fe_s, sgm_c, mat_o )
      class(sgmfun_t), intent(inout) :: sgmfun_t_in
      integer(int32),                    intent(in) :: iw
      integer(int32),                    intent(in) :: n_fe
      integer(int32),                    intent(in) :: n_fe_u
      integer(int32),                    intent(in) :: n_fe_d
      integer(int32),                    intent(in) :: n_fe_s
      type(sgmcfs_t),                    intent(in) :: sgm_c
      real(dp), dimension(n_forbs,n_fe), intent(in) :: mat_o
      call cmp_L_mat( n_fe, n_fe_u, n_fe_d, n_fe_s, sgm_c%mat_Lmbd, mat_o, sgmfun_t_in%mat_L(:,:,iw) )
      if( n_fe_u.eq.1_int32 ) then
         sgmfun_t_in%mat_G(1,1,iw) = dot_product( mat_o(:,1), sgmfun_t_in%mat_L(:,2,iw) )
         sgmfun_t_in%det_G(iw) = sgmfun_t_in%mat_G(1,1,iw)
         sgmfun_t_in%inv_G(1,1,iw) = 1.0_dp / sgmfun_t_in%mat_G(1,1,iw)
      else
         call cmp_G_mat( n_fe_u, mat_o(:,1:n_fe_u), sgmfun_t_in%mat_L(:,n_fe_u+1:2*n_fe_u,iw), sgmfun_t_in%mat_G(:,:,iw) )
         call lu_dec_det_inv( n_fe_u, sgmfun_t_in%mat_G(:,:,iw), sgmfun_t_in%det_G(iw), sgmfun_t_in%inv_G(:,:,iw) )
      endif
      sgmfun_t_in%n_mat_chngs(iw) = 0_int32
   end subroutine cmp_sgm_fun
   subroutine ratio_sgm_fun( sgmfun_t_in, iw, i_fe, n_fe_u, vec_o_new, g )
      class(sgmfun_t), intent(inout) :: sgmfun_t_in
      integer(int32),               intent(in) :: iw
      integer(int32),               intent(in)  :: i_fe
      integer(int32),               intent(in)  :: n_fe_u
      real(dp), dimension(n_forbs), intent(in)  :: vec_o_new
      real(dp),                     intent(out) :: g
      if(n_fe_u.eq.1_int32) then
         if(i_fe.le.n_fe_u) then
            sgmfun_t_in%vec_G(1,iw) = dot_product( sgmfun_t_in%mat_L(:,2,iw), vec_o_new )
         else
            sgmfun_t_in%vec_G(1,iw) = dot_product( sgmfun_t_in%mat_L(:,1,iw), vec_o_new )
         endif
         sgmfun_t_in%g(iw) = sgmfun_t_in%vec_G(1,iw) / sgmfun_t_in%det_G(iw)
      else
         if(i_fe.le.n_fe_u) then
            call cmp_G_vec( n_fe_u, sgmfun_t_in%mat_L(:,n_fe_u+1:2*n_fe_u,iw), vec_o_new, sgmfun_t_in%vec_G(:,iw) )
            sgmfun_t_in%vec_Gi(:,iw) = sgmfun_t_in%inv_G(:,i_fe,iw)
         else
            call cmp_G_vec( n_fe_u, sgmfun_t_in%mat_L(:,1:n_fe_u,iw), vec_o_new, sgmfun_t_in%vec_G(:,iw) )
            sgmfun_t_in%vec_Gi(:,iw) = sgmfun_t_in%inv_G(i_fe-n_fe_u,:,iw)
         endif
         sgmfun_t_in%g(iw) = dot_product( sgmfun_t_in%vec_Gi(:,iw), sgmfun_t_in%vec_G(:,iw) )
      endif
      g = sgmfun_t_in%g(iw)
   end subroutine ratio_sgm_fun
   subroutine upd_sgm_fun( sgmfun_t_in, iw, i_fe, n_fe, n_fe_u, sgm_c, vec_o_new )
      class(sgmfun_t), intent(inout) :: sgmfun_t_in
      integer(int32),               intent(in) :: iw
      integer(int32),               intent(in) :: i_fe
      integer(int32),               intent(in) :: n_fe
      integer(int32),               intent(in) :: n_fe_u
      type(sgmcfs_t),               intent(in) :: sgm_c
      real(dp), dimension(n_forbs), intent(in) :: vec_o_new
      sgmfun_t_in%det_G(iw) = sgmfun_t_in%det_G(iw) * sgmfun_t_in%g(iw)
      call upd_L_mat( i_fe, n_fe_u, sgm_c%mat_Lmbd(:,1:n_forbs), vec_o_new, sgmfun_t_in%mat_L(:,:,iw) )
      if(n_fe_u.eq.1_int32) then
         sgmfun_t_in%mat_G(1,1,iw) = sgmfun_t_in%vec_G(1,iw)
         sgmfun_t_in%inv_G(1,1,iw) = 1.0_dp / sgmfun_t_in%vec_G(1,iw)
      else
         if(i_fe.le.n_fe_u) then
            sgmfun_t_in%mat_G(i_fe,:,iw) = sgmfun_t_in%vec_G(:,iw)
         else
            sgmfun_t_in%mat_G(:,i_fe-n_fe_u,iw) = sgmfun_t_in%vec_G(:,iw)
         endif
         if( sgmfun_t_in%n_mat_chngs(iw).eq.n_fe .or.abs(sgmfun_t_in%g(iw)).lt.10d-14) then
            call lu_dec_det_inv( n_fe_u, sgmfun_t_in%mat_G(:,:,iw), sgmfun_t_in%det_G(iw), sgmfun_t_in%inv_G(:,:,iw) )
            sgmfun_t_in%n_mat_chngs(iw) = 0_int32
         else
            call upd_G_inv( i_fe, n_fe_u, sgmfun_t_in%g(iw), sgmfun_t_in%vec_G(:,iw), sgmfun_t_in%vec_Gi(:,iw), sgmfun_t_in%vec_Gig(:,iw), sgmfun_t_in%inv_G(:,:,iw) )
            sgmfun_t_in%n_mat_chngs(iw) = sgmfun_t_in%n_mat_chngs(iw) + 1_int32
         endif
      endif
   end subroutine upd_sgm_fun
   subroutine cmp_DD2_sgm_fun( sgmfun_t_in, iw, n_fe, n_fe_u, DD2_o )
      class(sgmfun_t), intent(inout) :: sgmfun_t_in
      integer(int32),                      intent(in) :: iw
      integer(int32),                      intent(in) :: n_fe
      integer(int32),                      intent(in) :: n_fe_u
      real(dp), dimension(n_forbs,4*n_fe), intent(in) :: DD2_o
      if( n_fe_u .eq. 1_int32 ) then
         sgmfun_t_in%mat_D(:,1,iw) = sgmfun_t_in%inv_G(1,1,iw) * sgmfun_t_in%mat_L(:,2,iw)
         sgmfun_t_in%mat_D(:,2,iw) = sgmfun_t_in%inv_G(1,1,iw) * sgmfun_t_in%mat_L(:,1,iw)
      else
         call cmp_D_mat( n_fe_u, sgmfun_t_in%mat_L(:,:,iw), sgmfun_t_in%inv_G(:,:,iw), sgmfun_t_in%mat_D(:,:,iw) )
      endif
      call cmp_DD2_vec( n_fe, n_fe_u, DD2_o, sgmfun_t_in%mat_D(:,:,iw), sgmfun_t_in%DD2_ln_det_G(:,iw) )
   end subroutine cmp_DD2_sgm_fun
   subroutine new_DD2_sgm_fun( sgmfun_t_in, iw, i_fe, vec_DD2_o_new )
      class(sgmfun_t), intent(inout) :: sgmfun_t_in
      integer(int32),                 intent(in) :: iw
      integer(int32),                 intent(in) :: i_fe
      real(dp), dimension(n_forbs,4), intent(in) :: vec_DD2_o_new
      call dgemv('T', n_forbs, 4, 1.0_dp/sgmfun_t_in%g(iw), vec_DD2_o_new, n_forbs, &
      & sgmfun_t_in%mat_D(:,i_fe,iw), 1_int32, 0.0_dp, sgmfun_t_in%DD2_ln_det_G_new(:,iw), 1_int32 )
      sgmfun_t_in%DD2_ln_det_G_new(4,iw) = sgmfun_t_in%DD2_ln_det_G_new(4,iw) - sum(sgmfun_t_in%DD2_ln_det_G_new(1:3,iw)**2)
   end subroutine new_DD2_sgm_fun
   subroutine upd_DD2_sgm_fun( sgmfun_t_in, iw, i_fe, n_fe, n_fe_u, DD2_o )
      class(sgmfun_t), intent(inout) :: sgmfun_t_in
      integer(int32),                      intent(in) :: iw
      integer(int32),                      intent(in) :: i_fe
      integer(int32),                      intent(in) :: n_fe
      integer(int32),                      intent(in) :: n_fe_u
      real(dp), dimension(n_forbs,4*n_fe), intent(in) :: DD2_o
      if( n_fe_u .eq. 1_int32 ) then
         sgmfun_t_in%mat_D(:,1,iw) = sgmfun_t_in%inv_G(1,1,iw) * sgmfun_t_in%mat_L(:,2,iw)
         sgmfun_t_in%mat_D(:,2,iw) = sgmfun_t_in%inv_G(1,1,iw) * sgmfun_t_in%mat_L(:,1,iw)
      else
         if( sgmfun_t_in%n_mat_chngs(iw).eq.0_int32 ) then
            call cmp_D_mat( n_fe_u, sgmfun_t_in%mat_L(:,:,iw), sgmfun_t_in%inv_G(:,:,iw), sgmfun_t_in%mat_D(:,:,iw) )
         else
            call upd_D_mat( i_fe, n_fe_u, sgmfun_t_in%g(iw), sgmfun_t_in%mat_L(:,:,iw), sgmfun_t_in%vec_G(:,iw), sgmfun_t_in%vec_Gi(:,iw),&
            & sgmfun_t_in%vec_Gig(:,iw), sgmfun_t_in%vec_Dg(:,iw), sgmfun_t_in%vec_D(:,iw), sgmfun_t_in%mat_D(:,:,iw) )
         endif
      endif
      call cmp_DD2_vec( n_fe, n_fe_u, DD2_o, sgmfun_t_in%mat_D(:,:,iw), sgmfun_t_in%DD2_ln_det_G(:,iw) )
   end subroutine upd_DD2_sgm_fun
   subroutine cmp_dl_sgm_fun( sgmfun_t_in, iw, n_fe, n_fe_u, n_fe_d, n_fe_s, sgm_c, &
      & mat_o, n_par_det, dl_ln_det_G )
      class(sgmfun_t), intent(inout) :: sgmfun_t_in
      integer(int32),                    intent(in)    :: iw
      integer(int32),                    intent(in)    :: n_fe
      integer(int32),                    intent(in)    :: n_fe_u
      integer(int32),                    intent(in)    :: n_fe_d
      integer(int32),                    intent(in)    :: n_fe_s
      type(sgmcfs_t),                    intent(in)    :: sgm_c
      real(dp), dimension(n_forbs,n_fe), intent(in)    :: mat_o
      integer(int32),                    intent(in)    :: n_par_det
      real(dp), dimension(n_par_det),    intent(inout) :: dl_ln_det_G
      integer(int32) :: i1
      if( n_fe_u .eq. 1_int32 ) then
         sgmfun_t_in%mat_K(:,1,iw) = mat_o(:,1) * sgmfun_t_in%inv_G(1,1,iw)
         sgmfun_t_in%mat_H(:,:,iw) = 0.0_dp
         call dger( n_forbs, n_forbs, 1.0_dp, mat_o(:,2), 1_int32, &
         & sgmfun_t_in%mat_K(:,1,iw), 1_int32, sgmfun_t_in%mat_H(:,:,iw), n_forbs )
      else
         call cmp_KH_mat( n_fe, n_fe_u, n_fe_d, mat_o, sgmfun_t_in%inv_G(:,:,iw), sgmfun_t_in%mat_K(:,:,iw), sgmfun_t_in%mat_H(:,:,iw) )
      endif
      do i1 = 1_int32, sgm_c%n_nz_Lmbd
         dl_ln_det_G(i1) = sgmfun_t_in%mat_H(sgm_c%nz_Lmbd(i1)%o1,sgm_c%nz_Lmbd(i1)%o2,iw)
      enddo
      if ( n_fe_s.gt.0_int32 ) then
         do i1 = 1_int32, sgm_c%n_nz_Lmbd_u
            dl_ln_det_G(sgm_c%n_nz_Lmbd+i1) = sgmfun_t_in%mat_K(sgm_c%nz_Lmbd_u(i1)%o1,n_fe_d+sgm_c%nz_Lmbd_u(i1)%o2-n_forbs,iw)
         enddo
      endif 
   end subroutine cmp_dl_sgm_fun
   subroutine cmp_da_sgm_fun( sgmfun_t_in, iw, n_fe, mat_da_o, da_ln_det_G )
      class(sgmfun_t), intent(inout) :: sgmfun_t_in
      integer(int32),                        intent(in)    :: iw
      integer(int32),                        intent(in)    :: n_fe
      real(dp), dimension(n_par_forbs,n_fe), intent(in)    :: mat_da_o
      real(dp), dimension(n_par_forbs),      intent(inout) :: da_ln_det_G
      integer(int32) :: i1, i2, i3
      integer(int32) :: io, pi, pf
      sgmfun_t_in%da_ln_det_G(:,iw) = 0.0_dp
      do i1 = 1_int32, n_fe ; do i2 = 1_int32, n_at ; do i3 = 1_int32, forbs(i2)%n_orbs
               io = forbs(i2)%orb(i3)%i_orb
               if( forbs(i2)%orb(i3)%n_par.gt.0_int32 ) then
                  pi = forbs(i2)%orb(i3)%i_par
                  pf = pi + forbs(i2)%orb(i3)%n_par - 1_int32
                  sgmfun_t_in%da_ln_det_G(pi:pf,iw) = sgmfun_t_in%da_ln_det_G(pi:pf,iw) + sgmfun_t_in%mat_D(io,i1,iw) * mat_da_o(pi:pf,i1)
               endif
            enddo; enddo ; enddo 
      da_ln_det_G(:) = da_ln_det_G(:) + sgmfun_t_in%da_ln_det_G(:,iw)
   end subroutine cmp_da_sgm_fun
end module singlet_geminal_c
