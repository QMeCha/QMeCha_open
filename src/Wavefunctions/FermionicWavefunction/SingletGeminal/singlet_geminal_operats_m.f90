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
module singlet_geminal_operats_m
   use fortran_kinds_v,          only: int32, dp
   use fermionic_orbitals_m,     only: n_forbs
   use fermionic_wavefunction_v, only: spin
   implicit none
   public :: cmp_L_mat, upd_L_mat, cmp_G_mat, cmp_G_vec, cmp_D_mat, upd_G_inv, &
   & cmp_KH_mat, upd_D_mat, cmp_DD2_vec
contains
   subroutine cmp_L_mat( n_fe, n_fe_u, n_fe_d, n_fe_s, mat_Lmbd, mat_o, mat_L )
      integer(int32),                              intent(in)    :: n_fe
      integer(int32),                              intent(in)    :: n_fe_u
      integer(int32),                              intent(in)    :: n_fe_d
      integer(int32),                              intent(in)    :: n_fe_s
      real(dp), dimension(n_forbs,n_forbs+n_fe_s), intent(in)    :: mat_Lmbd
      real(dp), dimension(n_forbs,n_fe),           intent(in)    :: mat_o
      real(dp), dimension(n_forbs,2*n_fe_u),       intent(inout) :: mat_L
      if( spin.eq.'R' ) then
         call dsymm('L', 'L', n_forbs, n_fe, 1.0_dp, mat_Lmbd(:,1:n_forbs), n_forbs, &
         & mat_o, n_forbs, 0.0_dp, mat_L(:, 1:n_fe), n_forbs )
      else
         call dgemm('N', 'N', n_forbs, n_fe_u, n_forbs, 1.0_dp, mat_Lmbd(:,1:n_forbs), n_forbs, &
         & mat_o(:,1:n_fe_u), n_forbs, 0.0_dp, mat_L(:,1:n_fe_u), n_forbs )
         call dgemm('T', 'N', n_forbs, n_fe_d, n_forbs, 1.0_dp, mat_Lmbd(:,1:n_forbs), n_forbs, &
         & mat_o(:,n_fe_u+1:n_fe), n_forbs, 0.0_dp, mat_L(:,n_fe_u+1:n_fe), n_forbs )
      endif
      if ( n_fe_s.gt.0_int32 ) then
         mat_L(:,n_fe+1:n_fe+n_fe_s) = mat_Lmbd(:,n_forbs+1:n_forbs+n_fe_s)
      endif
   end subroutine cmp_L_mat
   subroutine upd_L_mat( i_fe, n_fe_u, mat_Lmbd, vec_o_new, mat_L )
      integer(int32),                        intent(in)    :: i_fe
      integer(int32),                        intent(in)    :: n_fe_u
      real(dp), dimension(n_forbs,n_forbs),  intent(in)    :: mat_Lmbd
      real(dp), dimension(n_forbs),          intent(in)    :: vec_o_new
      real(dp), dimension(n_forbs,2*n_fe_u), intent(inout) :: mat_L
      if ( spin.eq.'R') then
         call dsymv( 'L', n_forbs, 1.0_dp, mat_Lmbd(:,1:n_forbs), n_forbs,&
         & vec_o_new, 1_int32, 0.0_dp, mat_L(:,i_fe), 1_int32 )
      else
         if( i_fe.le.n_fe_u ) then
            call dgemv( 'N', n_forbs, n_forbs, 1.0_dp, mat_Lmbd(:,1:n_forbs), n_forbs,&
            & vec_o_new, 1_int32, 0.0_dp, mat_L(:,i_fe), 1_int32 )
         else
            call dgemv( 'T', n_forbs, n_forbs, 1.0_dp, mat_Lmbd(:,1:n_forbs), n_forbs,&
            & vec_o_new, 1_int32, 0.0_dp, mat_L(:,i_fe), 1_int32 )
         endif
      endif
   end subroutine upd_L_mat
   subroutine cmp_G_mat( n_fe_u, mat_o, mat_L, mat_G )
      integer(int32),                      intent(in)    :: n_fe_u
      real(dp), dimension(n_forbs,n_fe_u), intent(in)    :: mat_o
      real(dp), dimension(n_forbs,n_fe_u), intent(in)    :: mat_L
      real(dp), dimension(n_fe_u,n_fe_u),  intent(inout) :: mat_G
      call dgemm('T','N', n_fe_u, n_fe_u, n_forbs, 1.0_dp, mat_o, n_forbs, &
      & mat_L, n_forbs, 0.0_dp, mat_G, n_fe_u )
   end subroutine cmp_G_mat
   subroutine cmp_G_vec( n_fe_u, mat_L, vec_o_new, vec_G_new )
      integer(int32),                      intent(in)    :: n_fe_u
      real(dp), dimension(n_forbs,n_fe_u), intent(in)    :: mat_L
      real(dp), dimension(n_forbs),        intent(in)    :: vec_o_new
      real(dp), dimension(n_fe_u,n_fe_u),  intent(inout) :: vec_G_new
      call dgemv('T', n_forbs, n_fe_u, 1.0_dp, mat_L, n_forbs,&
      & vec_o_new, 1_int32, 0.0_dp, vec_G_new, 1_int32 )
   end subroutine cmp_G_vec
   subroutine cmp_D_mat( n_fe_u, mat_L, inv_G, mat_D )
      integer(int32),                            intent(in)    :: n_fe_u
      real(dp), dimension(n_forbs,2*n_fe_u), intent(in)    :: mat_L
      real(dp), dimension(n_fe_u,n_fe_u),    intent(in)    :: inv_G
      real(dp), dimension(n_forbs,2*n_fe_u), intent(inout) :: mat_D
      call dgemm('N', 'N', n_forbs, n_fe_u, n_fe_u, 1.0_dp, mat_L(:,n_fe_u+1:2*n_fe_u), n_forbs, &
      & inv_G, n_fe_u, 0.0_dp, mat_D(:,1:n_fe_u), n_forbs )
      call dgemm('N', 'T', n_forbs, n_fe_u, n_fe_u, 1.0_dp, mat_L(:,1:n_fe_u), n_forbs, &
      & inv_G, n_fe_u, 0.0_dp, mat_D(:,n_fe_u+1:2*n_fe_u), n_forbs )
   end subroutine cmp_D_mat
   subroutine upd_G_inv( i_fe, n_fe_u, g, vec_G, vec_Gi, vec_Gig, inv_G )
      integer(int32),                     intent(in)    :: i_fe
      integer(int32),                     intent(in)    :: n_fe_u
      real(dp),                           intent(in)    :: g
      real(dp), dimension(n_fe_u),        intent(in)    :: vec_G
      real(dp), dimension(n_fe_u),        intent(in)    :: vec_Gi
      real(dp), dimension(n_fe_u),        intent(inout) :: vec_Gig
      real(dp), dimension(n_fe_u,n_fe_u), intent(inout) :: inv_G
      integer(int32) :: i_d
      if( i_fe.le.n_fe_u ) then
         i_d = i_fe
         call sm_upd_inv('R', n_fe_u, i_d, vec_G, vec_Gi, vec_Gig, g, inv_G)
      else
         i_d = i_fe - n_fe_u
         call sm_upd_inv('C', n_fe_u, i_d, vec_G, vec_Gi, vec_Gig, g, inv_G)
      endif
   end subroutine upd_G_inv
   subroutine upd_D_mat( i_fe, n_fe_u, g, mat_L, vec_G, vec_Gi, vec_Gig, vec_Dg, vec_D, mat_D )      
      integer(int32),                        intent(in)    :: i_fe
      integer(int32),                        intent(in)    :: n_fe_u
      real(dp),                              intent(in)    :: g
      real(dp), dimension(n_forbs,2*n_fe_u), intent(in)    :: mat_L
      real(dp), dimension(n_fe_u),           intent(in)    :: vec_G
      real(dp), dimension(n_fe_u),           intent(in)    :: vec_Gi
      real(dp), dimension(n_fe_u),           intent(in)    :: vec_Gig
      real(dp), dimension(n_forbs),          intent(inout) :: vec_Dg
      real(dp), dimension(n_forbs),          intent(inout) :: vec_D
      real(dp), dimension(n_forbs,2*n_fe_u), intent(inout) :: mat_D
      vec_D = mat_D(:,i_fe)
      vec_Dg = - mat_L(:,i_fe)
      if( i_fe.le.n_fe_u ) then
         call dger ( n_forbs, n_fe_u, -1.0_dp / g, vec_D, 1_int32,&
         & vec_Gig, 1_int32, mat_D(:,1:n_fe_u), n_forbs )
         call dgemv('N', n_forbs, n_fe_u, 1.0_dp, mat_D(:,n_fe_u+1:2*n_fe_u), n_forbs, &
         & vec_G, 1_int32, 1.0_dp, vec_Dg, 1_int32 )
         call dger ( n_forbs, n_fe_u, -1.0_dp/g, vec_Dg, 1_int32,&
         & vec_Gi, 1_int32, mat_D(:,n_fe_u+1:2*n_fe_u), n_forbs )
      else
         call dgemv('N', n_forbs, n_fe_u, 1.0_dp, mat_D(:,1:n_fe_u), n_forbs,&
         & vec_G, 1_int32, 1.0_dp, vec_Dg, 1_int32 )
         call dger ( n_forbs, n_fe_u, -1.0_dp/ g, vec_Dg, 1_int32,&
         & vec_Gi, 1_int32, mat_D(:,1:n_fe_u), n_forbs )
         call dger ( n_forbs, n_fe_u, -1.0_dp/g, vec_D, 1_int32,&
         & vec_Gig, 1_int32, mat_D(:,n_fe_u+1:2*n_fe_u), n_forbs )
      endif
   end subroutine upd_D_mat
   subroutine cmp_DD2_vec( n_fe, n_fe_u, DD2_o, mat_D, DD2_ln_det_G )
      integer(int32),                        intent(in)    :: n_fe
      integer(int32),                        intent(in)    :: n_fe_u
      real(dp), dimension(n_forbs,4*n_fe),   intent(in)    :: DD2_o
      real(dp), dimension(n_forbs,2*n_fe_u), intent(in)    :: mat_D
      real(dp), dimension(4*n_fe),           intent(inout) :: DD2_ln_det_G
      integer(int32) :: i1
      do i1 = 1_int32, n_fe
         call dgemv('T', n_forbs, 4, 1.0_dp, DD2_o(:,4*i1-3:4*i1), n_forbs, &
         & mat_D(:,i1), 1_int32, 0.0_dp, DD2_ln_det_G(4*i1-3:4*i1), 1_int32)
         DD2_ln_det_G(4*i1) = DD2_ln_det_G(4*i1) - sum(DD2_ln_det_G(4*i1-3:4*i1-1)**2)
      enddo
   end subroutine cmp_DD2_vec
   subroutine cmp_KH_mat( n_fe, n_fe_u, n_fe_d, mat_o, inv_G, mat_K, mat_H )
      integer(int32),                       intent(in)    :: n_fe
      integer(int32),                       intent(in)    :: n_fe_u
      integer(int32),                       intent(in)    :: n_fe_d
      real(dp), dimension(n_forbs,n_fe),    intent(in)    :: mat_o
      real(dp), dimension(n_fe_u,n_fe_u),   intent(in)    :: inv_G
      real(dp), dimension(n_forbs,n_fe_u),  intent(inout) :: mat_K
      real(dp), dimension(n_forbs,n_forbs), intent(inout) :: mat_H
      call dgemm('N', 'T', n_forbs, n_fe_u, n_fe_u, 1.0_dp, mat_o(:,1:n_fe_u), n_forbs, &
      & inv_G, n_fe_u, 0.0_dp, mat_K(:,1:n_fe_u), n_forbs )
      call dgemm('N', 'T', n_forbs, n_forbs, n_fe_d, 1.0_dp, mat_o(:,n_fe_u+1:n_fe), n_forbs, &
      & mat_K(:,1:n_fe_d), n_forbs, 0.0_dp, mat_H, n_forbs )
   end subroutine cmp_KH_mat
end module singlet_geminal_operats_m
