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
module pfaffian_operats_m
   use fortran_kinds_v, only: int32, dp
   use pfaffian_params_m, only: pffcfs_t
   use fermionic_orbitals_m, only: n_forbs
   use fermionic_wavefunction_v, only: spin
   implicit none
   public :: cmp_L_mat, cmp_P_mat, cmp_P_vec, cmp_D_mat, upd_L_mat, upd_P_inv, &
   & cmp_Lp_mat, upd_D_mat, cmp_DD2_vec, cmp_H_mat
contains
   subroutine cmp_L_mat( pff_c, n_fe, n_fe_u, n_fe_d, d_P, mat_o, mat_L_u, mat_L_d )
      type(pffcfs_t),                    intent(in)    :: pff_c
      integer(int32),                    intent(in)    :: n_fe
      integer(int32),                    intent(in)    :: n_fe_u
      integer(int32),                    intent(in)    :: n_fe_d
      integer(int32),                    intent(in)    :: d_P
      real(dp), dimension(n_forbs,n_fe), intent(in)    :: mat_o
      real(dp), dimension(n_forbs,d_P),  intent(inout) :: mat_L_u
      real(dp), dimension(n_forbs,d_P),  intent(inout) :: mat_L_d
      call dgemm('T', 'N', n_forbs, n_fe_u, n_forbs, 1.0_dp, pff_c%mat_Z_u, n_forbs, &
      & mat_o(:,1:n_fe_u), n_forbs, 0.0_dp, mat_L_u(:,1:n_fe_u), n_forbs )
      if( spin.eq.'R' ) then
         call dsymm('L', 'L', n_forbs, n_fe_d, 1.0_dp, pff_c%mat_Lmbd(:,1:n_forbs), n_forbs, &
         & mat_o(:,n_fe_u+1:n_fe), n_forbs, 0.0_dp, mat_L_u(:,n_fe_u+1:n_fe), n_forbs )
      else
         call dgemm('T', 'N', n_forbs, n_fe_d, n_forbs, 1.0_dp, pff_c%mat_Lmbd(:,1:n_forbs), n_forbs, &
         & mat_o(:,n_fe_u+1:n_fe), n_forbs, 0.0_dp, mat_L_u(:,n_fe_u+1:n_fe), n_forbs )
      endif
      if( spin.eq.'R' ) then
         call dsymm('L', 'L', n_forbs, n_fe_u, -1.0_dp, pff_c%mat_Lmbd(:,1:n_forbs), n_forbs, &
         & mat_o(:,1:n_fe_u), n_forbs, 0.0_dp, mat_L_d(:,1:n_fe_u), n_forbs )
      else
         call dgemm('N', 'N', n_forbs, n_fe_u, n_forbs, -1.0_dp, pff_c%mat_Lmbd(:,1:n_forbs), n_forbs, &
         & mat_o(:,1:n_fe_u), n_forbs, 0.0_dp, mat_L_d(:,1:n_fe_u), n_forbs )
      endif
      if ( n_fe_d.gt.1_int32 ) then
         call dgemm('T', 'N', n_forbs, n_fe_d, n_forbs, 1.0_dp, pff_c%mat_Z_d, n_forbs, &
         & mat_o(:,n_fe_u+1:n_fe), n_forbs, 0.0_dp, mat_L_d(:,n_fe_u+1:n_fe), n_forbs )
      endif
      if( d_P.gt.n_fe ) then
         mat_L_u(:,d_P) = pff_c%vec_u(:)
         mat_L_d(:,d_P) = pff_c%vec_d(:)
      endif
   end subroutine cmp_L_mat
   subroutine cmp_P_mat( n_fe, n_fe_u, n_fe_d, d_P, mat_o, mat_L_u, mat_L_d, mat_P )
      integer(int32),                    intent(in)    :: n_fe
      integer(int32),                    intent(in)    :: n_fe_u
      integer(int32),                    intent(in)    :: n_fe_d
      integer(int32),                    intent(in)    :: d_P
      real(dp), dimension(n_forbs,n_fe), intent(in)    :: mat_o
      real(dp), dimension(n_forbs,d_P),  intent(in)    :: mat_L_u
      real(dp), dimension(n_forbs,d_P),  intent(in)    :: mat_L_d
      real(dp), dimension(d_P,d_P),      intent(inout) :: mat_P
      integer(int32) :: i1
      call dgemm('T', 'N', d_P, n_fe_u, n_forbs, -1.0_dp, mat_L_u, n_forbs, &
      & mat_o(:,1:n_fe_u), n_forbs, 0.0_dp, mat_P(:,1:n_fe_u), d_P )
      call dgemm('T', 'N', d_P, n_fe_d, n_forbs, -1.0_dp, mat_L_d, n_forbs, &
      & mat_o(:,n_fe_u+1:n_fe), n_forbs, 0.0_dp, mat_P(:,n_fe_u+1:n_fe), d_P )
      if (d_P.gt.n_fe) mat_P(:,d_P) = - mat_P(d_P,:)
      do i1 = 1_int32, d_P
         mat_P(i1,i1) = 0.0_dp
      enddo
   end subroutine cmp_P_mat
   subroutine cmp_P_vec( i_fe, d_P, o_vec, mat_L, vec_P )
      integer(int32),                   intent(in)    :: i_fe
      integer(int32),                   intent(in)    :: d_P
      real(dp), dimension(n_forbs),     intent(in)    :: o_vec
      real(dp), dimension(n_forbs,d_P), intent(in)    :: mat_L
      real(dp), dimension(d_P),         intent(inout) :: vec_P
      call dgemv( 'T', n_forbs, d_P, -1.0_dp, mat_L, n_forbs, o_vec, 1_int32, &
      & 0.0_dp, vec_P, 1_int32 )
      vec_P(i_fe) = 0.0_dp
   end subroutine cmp_P_vec
   subroutine cmp_D_mat( n_fe, n_fe_u, n_fe_d, d_P, mat_L_u, mat_L_d, inv_P, mat_D )
      integer(int32),                    intent(in)    :: n_fe
      integer(int32),                    intent(in)    :: n_fe_u
      integer(int32),                    intent(in)    :: n_fe_d
      integer(int32),                    intent(in)    :: d_P
      real(dp), dimension(n_forbs,d_P),  intent(in)    :: mat_L_u
      real(dp), dimension(n_forbs,d_P),  intent(in)    :: mat_L_d
      real(dp), dimension(d_P,d_P),      intent(in)    :: inv_P
      real(dp), dimension(n_forbs,n_fe), intent(inout) :: mat_D
      call dgemm('N', 'N', n_forbs, n_fe_u, d_P, 1.0_dp, mat_L_u, n_forbs, &
      & inv_P(:,1:n_fe_u), d_P, 0.0_dp, mat_D(:,1:n_fe_u), n_forbs )
      call dgemm('N', 'N', n_forbs, n_fe_d, d_P, 1.0_dp, mat_L_d, n_forbs, &
      & inv_P(:,n_fe_u+1:n_fe), d_P, 0.0_dp, mat_D(:,n_fe_u+1:n_fe), n_forbs )
   end subroutine cmp_D_mat
   subroutine upd_L_mat( pff_c, i_fe, n_fe, n_fe_u, n_fe_d, d_P, mat_o, vec_dl, mat_L_u, mat_L_d )
      type(pffcfs_t),                    intent(in)    :: pff_c
      integer(int32),                    intent(in)    :: i_fe
      integer(int32),                    intent(in)    :: n_fe
      integer(int32),                    intent(in)    :: n_fe_u
      integer(int32),                    intent(in)    :: n_fe_d
      integer(int32),                    intent(in)    :: d_P
      real(dp), dimension(n_forbs,n_fe), intent(in)    :: mat_o
      real(dp), dimension(2*n_forbs),    intent(inout) :: vec_dl
      real(dp), dimension(n_forbs,d_P),  intent(inout) :: mat_L_u
      real(dp), dimension(n_forbs,d_P),  intent(inout) :: mat_L_d
      vec_dl(1:n_forbs)           = mat_L_u(:,i_fe)
      vec_dl(n_forbs+1:2*n_forbs) = mat_L_d(:,i_fe)
      if ( i_fe.le.n_fe_u ) then
         call dgemv('T', n_forbs, n_forbs, 1.0_dp, pff_c%mat_Z_u, n_forbs, &
         & mat_o(:,i_fe), 1_int32, 0.0_dp, mat_L_u(:,i_fe), 1_int32  )
         if( spin.eq.'R' ) then
            call dsymv('L', n_forbs, -1.0_dp, pff_c%mat_Lmbd(:,1:n_forbs), n_forbs, &
            & mat_o(:,i_fe), 1_int32, 0.0_dp, mat_L_d(:,i_fe), 1_int32  )
         else
            call dgemv('N', n_forbs, n_forbs, -1.0_dp, pff_c%mat_Lmbd(:,1:n_forbs), n_forbs, &
            & mat_o(:,i_fe), 1_int32, 0.0_dp, mat_L_d(:,i_fe), 1_int32  )
         endif
      else
         if( spin.eq.'R' ) then
            call dsymv('L', n_forbs, 1.0_dp, pff_c%mat_Lmbd(:,1:n_forbs), n_forbs, &
            & mat_o(:,i_fe), 1_int32, 0.0_dp, mat_L_u(:,i_fe), 1_int32  )
         else
            call dgemv('T', n_forbs, n_forbs, 1.0_dp, pff_c%mat_Lmbd(:,1:n_forbs), n_forbs, &
            & mat_o(:,i_fe), 1_int32, 0.0_dp, mat_L_u(:,i_fe), 1_int32  )
         endif
         if ( n_fe_d.gt.1_int32 ) then
            call dgemv('T', n_forbs, n_forbs, 1.0_dp, pff_c%mat_Z_d, n_forbs, &
            & mat_o(:,i_fe), 1_int32, 0.0_dp, mat_L_d(:,i_fe), 1_int32  )
         endif
      endif
      vec_dl(1:n_forbs)           = vec_dl(1:n_forbs)           - mat_L_u(:,i_fe)
      vec_dl(n_forbs+1:2*n_forbs) = vec_dl(n_forbs+1:2*n_forbs) - mat_L_d(:,i_fe)
   end subroutine upd_L_mat
   subroutine upd_P_inv( i_fe, d_P, g, vec_P, vec_Pi, vec_Pip, inv_P )
      integer(int32),               intent(in)    :: i_fe
      integer(int32),               intent(in)    :: d_P
      real(dp),                     intent(in)    :: g
      real(dp), dimension(d_P),     intent(in)    :: vec_P
      real(dp), dimension(d_P),     intent(in)    :: vec_Pi
      real(dp), dimension(d_P),     intent(inout) :: vec_Pip
      real(dp), dimension(d_P,d_P), intent(inout) :: inv_P
      integer(int32) :: i1
      call dgemv( 'T', d_P, d_P, -1.0_dp, inv_P, d_P, vec_P(1:d_P), 1_int32, &
      & 0.0_dp, vec_Pip(1:d_P), 1_int32)
      vec_Pip(i_fe) = vec_Pip(i_fe) - 1.0_dp
      call dger( d_P, d_P, -1.0_dp/g, vec_Pip(1:d_P), 1_int32, vec_Pi(1:d_P), 1_int32, &
      &  inv_P, d_P )
      call dger( d_P, d_P, +1.0_dp/g, vec_Pi(1:d_P), 1_int32, vec_Pip(1:d_P), 1_int32, &
      &  inv_P, d_P )
      do i1 = 1_int32, d_P
         inv_P(i1,i1) = 0.0_dp
      enddo
   end subroutine upd_P_inv
   subroutine upd_D_mat( n_fe, n_fe_u, n_fe_d, d_P, mat_L_u, mat_L_d, g, &
      & vec_Pi, vec_Pip, vec_LPi, vec_LPip, mat_D )
      integer(int32),                    intent(in)    :: n_fe
      integer(int32),                    intent(in)    :: n_fe_u
      integer(int32),                    intent(in)    :: n_fe_d
      integer(int32),                    intent(in)    :: d_P
      real(dp), dimension(n_forbs,d_P),  intent(in)    :: mat_L_u
      real(dp), dimension(n_forbs,d_P),  intent(in)    :: mat_L_d
      real(dp),                          intent(in)    :: g
      real(dp), dimension(d_P),          intent(in)    :: vec_Pi
      real(dp), dimension(d_P),          intent(in)    :: vec_Pip
      real(dp), dimension(2*n_forbs),    intent(inout) :: vec_LPi
      real(dp), dimension(2*n_forbs),    intent(inout) :: vec_LPip
      real(dp), dimension(n_forbs,n_fe), intent(inout) :: mat_D
      call dgemv( 'N', n_forbs, d_P, 1.0_dp, mat_L_u, n_forbs, vec_Pip(1:d_P), 1_int32, &
      & g, vec_LPip(1:n_forbs), 1_int32)
      call dgemv( 'N', n_forbs, d_P, 1.0_dp, mat_L_d, n_forbs, vec_Pip(1:d_P), 1_int32, &
      & g, vec_LPip(n_forbs+1:2*n_forbs), 1_int32)
      call dgemv( 'N', n_forbs, d_P, 1.0_dp, mat_L_u, n_forbs, vec_Pi(1:d_P), 1_int32, &
      & 0.0_dp, vec_LPi(1:n_forbs), 1_int32)
      call dgemv( 'N', n_forbs, d_P, 1.0_dp, mat_L_d, n_forbs, vec_Pi(1:d_P), 1_int32, &
      & 0.0_dp, vec_LPi(n_forbs+1:2*n_forbs), 1_int32)
      call dger (n_forbs, n_fe_u, -1.0_dp / g, vec_LPip(1:n_forbs), 1_int32, &
      & vec_Pi(1:n_fe_u), 1_int32, mat_D(:,1:n_fe_u), n_forbs)
      call dger (n_forbs, n_fe_u, 1.0_dp / g, vec_LPi(1:n_forbs), 1_int32, &
      & vec_Pip(1:n_fe_u), 1_int32, mat_D(:,1:n_fe_u), n_forbs)
      call dger (n_forbs, n_fe_d, -1.0_dp / g, vec_LPip(n_forbs+1:2*n_forbs), 1_int32, &
      & vec_Pi(n_fe_u+1:n_fe), 1_int32, mat_D(:,n_fe_u+1:n_fe), n_forbs)
      call dger (n_forbs, n_fe_d, 1.0_dp / g, vec_LPi(n_forbs+1:2*n_forbs), 1_int32, &
      & vec_Pip(n_fe_u+1:n_fe), 1_int32, mat_D(:,n_fe_u+1:n_fe), n_forbs)
   end subroutine upd_D_mat
   subroutine cmp_DD2_vec( n_fe, mat_D, DD2_o, DD2_ln_pff_P )
      integer(int32),                      intent(in)    :: n_fe
      real(dp), dimension(n_forbs,n_fe),   intent(in)    :: mat_D
      real(dp), dimension(n_forbs,4*n_fe), intent(in)    :: DD2_o
      real(dp), dimension(4*n_fe),         intent(inout) :: DD2_ln_pff_P
      integer(int32) :: i1
      do i1 = 1_int32, n_fe
         call dgemv('T', n_forbs, 4, 1.0_dp, DD2_o(:,4*i1-3:4*i1), n_forbs, &
         & mat_D(:,i1), 1_int32, 0.0_dp, DD2_ln_pff_P(4*i1-3:4*i1), 1_int32)
         DD2_ln_pff_P(4*i1) = DD2_ln_pff_P(4*i1) - sum(DD2_ln_pff_P(4*i1-3:4*i1-1)**2)
      enddo
   end subroutine cmp_DD2_vec
   subroutine cmp_H_mat( n_fe, n_fe_u, n_fe_d, d_P, mat_o, inv_P, mat_OPi, mat_H )
      integer(int32),                      intent(in)    :: n_fe
      integer(int32),                      intent(in)    :: n_fe_u
      integer(int32),                      intent(in)    :: n_fe_d
      integer(int32),                      intent(in)    :: d_P
      real(dp), dimension(n_forbs,n_fe),   intent(in)    :: mat_o
      real(dp), dimension(d_P,d_P),        intent(in)    :: inv_P
      real(dp), dimension(n_forbs,n_fe_u), intent(inout) :: mat_OPi
      real(dp), dimension(n_forbs,*),      intent(inout) :: mat_H
      call dgemm('N', 'N', n_forbs, n_fe_u, n_fe_u, 1.0_dp, mat_o(:,1:n_fe_u), n_forbs, &
      & inv_P(1:n_fe_u,1:n_fe_u), n_fe_u, 0.0_dp, mat_OPi, n_forbs )
      call dgemm('N', 'T', n_forbs, n_forbs, n_fe_u, 0.5_dp, mat_OPi(:,1:n_fe_u), n_forbs, &
      & mat_o(:,1:n_fe_u), n_forbs, 0.0_dp, mat_H(:,1:n_forbs), n_forbs )
      if ( n_fe_d.gt.1_int32 ) then
         call dgemm('N', 'N', n_forbs, n_fe_d, n_fe_d, 1.0_dp, mat_o(:,n_fe_u+1:n_fe), n_forbs, &
         & inv_P(n_fe_u+1:n_fe,n_fe_u+1:n_fe), n_fe_d, 0.0_dp, mat_OPi(:,1:n_fe_d), n_forbs )
         if ( spin.eq.'U') then
            call dgemm('N', 'T', n_forbs, n_forbs, n_fe_d, 0.5_dp, mat_OPi(:,1:n_fe_d), n_forbs, &
            & mat_o(:,n_fe_u+1:n_fe), n_forbs, 0.0_dp, mat_H(:,n_forbs+1:2*n_forbs), n_forbs )
         else
            call dgemm('N', 'T', n_forbs, n_forbs, n_fe_d, 0.5_dp, mat_OPi(:,1:n_fe_d), n_forbs, &
            & mat_o(:,n_fe_u+1:n_fe), n_forbs, 1.0_dp, mat_H(:,1:n_forbs), n_forbs )
         endif
      endif
   end subroutine cmp_H_mat
   subroutine cmp_Lp_mat( n_fe, n_fe_u, n_fe_d, d_P, mat_o, inv_P, mat_OPi, mat_Lp, vec_lu, vec_ld )
      integer(int32),                      intent(in)    :: n_fe
      integer(int32),                      intent(in)    :: n_fe_u
      integer(int32),                      intent(in)    :: n_fe_d
      integer(int32),                      intent(in)    :: d_P
      real(dp), dimension(n_forbs,n_fe),   intent(in)    :: mat_o
      real(dp), dimension(d_P,d_P),        intent(in)    :: inv_P
      real(dp), dimension(n_forbs,n_fe_u), intent(inout) :: mat_OPi
      real(dp), dimension(n_forbs,n_forbs),intent(inout) :: mat_Lp
      real(dp), dimension(n_forbs),        intent(inout) :: vec_lu, vec_ld
      call dgemm('N', 'N', n_forbs, n_fe_u, n_fe_d, 1.0_dp, mat_o(:,n_fe_u+1:n_fe), n_forbs, &
      & inv_P(n_fe_u+1:n_fe,1:n_fe_u), n_fe_d, 0.0_dp, mat_OPi(:,1:n_fe_u), n_forbs )
      call dgemm('N', 'T', n_forbs, n_forbs, n_fe_u, 1.0_dp, mat_OPi(:,1:n_fe_u), n_forbs, &
      & mat_o(:,1:n_fe_u), n_forbs, 0.0_dp, mat_Lp(:,1:n_forbs), n_forbs )
      if ( d_P.gt.n_fe ) then
         call dgemv( 'N', n_forbs, n_fe_u, -1.0_dp, mat_o(1:n_forbs,1:n_fe_u), n_forbs, inv_P(1:n_fe_u,d_P), &
         & 1_int32, 0.0_dp, vec_lu(:), 1_int32)
         call dgemv( 'N', n_forbs, n_fe_d, -1.0_dp, mat_o(1:n_forbs,n_fe_u+1:n_fe), n_forbs, inv_P(n_fe_u+1:n_fe,d_P), &
         & 1_int32, 0.0_dp, vec_ld(:), 1_int32)
      endif
   end subroutine cmp_Lp_mat
end module pfaffian_operats_m
