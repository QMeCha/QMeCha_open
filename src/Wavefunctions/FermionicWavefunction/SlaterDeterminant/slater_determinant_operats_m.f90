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
module slater_determinant_operats_m
   use fortran_kinds_v,      only: int32, dp
   use molecular_system_v,   only: n_at
   use fermionic_orbitals_m, only: n_forbs, forbs, n_par_forbs
   implicit none
   public :: cmp_S_mat, cmp_S_vec, upd_S_inv_det, cmp_D_mat, upd_D_mat, &
   & cmp_K_mat, cmp_da_sld
contains
   subroutine cmp_S_mat( n_orbs, n_e, mat_o, mat_L, mat_S )
      integer(int32),                  intent(in)  :: n_orbs
      integer(int32),                  intent(in)  :: n_e
      real(dp), dimension(n_orbs,n_e), intent(in)  :: mat_o
      real(dp), dimension(n_orbs,n_e), intent(in)  :: mat_L
      real(dp), dimension(n_e,n_e),    intent(out) :: mat_S
      if(n_e.eq.1_int32) then
         mat_S(1,1) = dot_product( mat_L(:,1), mat_o(:,1) )
      else
         call dgemm('T', 'N', n_e, n_e, n_orbs, 1.0_dp, mat_L, n_orbs, &
         & mat_o, n_orbs, 0.0_dp, mat_S, n_e)
      endif
   end subroutine cmp_S_mat
   subroutine cmp_S_vec( n_orbs, n_e, o_new, mat_L, vec_S_new )
      integer(int32),                  intent(in)    :: n_orbs
      integer(int32),                  intent(in)    :: n_e
      real(dp), dimension(n_orbs),     intent(in)    :: o_new
      real(dp), dimension(n_orbs,n_e), intent(in)    :: mat_L
      real(dp), dimension(n_e),        intent(inout) :: vec_S_new
      if( n_e.eq.1_int32 ) then
         vec_S_new(1) = vec_S_new(1) + dot_product( mat_L(:,1), o_new )
      else
         call dgemv('T', n_orbs, n_e, 1.0_dp, mat_L(:,1:n_e), n_orbs, &
         & o_new, 1_int32, 1.0_dp, vec_S_new(1:n_e), 1_int32)
      endif
   end subroutine cmp_S_vec
   subroutine upd_S_inv_det( i_fe, n_e, n_mat_chngs, mat_S, inv_S, det_S, &
   & vec_S_new, vec_v_s, vec_v, g )
      integer(int32),               intent(in)    :: i_fe
      integer(int32),               intent(in)    :: n_e
      integer(int32),               intent(inout) :: n_mat_chngs
      real(dp), dimension(n_e,n_e), intent(inout) :: mat_S
      real(dp), dimension(n_e,n_e), intent(inout) :: inv_S
      real(dp),                     intent(inout) :: det_S
      real(dp), dimension(n_e),     intent(inout) :: vec_S_new
      real(dp), dimension(n_e),     intent(inout) :: vec_v_s
      real(dp), dimension(n_e),     intent(inout) :: vec_v
      real(dp),                     intent(inout) :: g
      det_S = det_S * g
      if( n_e.eq.1_int32 ) then
         mat_S(1,1) = vec_S_new(1)
         inv_S(1,1) = 1.0_dp / mat_S(1,1)
      else
         mat_S(1:n_e,i_fe) = vec_S_new(1:n_e)
         if( n_mat_chngs.eq.n_e .or. abs(g).lt.10d-14 ) then
            call lu_dec_det_inv( n_e, mat_S, det_S, inv_S )
            n_mat_chngs = 0_int32
         else
            call sm_upd_inv('C', n_e, i_fe, vec_S_new, vec_v_s,&
            & vec_v, g, inv_S)
            n_mat_chngs = n_mat_chngs + 1_int32
         endif
      endif
   end subroutine upd_S_inv_det
   subroutine cmp_D_mat( n_e, mat_L, inv_S, mat_D )
      integer(int32),                   intent(in)    :: n_e
      real(dp), dimension(n_forbs,n_e), intent(in)    :: mat_L
      real(dp), dimension(n_e,n_e),     intent(inout) :: inv_S
      real(dp), dimension(n_forbs,n_e), intent(inout) :: mat_D
      if( n_e.eq.1_int32 ) then
         mat_D(:,1) = inv_S(1,1) * mat_L(:,1)
      else
         call dgemm('N', 'T', n_forbs, n_e, n_e, 1.0_dp, mat_L, n_forbs,&
         & inv_S, n_e, 0.0_dp, mat_D(:,:), n_forbs )
      endif
   end subroutine cmp_D_mat
   subroutine upd_D_mat( n_e, i_fe, n_mat_chngs, mat_L, inv_S, mat_D, vec_v_d, vec_v, g)
      integer(int32),                   intent(in)    :: n_e
      integer(int32),                   intent(in)    :: i_fe
      integer(int32),                   intent(inout) :: n_mat_chngs
      real(dp), dimension(n_forbs,n_e), intent(in)    :: mat_L
      real(dp), dimension(n_e,n_e),     intent(inout) :: inv_S
      real(dp), dimension(n_forbs,n_e), intent(inout) :: mat_D
      real(dp), dimension(n_forbs),     intent(inout) :: vec_v_d
      real(dp), dimension(n_e),         intent(in)    :: vec_v
      real(dp),                         intent(in)    :: g
      if(n_e.eq.1_int32) then
         mat_D(:,1) = inv_S(1,1) * mat_L(:,1)
      elseif(n_mat_chngs.eq.0_int32) then
         call dgemm('N', 'T', n_forbs, n_e, n_e, 1.0_dp, mat_L, n_forbs,&
         & inv_S, n_e, 0.0_dp, mat_D(:,:), n_forbs )
      else
         vec_v_d(:) = mat_D(:,i_fe)
         call dger ( n_forbs, n_e, -1.0_dp/g, vec_v_d, 1_int32,&
         & vec_v, 1_int32, mat_D(:,:), n_forbs )
      endif
   end subroutine upd_D_mat
   subroutine cmp_K_mat( n_e, mat_o, inv_S, mat_K )
      integer(int32),                   intent(in)    :: n_e
      real(dp), dimension(n_forbs,n_e), intent(in)    :: mat_o
      real(dp), dimension(n_e,n_e),     intent(in)    :: inv_S
      real(dp), dimension(n_forbs,n_e), intent(inout) :: mat_K
      if( n_e .eq. 1_int32 ) then
         mat_K(:,1) = mat_o(:,1) * inv_S(1,1)
      else
         call dgemm('N', 'N', n_forbs, n_e, n_e, 1.0_dp, mat_o(:,1:n_e), n_forbs, &
         & inv_S, n_e, 0.0_dp, mat_K, n_forbs )
      endif
   end subroutine cmp_K_mat
   subroutine cmp_da_sld( n_e, da_o, mat_D, da_ln_det_S )
      integer(int32),                       intent(in)  :: n_e
      real(dp), dimension(n_par_forbs,n_e), intent(in)  :: da_o
      real(dp), dimension(n_forbs,n_e),     intent(in)  :: mat_D
      real(dp), dimension(n_par_forbs),     intent(out) :: da_ln_det_S
      integer(int32) :: i1, i2, i3
      integer(int32) :: io, pi, pf
      da_ln_det_S = 0.0_dp
      do i1 = 1_int32, n_e
         do i2 = 1_int32, n_at ; do i3 = 1_int32, forbs(i2)%n_orbs
               io = forbs(i2)%orb(i3)%i_orb
               if( forbs(i2)%orb(i3)%n_par.gt.0_int32 ) then
                  pi = forbs(i2)%orb(i3)%i_par
                  pf = pi + forbs(i2)%orb(i3)%n_par - 1_int32
                  da_ln_det_S(pi:pf) = da_ln_det_S(pi:pf) + mat_D(io,i1) * da_o(pi:pf,i1)
               endif
            enddo ; enddo  !  i2 (n_at) i3 (orbs(i2)%n_orbs)
      enddo ! i1 (n_el_u)
   end subroutine cmp_da_sld
end module slater_determinant_operats_m
