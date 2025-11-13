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
module triplet_geminal_operats_m
   use fortran_kinds_v,              only: int32, dp
   use molecular_system_v,           only: n_at
   use fermionic_orbitals_m,         only: n_forbs, n_par_forbs, forbs
   use triplet_geminal_params_m, only: tgmcfs_t
   implicit none
   public :: cmp_Z_mat, upd_Z_mat, cmp_T_mat, cmp_T_vec, upd_T_inv, cmp_D_mat, &
   & upd_D_mat, cmp_H_mat, upd_H_mat, cmp_DD2_vec, upd_DD2_vec, cmp_da_vec
contains
   subroutine cmp_Z_mat( d_T, d_Z, n_fe, mat_Z_p, mat_o, mat_Z )
      integer(int32),                    intent(in)    :: d_T
      integer(int32),                    intent(in)    :: d_Z
      integer(int32),                    intent(in)    :: n_fe
      real(dp), dimension(n_forbs,d_Z),  intent(in)    :: mat_Z_p
      real(dp), dimension(n_forbs,n_fe), intent(in)    :: mat_o
      real(dp), dimension(n_forbs,d_T),  intent(inout) :: mat_Z
      call dgemm('T', 'N', n_forbs, n_fe, n_forbs, 1.0_dp, mat_Z_p(:,1:n_forbs), n_forbs, &
      & mat_o, n_forbs, 0.0_dp, mat_Z(:,1:n_fe), n_forbs )
      if( d_T.gt.n_fe ) mat_Z(:,d_T) = mat_Z_p(:,d_Z)
   end subroutine cmp_Z_mat
   subroutine upd_Z_mat( i_fe, n_fe, n_fe_u, n_fe_d, d_Tu, d_Td, tgm_c, mat_o, mat_Z )
      integer(int32),                         intent(in)    :: i_fe
      integer(int32),                         intent(in)    :: n_fe
      integer(int32),                         intent(in)    :: n_fe_u
      integer(int32),                         intent(in)    :: n_fe_d
      integer(int32),                         intent(in)    :: d_Tu
      integer(int32),                         intent(in)    :: d_Td
      type(tgmcfs_t),                         intent(in)    :: tgm_c
      real(dp), dimension(n_forbs,n_fe),      intent(in)    :: mat_o
      real(dp), dimension(n_forbs,d_Td+d_Tu), intent(inout) :: mat_Z
      if (i_fe.le.n_fe_u) then
         call dgemv( 'T', n_forbs, n_forbs, 1.0_dp, tgm_c%mat_Zu(:,1:n_forbs), n_forbs, &
         & mat_o(:,i_fe), 1_int32, 0.0_dp, mat_Z(:,i_fe), 1_int32 )
      else
         if(n_fe_d.gt.1_int32) then
            call dgemv( 'T', n_forbs, n_forbs, 1.0_dp, tgm_c%mat_Zd(:,1:n_forbs), n_forbs, &
            & mat_o(:,i_fe), 1_int32, 0.0_dp, mat_Z(:,d_Tu+i_fe-n_fe_u), 1_int32 )
         endif
      endif
   end subroutine upd_Z_mat
   subroutine cmp_T_mat( d_T, n_fe, mat_o, mat_Z, mat_T )
      integer(int32),                    intent(in)    :: d_T
      integer(int32),                    intent(in)    :: n_fe
      real(dp), dimension(n_forbs,n_fe), intent(in)    :: mat_o
      real(dp), dimension(n_forbs,d_T),  intent(in)    :: mat_Z
      real(dp), dimension(d_T,d_T),      intent(inout) :: mat_T
      integer(int32) :: i1
      call dgemm('T','N', d_T, n_fe, n_forbs, -1.0_dp, mat_Z, n_forbs, &
      & mat_o, n_forbs, 0.0_dp, mat_T(:,1:n_fe), d_T )
      if( d_T.gt.n_fe ) mat_T(:,d_T) = - mat_T(d_T,:)
      do i1 = 1_int32, d_T
         mat_T(i1,i1) = 0.0_dp
      enddo ! i1 (d_T)
   end subroutine cmp_T_mat
   subroutine cmp_T_vec( i_fe, d_T, vec_o, mat_Z, vec_T )
      integer(int32),                   intent(in)    :: i_fe
      integer(int32),                   intent(in)    :: d_T
      real(dp), dimension(n_forbs),     intent(in)    :: vec_o
      real(dp), dimension(n_forbs,d_T), intent(in)    :: mat_Z
      real(dp), dimension(d_T),         intent(inout) :: vec_T
      call dgemv( 'T', n_forbs, d_T, -1.0_dp, mat_Z(:,1:d_T), n_forbs,&
      & vec_o, 1_int32, 0.0_dp, vec_T(1:d_T), 1_int32 )
      vec_T(i_fe) = 0.0_dp
   end subroutine cmp_T_vec
   subroutine upd_T_inv( i_fe, d_T, g, vec_T, vec_Ti, vec_Tit, inv_T )
      integer(int32),               intent(in)    :: i_fe
      integer(int32),               intent(in)    :: d_T
      real(dp),                     intent(in)    :: g
      real(dp), dimension(d_T),     intent(in)    :: vec_T
      real(dp), dimension(d_T),     intent(in)    :: vec_Ti
      real(dp), dimension(d_T),     intent(inout) :: vec_Tit
      real(dp), dimension(d_T,d_T), intent(inout) :: inv_T
      call dgemv( 'T', d_T, d_T, -1.0_dp, inv_T, d_T, vec_T(1:d_T), 1_int32, &
      & 0.0_dp, vec_Tit(1:d_T), 1_int32)
      vec_Tit(i_fe) = vec_Tit(i_fe) - 1.0_dp
      call dger( d_T, d_T, -1.0_dp/g, vec_Tit(1:d_T), 1_int32, vec_Ti(1:d_T), 1_int32, &
      &  inv_T, d_T )
      call dger( d_T, d_T, +1.0_dp/g, vec_Ti(1:d_T), 1_int32, vec_Tit(1:d_T), 1_int32, &
      &  inv_T, d_T )
   end subroutine upd_T_inv
   subroutine cmp_D_mat( d_T, mat_Z, inv_T, mat_D )
      integer(int32),                   intent(in)  :: d_T
      real(dp), dimension(n_forbs,d_T), intent(in)  :: mat_Z
      real(dp), dimension(d_T,d_T),     intent(in)  :: inv_T
      real(dp), dimension(n_forbs,d_T), intent(out) :: mat_D
      call dgemm('N', 'N', n_forbs, d_T, d_T, 1.0_dp, mat_Z(:,1:d_T), n_forbs, &
      & inv_T, d_T, 0.0_dp, mat_D(:,1:d_T), n_forbs )
   end subroutine cmp_D_mat
   subroutine upd_D_mat ( i_fe, d_T, vec_Z, vec_T, vec_Tit, vec_Ti, g, &
   & vec_D, vec_Dt, mat_D )
      integer(int32),                   intent(in)  :: i_fe
      integer(int32),                   intent(in)  :: d_T
      real(dp), dimension(n_forbs),     intent(in)  :: vec_Z
      real(dp), dimension(d_T),         intent(in)  :: vec_T
      real(dp), dimension(d_T),         intent(in)  :: vec_Tit
      real(dp), dimension(d_T),         intent(in)  :: vec_Ti
      real(dp),                         intent(in)  :: g
      real(dp), dimension(n_forbs),     intent(inout) :: vec_D
      real(dp), dimension(n_forbs),     intent(inout) :: vec_Dt
      real(dp), dimension(n_forbs,d_T), intent(inout) :: mat_D
      call dgemv( 'N', n_forbs, d_T, 1.0_dp, mat_D, n_forbs, &
      & vec_T(1:d_T), 1_int32, 0.0_dp, vec_Dt, 1_int32)
      vec_Dt = vec_Dt - vec_Z
      vec_D  = mat_D(:,i_fe)
      call dger ( n_forbs, d_T, -1.0_dp/g, vec_Dt, 1_int32, vec_Ti(1:d_T), 1_int32, &
      & mat_D, n_forbs )
      call dger ( n_forbs, d_T, -1.0_dp/g, vec_D, 1_int32, vec_Tit(1:d_T), 1_int32, &
      & mat_D, n_forbs )
   end subroutine upd_D_mat
   subroutine cmp_H_mat( n_fe, d_T, d_Z, mat_o, inv_T, mat_OTi, mat_H )
      integer(int32),                    intent(in)    :: n_fe
      integer(int32),                    intent(in)    :: d_T
      integer(int32),                    intent(in)    :: d_Z
      real(dp), dimension(n_forbs,n_fe), intent(in)    :: mat_o
      real(dp), dimension(d_T,d_T),      intent(in)    :: inv_T
      real(dp), dimension(n_forbs,d_T),  intent(inout) :: mat_OTi
      real(dp), dimension(n_forbs,d_Z),  intent(inout) :: mat_H
      call dgemm('N','T', n_forbs, d_T, n_fe, -1.0_dp, mat_o, n_forbs, inv_T(:,1:n_fe), d_T, &
      & 0.0_dp, mat_OTi, n_forbs )
      call dgemm('N','T', n_forbs, n_forbs, n_fe, 0.5_dp, mat_OTi(:,1:n_fe), n_forbs, mat_o, n_forbs,&
      & 0.0_dp, mat_H(:,1:n_forbs), n_forbs )
      if ( d_Z .gt. n_forbs ) then
         mat_H(:,d_Z) = - mat_OTi(:,d_T)
      endif
   end subroutine cmp_H_mat
   subroutine upd_H_mat( n_fe, d_T, d_Z, mat_o, vec_o, vec_Tit, vec_Ti, g, vec_OTi, &
   & vec_OTit, mat_H )
      integer(int32),                    intent(in)  :: n_fe
      integer(int32),                    intent(in)  :: d_T
      integer(int32),                    intent(in)  :: d_Z
      real(dp), dimension(n_forbs,n_fe), intent(in)  :: mat_o
      real(dp), dimension(n_forbs),      intent(in)  :: vec_o
      real(dp), dimension(d_T),          intent(in)  :: vec_Tit
      real(dp), dimension(d_T),          intent(in)  :: vec_Ti
      real(dp),                          intent(in)  :: g
      real(dp), dimension(n_forbs),      intent(inout) :: vec_OTi
      real(dp), dimension(n_forbs),      intent(inout) :: vec_OTit
      real(dp), dimension(n_forbs,d_Z),  intent(inout) :: mat_H
      call dgemv( 'N', n_forbs, n_fe, 1.0_dp, mat_o, n_forbs, vec_Tit(1:n_fe), 1_int32,&
      & 0.0_dp, vec_OTit(1:n_forbs), 1_int32 )
      vec_OTit = vec_OTit - g * vec_o
      call dgemv( 'N', n_forbs, n_fe, 1.0_dp, mat_o, n_forbs, vec_Ti(1:n_fe), 1_int32,&
      & 0.0_dp, vec_OTi(1:n_forbs), 1_int32 )
      call dger( n_forbs, n_forbs,-0.5_dp / g, vec_OTi, 1_int32, vec_OTit,1_int32, &
      & mat_H(:,1:n_forbs), n_forbs )
      call dger( n_forbs, n_forbs, 0.5_dp / g, vec_OTit,1_int32, vec_OTi, 1_int32, &
      & mat_H(:,1:n_forbs), n_forbs )
      if ( d_T .gt. n_fe ) then
         mat_H(:,d_Z) = mat_H(:,d_Z) + vec_OTit(:) * vec_Ti(d_T) / g &
         & - vec_OTi(:) * vec_Tit(d_T) / g
      endif
   end subroutine upd_H_mat
   subroutine cmp_DD2_vec( d_Tu, n_fe, n_fe_u, n_fe_d, mat_D, DD2_o, DD2_ln_pf_T )
      integer(int32),                           intent(in)  :: d_Tu
      integer(int32),                           intent(in)  :: n_fe
      integer(int32),                           intent(in)  :: n_fe_u
      integer(int32),                           intent(in)  :: n_fe_d
      real(dp), dimension(n_forbs,d_Tu+n_fe_d), intent(in)  :: mat_D
      real(dp), dimension(n_forbs,4*n_fe),      intent(in)  :: DD2_o
      real(dp), dimension(4*n_fe),              intent(out) :: DD2_ln_pf_T
      integer(int32) :: i1
      do i1 = 1_int32, n_fe
         if(i1.le.n_fe_u) then
            call dgemv('T', n_forbs, 4, 1.0_dp, DD2_o(:,4*i1-3:4*i1), n_forbs, &
            & mat_D(:,i1), 1_int32, 0.0_dp, DD2_ln_pf_T(4*i1-3:4*i1), 1_int32)
         else
            call dgemv('T', n_forbs, 4, 1.0_dp, DD2_o(:,4*i1-3:4*i1), n_forbs, &
            & mat_D(:,d_Tu+i1-n_fe_u), 1_int32, 0.0_dp, DD2_ln_pf_T(4*i1-3:4*i1), 1_int32)
         endif
         DD2_ln_pf_T(4*i1) = DD2_ln_pf_T(4*i1) - sum(DD2_ln_pf_T(4*i1-3:4*i1-1)**2)
      enddo
   end subroutine cmp_DD2_vec
   subroutine upd_DD2_vec( n_fe, mat_D, DD2_o, DD2_ln_pf_T )
      integer(int32),                      intent(in)    :: n_fe
      real(dp), dimension(n_forbs,n_fe),   intent(in)    :: mat_D
      real(dp), dimension(n_forbs,4*n_fe), intent(in)    :: DD2_o
      real(dp), dimension(4*n_fe),         intent(inout) :: DD2_ln_pf_T
      integer(int32) :: i1
      do i1 = 1_int32, n_fe
         call dgemv('T', n_forbs, 4, 1.0_dp, DD2_o(:,4*i1-3:4*i1), n_forbs, &
         & mat_D(:,i1), 1_int32, 0.0_dp, DD2_ln_pf_T(4*i1-3:4*i1), 1_int32)
         DD2_ln_pf_T(4*i1) = DD2_ln_pf_T(4*i1) - sum(DD2_ln_pf_T(4*i1-3:4*i1-1)**2)
      enddo
   end subroutine upd_DD2_vec
   subroutine cmp_da_vec( n_fe, da_o, mat_D, da_ln_pf_T )
      integer(int32),                        intent(in)    :: n_fe
      real(dp), dimension(n_par_forbs,n_fe), intent(in)    :: da_o
      real(dp), dimension(n_forbs,n_fe),     intent(in)    :: mat_D
      real(dp), dimension(n_par_forbs),      intent(inout) :: da_ln_pf_T
      integer(int32) :: i1, i2, i3
      integer(int32) :: io, pi, pf
      da_ln_pf_T = 0.0_dp
      do i1 = 1_int32, n_fe
         do i2 = 1_int32, n_at ; do i3 = 1_int32, forbs(i2)%n_orbs
               io = forbs(i2)%orb(i3)%i_orb
               if( forbs(i2)%orb(i3)%n_par.gt.0_int32 ) then
                  pi = forbs(i2)%orb(i3)%i_par
                  pf = pi + forbs(i2)%orb(i3)%n_par - 1_int32
                  da_ln_pf_T(pi:pf) = da_ln_pf_T(pi:pf) + mat_D(io,i1) * da_o(pi:pf,i1)
               endif
            enddo ; enddo
      enddo
   end subroutine cmp_da_vec
end module triplet_geminal_operats_m
