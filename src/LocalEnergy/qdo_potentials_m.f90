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
module qdo_potentials_m
   use fortran_kinds_v,      only: dp, int32
   use physical_constants_v, only: pi
   use qdo_system_v,         only: n_qdo, qdos, d_qq, n_pchrg, d_pchpch, d_qpch, &
   & pchrgs, qdoham
   use drudonic_config_c,    only: drd_cnf
   use qdo_jastrow_m,        only: jdd_prs, jqd_prs
   implicit none
   public :: cmp_v_qq, cmp_v_dq, cmp_v_dd, upd_v_dq, cmp_v_pchpch, cmp_v_qpch, &
   & cmp_v_dpch, upd_v_dpch
contains
   subroutine cmp_v_qq( v_qq )
      real(dp), intent(inout) :: v_qq
      integer(int32) :: i1, i2, i3
      real(dp)       :: sigma
      v_qq = 0.0_dp
      if ( n_qdo.gt.1_int32 ) then
         i3 = 0_int32
         do i1 = 1_int32, n_qdo ;do i2 = i1 + 1_int32, n_qdo
            i3 = i3 + 1_int32
            if ( qdos(i1)%qdo_q*qdos(i2)%qdo_q.gt.10.0d-16 ) then
               sigma = sqrt( 2.0_dp*(qdos(i1)%qdo_sigma_c**2 + qdos(i2)%qdo_sigma_c**2 ) )
               if (sigma.le.1.0d-13) then
                  v_qq = v_qq + dble( qdos(i1)%qdo_q * qdos(i2)%qdo_q) / d_qq(i3)%m
               else
                  v_qq = v_qq + dble( qdos(i1)%qdo_q * qdos(i2)%qdo_q) * erf( d_qq(i3)%m / sigma ) / d_qq(i3)%m
               endif
            endif
         enddo ; enddo 
      endif 
   end subroutine cmp_v_qq
   subroutine cmp_v_dq( iw, v_dq )
      integer(int32),             intent(in)    :: iw
      real(dp), dimension(n_qdo), intent(inout) :: v_dq
      integer(int32) :: i1, i2
      real(dp)       :: sigma
      v_dq = 0.0_dp
      do i1 = 1_int32, n_qdo ; do i2 = 1_int32, n_qdo 
         if (i1.ne.i2) then
            if (qdoham.eq.'cou') then
               if (jqd_prs) then
                  v_dq(i1) = v_dq(i1) - qdos(i2)%qdo_q * qdos(i1)%qdo_q / drd_cnf(iw)%d_dq(i2,i1)%m
               else
                  sigma = sqrt(2.0_dp*(qdos(i1)%qdo_sigma_d**2 + qdos(i2)%qdo_sigma_c**2 )  )
                  if (sigma.gt.10d-10) then 
                     v_dq(i1) = v_dq(i1) - qdos(i2)%qdo_q * qdos(i1)%qdo_q * erf( drd_cnf(iw)%d_dq(i2,i1)%m / sigma ) / drd_cnf(iw)%d_dq(i2,i1)%m
                  else 
                     v_dq(i1) = v_dq(i1) - qdos(i2)%qdo_q * qdos(i1)%qdo_q / drd_cnf(iw)%d_dq(i2,i1)%m
                  endif
               endif
            endif
         else
            v_dq(i1) = v_dq(i1) + &
            & 0.5_dp * qdos(i1)%qdo_m * ( qdos(i1)%qdo_omega**2 ) * ( drd_cnf(iw)%d_dq(i1,i1)%m**2 )
         endif
      enddo ; enddo 
   end subroutine cmp_v_dq
   subroutine cmp_v_dd( iw, v_dd )
      integer(int32),                         intent(in)    :: iw
      real(dp), dimension(n_qdo*(n_qdo-1)/2), intent(inout) :: v_dd
      integer(int32) :: i1, i2, i3
      real(dp)       :: sigma
      if (qdoham.eq.'cou') then
         i3 = 0_int32
         do i1 = 1_int32, n_qdo ; do i2 = i1 + 1_int32, n_qdo
            i3 = i3 + 1_int32
            if (jdd_prs) then
               v_dd(i3) = qdos(i1)%qdo_q * qdos(i2)%qdo_q / drd_cnf(iw)%d_dd(i3)%m
            else
               sigma = sqrt( 2.0_dp *( qdos(i1)%qdo_sigma_d**2 + qdos(i1)%qdo_sigma_d**2 )  )
               if (sigma.gt.10d-10) then 
                  v_dd(i3) = qdos(i1)%qdo_q * qdos(i2)%qdo_q * erf( drd_cnf(iw)%d_dd(i3)%m / sigma ) / drd_cnf(iw)%d_dd(i3)%m
               else 
                  v_dd(i3) = qdos(i1)%qdo_q * qdos(i2)%qdo_q / drd_cnf(iw)%d_dd(i3)%m
               endif
            endif
         enddo ; enddo 
      else
         i3 = 0_int32
         do i1 = 1_int32, n_qdo ; do i2 = i1 + 1_int32, n_qdo
            i3 = i3 + 1_int32
            v_dd(i3) = qdos(i1)%qdo_q * qdos(i2)%qdo_q / d_qq(i3)%m**3 * &
            & ( dot_product(drd_cnf(iw)%d_dq(i1,i1)%v(1:3), drd_cnf(iw)%d_dq(i2,i2)%v(1:3)) - &
            &  3.0_dp * dot_product(drd_cnf(iw)%d_dq(i1,i1)%v(1:3), d_qq(i3)%v(1:3)) * dot_product(drd_cnf(iw)%d_dq(i2,i2)%v(1:3),d_qq(i3)%v(1:3)) / d_qq(i3)%m**2 )
         enddo ; enddo 
      endif
   end subroutine cmp_v_dd
   subroutine upd_v_dq( iw, v_dq )
      integer(int32),             intent(in)    :: iw
      real(dp), dimension(n_qdo), intent(inout) :: v_dq
      integer(int32) :: i1
      real(dp)       :: sigma
      v_dq(drd_cnf(iw)%i_drd) = 0.0_dp
      do i1 = 1_int32, n_qdo 
         if (drd_cnf(iw)%i_drd.ne.i1) then
            if (qdoham.eq.'cou') then
               if (jqd_prs) then
                  v_dq(drd_cnf(iw)%i_drd) = v_dq(drd_cnf(iw)%i_drd) - qdos(i1)%qdo_q * qdos(drd_cnf(iw)%i_drd)%qdo_q / drd_cnf(iw)%d_dq(i1,drd_cnf(iw)%i_drd)%m
               else
                  sigma = sqrt( 2.0_dp * (qdos(i1)%qdo_sigma_c**2 + qdos( drd_cnf(iw)%i_drd )%qdo_sigma_d**2 )  )
                  if (sigma.gt.10d-10) then 
                     v_dq(drd_cnf(iw)%i_drd) = v_dq(drd_cnf(iw)%i_drd) - qdos(i1)%qdo_q * qdos(drd_cnf(iw)%i_drd)%qdo_q * &
                     & erf( drd_cnf(iw)%d_dq(i1,drd_cnf(iw)%i_drd)%m / sigma ) / drd_cnf(iw)%d_dq(i1,drd_cnf(iw)%i_drd)%m
                  else 
                     v_dq(drd_cnf(iw)%i_drd) = v_dq(drd_cnf(iw)%i_drd) - qdos(i1)%qdo_q * qdos(drd_cnf(iw)%i_drd)%qdo_q &
                     & / drd_cnf(iw)%d_dq(i1,drd_cnf(iw)%i_drd)%m
                  endif
               endif
            endif
         else
            v_dq(drd_cnf(iw)%i_drd) = v_dq(drd_cnf(iw)%i_drd) + &
            & 0.5_dp * qdos(i1)%qdo_m * ( qdos(i1)%qdo_omega**2 ) * ( drd_cnf(iw)%d_dq(i1,i1)%m**2 )
         endif
      enddo 
   end subroutine upd_v_dq
   subroutine cmp_v_pchpch( v_pchpch )
      real(dp), intent(inout) :: v_pchpch
      integer(int32) :: i1, i2, i3
      real(dp)       :: sigma
      v_pchpch = 0.0_dp
      if ( n_pchrg.gt.1_int32 ) then
         i3 = 0_int32
         do i1 = 1_int32, n_pchrg ;do i2 = i1 + 1_int32, n_pchrg
            i3 = i3 + 1_int32
            sigma = sqrt( pchrgs(i1)%pchrg_sigma**2 + pchrgs(i2)%pchrg_sigma**2 )*sqrt( 2.0_dp ) 
            if (sigma.gt.10d-10) then
               v_pchpch = v_pchpch + dble( pchrgs(i1)%pchrg_q * pchrgs(i2)%pchrg_q ) * erf( d_pchpch(i3)%m   / sigma ) / d_pchpch(i3)%m
            else                
               v_pchpch = v_pchpch + dble( pchrgs(i1)%pchrg_q * pchrgs(i2)%pchrg_q ) / d_pchpch(i3)%m
            endif
         enddo ; enddo 
      endif 
   end subroutine cmp_v_pchpch
   subroutine cmp_v_qpch( v_qpch )
      real(dp), intent(inout) :: v_qpch
      real(dp)       :: sigma
      integer(int32) :: i1, i2
      v_qpch = 0.0_dp
      if (n_qdo.gt.0_int32.and.n_pchrg.gt.0_int32) then
         do i1 = 1, n_qdo ; do i2 = 1, n_pchrg
            sigma = sqrt( qdos(i1)%qdo_sigma_c**2 + pchrgs(i2)%pchrg_sigma**2 )*sqrt( 2.0_dp )   
            if (sigma.gt.10d-10) then
               if (d_qpch(i1,i2)%m.gt.1e-6) then
                  v_qpch = v_qpch + qdos(i1)%qdo_q * pchrgs(i2)%pchrg_q * erf( d_qpch(i1,i2)%m   / sigma ) / d_qpch(i1,i2)%m
               else
                  v_qpch = v_qpch + qdos(i1)%qdo_q * pchrgs(i2)%pchrg_q * 2.0_dp / ( sigma * sqrt( pi ) )
               endif
            else 
               v_qpch = v_qpch + qdos(i1)%qdo_q * pchrgs(i2)%pchrg_q / d_qpch(i1,i2)%m
            endif
         enddo ; enddo 
      endif
   end subroutine  cmp_v_qpch
   subroutine cmp_v_dpch( iw, v_dpch )
      integer(int32),             intent(in)     :: iw
      real(dp), dimension(n_qdo), intent(inout)  :: v_dpch
      real(dp)       :: sigma
      integer(int32) :: i1, i2
      v_dpch(:) = 0.0_dp
      do i1 = 1_int32, n_qdo ; do i2 = 1_int32, n_pchrg
         if ( qdos(i1)%qdo_name.ne.pchrgs(i2)%pchrg_name ) then
            sigma = sqrt( qdos(i1)%qdo_sigma_d**2 + pchrgs(i2)%pchrg_sigma**2 )*sqrt( 2.0_dp )
            if (sigma.gt.10d-10) then
               v_dpch(i1) = v_dpch(i1) - qdos(i1)%qdo_q * pchrgs(i2)%pchrg_q *&
               & erf( drd_cnf(iw)%d_dpch(i1,i2)%m / sigma ) / drd_cnf(iw)%d_dpch(i1,i2)%m
            else 
               v_dpch(i1) = v_dpch(i1) - qdos(i1)%qdo_q * pchrgs(i2)%pchrg_q / drd_cnf(iw)%d_dpch(i1,i2)%m
            endif
         endif
      enddo ; enddo 
   end subroutine cmp_v_dpch
   subroutine upd_v_dpch( iw, v_dpch )
      integer(int32),             intent(in)    :: iw
      real(dp), dimension(n_qdo), intent(inout) :: v_dpch
      real(dp)       :: sigma
      integer(int32) :: i1
      v_dpch(drd_cnf(iw)%i_drd) = 0.0_dp
      do i1 = 1_int32, n_pchrg
         if ( qdos(drd_cnf(iw)%i_drd)%qdo_name.ne.pchrgs(i1)%pchrg_name ) then
            sigma = sqrt( qdos(drd_cnf(iw)%i_drd)%qdo_sigma_d**2 + pchrgs(i1)%pchrg_sigma**2 )*sqrt( 2.0_dp )
            if (sigma.gt.10d-10) then
               v_dpch(drd_cnf(iw)%i_drd) = v_dpch(drd_cnf(iw)%i_drd) - qdos(drd_cnf(iw)%i_drd)%qdo_q * pchrgs(i1)%pchrg_q *&
               & erf( drd_cnf(iw)%d_dpch(drd_cnf(iw)%i_drd,i1)%m  / sigma ) / drd_cnf(iw)%d_dpch(drd_cnf(iw)%i_drd,i1)%m
            else 
               v_dpch(drd_cnf(iw)%i_drd) = v_dpch(drd_cnf(iw)%i_drd) - qdos(drd_cnf(iw)%i_drd)%qdo_q * pchrgs(i1)%pchrg_q &
               & / drd_cnf(iw)%d_dpch(drd_cnf(iw)%i_drd,i1)%m
            endif
         endif
      enddo
   end subroutine upd_v_dpch
end module qdo_potentials_m
