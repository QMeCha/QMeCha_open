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
module embedding_potential_m
   use fortran_kinds_v,    only: dp, int32
   use qdo_system_v,       only: n_qdo, qdos, d_qq, d_qn, el_sigma, n_pchrg, &
   & d_npch, pchrgs
   use molecular_system_v, only: n_at, atoms, n_fe, f_chrg
   use pseudopotentials_m, only: eval_local_cmp
   use pseudopotentials_v, only: psdpot
   use drudonic_config_c,  only: drd_cnf
   use fermionic_config_c, only: frm_cnf
   implicit none
   public :: cmp_v_qn, cmp_v_dn, cmp_v_fqd, upd_v_dn, cmp_v_npch, cmp_v_fpch, &
   & upd_v_fpch
contains
   subroutine cmp_v_qn( v_qn )
      real(dp), intent(inout) :: v_qn
      integer(int32) :: i1, i2
      v_qn = 0.0_dp
      if (n_at.gt.0_int32.and.n_qdo.gt.0_int32) then
         do i1 = 1, n_qdo ; do i2 = 1, n_at
               if ( abs(qdos(i1)%qdo_q * dble(atoms(i2)%atm_z)).gt.10d-15) then
                  v_qn = v_qn + qdos(i1)%qdo_q * dble(atoms(i2)%atm_z) / d_qn(i1,i2)%m
               endif
            enddo ; enddo
      endif
   end subroutine cmp_v_qn
   subroutine cmp_v_dn( iw, v_dn )
      integer(int32),             intent(in)    :: iw
      real(dp), dimension(n_qdo), intent(inout) :: v_dn
      integer(int32) :: i1, i2
      v_dn(:) = 0.0_dp
      do i1 = 1_int32, n_qdo ; do i2 = 1_int32, n_at
         v_dn(i1) = v_dn(i1) &
            & - qdos(i1)%qdo_q * dble( atoms(i2)%atm_z) / drd_cnf(iw)%d_dn(i1,i2)%m
      enddo ; enddo 
   end subroutine cmp_v_dn
   subroutine cmp_v_fqd( iw, v_fq, v_df )
      integer(int32),                  intent(in)    :: iw
      real(dp), dimension(n_fe),       intent(inout) :: v_fq
      real(dp), dimension(n_qdo,n_fe), intent(inout) :: v_df
      integer(int32) :: i1, i2
      real(dp)       :: sigma
      v_fq(:) = 0.0_dp
      do i1 = 1_int32, n_fe ; do i2 = 1_int32, n_qdo
         sigma = sqrt( el_sigma**2 + qdos(i2)%qdo_sigma_c**2 )*sqrt( 2.0_dp )
         if (sigma.gt.10d-16) then
            v_fq(i1) = v_fq(i1) + dble(f_chrg(i1))*qdos(i2)%qdo_q*erf( drd_cnf(iw)%d_fq(i1,i2)%m / sigma ) / drd_cnf(iw)%d_fq(i1,i2)%m
         else 
            v_fq(i1) = v_fq(i1) + dble(f_chrg(i1))*qdos(i2)%qdo_q/ drd_cnf(iw)%d_fq(i1,i2)%m
         endif
      enddo ; enddo
      do i1 = 1, n_qdo ; do i2 = 1, n_fe
         sigma = sqrt( qdos(i1)%qdo_sigma_d**2 + el_sigma**2 )*sqrt( 2.0_dp )
         v_df(i1,i2) = - dble(f_chrg(i2))*qdos(i1)%qdo_q*erf( drd_cnf(iw)%d_df(i1,i2)%m  / sigma ) / drd_cnf(iw)%d_df(i1,i2)%m
      enddo ; enddo 
   end subroutine cmp_v_fqd
   subroutine upd_v_dn( iw, v_dn )
      integer(int32),               intent(in)    :: iw
      real(dp), dimension(n_qdo), intent(inout) :: v_dn
      integer(int32) :: i1
      v_dn(drd_cnf(iw)%i_drd) = 0.0_dp
      do i1 = 1_int32, n_at
         v_dn(drd_cnf(iw)%i_drd) = v_dn(drd_cnf(iw)%i_drd) &
            - qdos(drd_cnf(iw)%i_drd)%qdo_q * dble( atoms(i1)%atm_z) / drd_cnf(iw)%d_dn(drd_cnf(iw)%i_drd,i1)%m
      enddo 
   end subroutine upd_v_dn
   subroutine upd_v_fqd( iw, v_fq, v_df, i_typ )
      integer(int32),                      intent(in)    :: iw
      real(dp), dimension(n_fe),         intent(inout) :: v_fq
      real(dp), dimension(n_qdo,n_fe), intent(inout) :: v_df
      integer(int32),                      intent(in)    :: i_typ
      integer(int32) :: i1
      real(dp)       :: sigma
      if (i_typ.eq.0_int32) then
         v_fq(frm_cnf(iw)%i_fe)  = 0.0_dp
         do i1 = 1_int32, n_qdo
            sigma = sqrt( el_sigma**2 + qdos(i1)%qdo_sigma_c**2 )*sqrt( 2.0_dp )
            if (sigma.gt.10d-16) then
               v_fq(frm_cnf(iw)%i_fe) = v_fq(frm_cnf(iw)%i_fe) + f_chrg(frm_cnf(iw)%i_fe)*qdos(i1)%qdo_q &
               & * erf( drd_cnf(iw)%d_fq(frm_cnf(iw)%i_fe,i1)%m  / sigma ) / drd_cnf(iw)%d_fq(frm_cnf(iw)%i_fe,i1)%m
            else 
               v_fq(frm_cnf(iw)%i_fe) = v_fq(frm_cnf(iw)%i_fe) + f_chrg(frm_cnf(iw)%i_fe)*qdos(i1)%qdo_q &
               & / drd_cnf(iw)%d_fq(frm_cnf(iw)%i_fe,i1)%m
            endif
         enddo
         do i1 = 1, n_qdo
            sigma = sqrt( qdos(i1)%qdo_sigma_d**2 + el_sigma**2 )*sqrt( 2.0_dp )
            if (sigma.gt.10d-16) then
               v_df(i1,frm_cnf(iw)%i_fe) = - f_chrg(frm_cnf(iw)%i_fe)*qdos(i1)%qdo_q*&
               &erf( drd_cnf(iw)%d_df(i1,frm_cnf(iw)%i_fe)%m  /sigma ) / drd_cnf(iw)%d_df(i1,frm_cnf(iw)%i_fe)%m
            else 
               v_df(i1,frm_cnf(iw)%i_fe) = - f_chrg(frm_cnf(iw)%i_fe)*qdos(i1)%qdo_q&
               &/ drd_cnf(iw)%d_df(i1,frm_cnf(iw)%i_fe)%m
            endif
         enddo
      else 
         do i1 = 1, n_fe
            sigma = sqrt( qdos(drd_cnf(iw)%i_drd)%qdo_sigma_d**2 + el_sigma**2 )*sqrt( 2.0_dp )
            if (sigma.gt.10d-16) then
               v_df(drd_cnf(iw)%i_drd,i1) = - f_chrg(i1)*qdos(drd_cnf(iw)%i_drd)%qdo_q&
               &*erf( drd_cnf(iw)%d_df(drd_cnf(iw)%i_drd,i1)%m /sigma ) / drd_cnf(iw)%d_df(drd_cnf(iw)%i_drd,i1)%m
            else 
               v_df(drd_cnf(iw)%i_drd,i1) = - f_chrg(i1)*qdos(drd_cnf(iw)%i_drd)%qdo_q&
               & / drd_cnf(iw)%d_df(drd_cnf(iw)%i_drd,i1)%m
            endif
         enddo
      endif
   end subroutine upd_v_fqd
   subroutine cmp_v_npch( v_npch )
      real(dp), intent(inout) :: v_npch
      integer(int32) :: i1, i2
      v_npch = 0.0_dp
      if (n_at.gt.0_int32.and.n_pchrg.gt.0_int32) then
         do i1 = 1, n_at ; do i2 = 1, n_pchrg
            v_npch = v_npch + dble(atoms(i1)%atm_z) * pchrgs(i2)%pchrg_q / d_npch(i1,i2)%m
         enddo ; enddo 
      endif
   end subroutine cmp_v_npch
   subroutine cmp_v_fpch( iw, v_fpch )
      integer(int32),              intent(in)     :: iw
      real(dp), dimension(n_fe), intent(inout) :: v_fpch
      integer(int32) :: i1, i2
      real(dp)       :: sigma
      v_fpch(:) = 0.0_dp
      if (n_pchrg.gt.0_int32) then
         do i1 = 1_int32, n_fe
            do i2 = 1_int32, n_pchrg
               sigma = sqrt( el_sigma**2 + pchrgs(i2)%pchrg_sigma**2 )*sqrt( 2.0_dp )
               if (sigma.gt.10d-16) then
                  v_fpch(i1) = v_fpch(i1) + dble(f_chrg(i1)) * pchrgs(i2)%pchrg_q *&
                  & erf( drd_cnf(iw)%d_fpch(i1,i2)%m / sigma ) / drd_cnf(iw)%d_fpch(i1,i2)%m
               else 
                  v_fpch(i1) = v_fpch(i1) + dble(f_chrg(i1)) * pchrgs(i2)%pchrg_q &
                  & / drd_cnf(iw)%d_fpch(i1,i2)%m               
               endif
            enddo 
         enddo
      endif
   end subroutine cmp_v_fpch
   subroutine upd_v_fpch( iw, v_fpch )
      integer(int32),             intent(in)     :: iw
      real(dp), dimension(n_fe), intent(inout) :: v_fpch
      integer(int32) :: i1
      real(dp)       :: sigma
      v_fpch(frm_cnf(iw)%i_fe)  = 0.0_dp
      if (n_pchrg.gt.0_int32) then
         do i1 = 1_int32, n_pchrg
            sigma = sqrt( el_sigma**2 + pchrgs(i1)%pchrg_sigma**2 )*sqrt( 2.0_dp )
            if (sigma.gt.10d-16) then
               v_fpch(frm_cnf(iw)%i_fe) = v_fpch(frm_cnf(iw)%i_fe) + f_chrg(frm_cnf(iw)%i_fe) * pchrgs(i1)%pchrg_q *&
               & erf( drd_cnf(iw)%d_fpch(frm_cnf(iw)%i_fe,i1)%m  / sigma ) / drd_cnf(iw)%d_fpch(frm_cnf(iw)%i_fe,i1)%m
            else 
               v_fpch(frm_cnf(iw)%i_fe) = v_fpch(frm_cnf(iw)%i_fe) + f_chrg(frm_cnf(iw)%i_fe) * pchrgs(i1)%pchrg_q &
               & / drd_cnf(iw)%d_fpch(frm_cnf(iw)%i_fe,i1)%m
            endif
         enddo
      endif
   end subroutine upd_v_fpch
end module embedding_potential_m
