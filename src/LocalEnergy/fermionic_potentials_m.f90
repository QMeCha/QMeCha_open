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
module fermionic_potentials_m
   use fortran_kinds_v,    only: dp, int32
   use molecular_system_v, only: n_at, n_fe, n_el, n_po, atoms, d_nn, f_chrg
   use fermionic_config_c, only: frm_cnf
   implicit none
   public :: cmp_v_nn, cmp_v_fn, cmp_v_ee, cmp_v_pp, cmp_v_ep, upd_v_fn, upd_v_ee,&
   & upd_v_pp, upd_v_ep
contains
   subroutine cmp_v_nn( v_nn )
      real(dp), intent(inout) :: v_nn
      integer(int32) :: i1, i2, i3
      v_nn = 0.0_dp
      if( n_at.gt.1_int32 ) then
         i3 = 0_int32
         do i1 = 1_int32, n_at ; do i2 = i1 + 1_int32, n_at
            i3 = i3 + 1_int32
            if (d_nn(i3)%m.gt.0.1_dp) v_nn = v_nn + dble( atoms(i1)%atm_z * atoms(i2)%atm_z ) / d_nn(i3)%m
         enddo ; enddo 
      endif 
   end subroutine cmp_v_nn
   subroutine cmp_v_fn( iw, v_fn )
      integer(int32),              intent(in)    :: iw
      real(dp), dimension(n_fe), intent(inout) :: v_fn
      integer(int32) :: i1, i2
      v_fn = 0.0_dp
      do i1 = 1_int32, n_fe
         if (atoms(1)%atm_z.ne.0_int32) &
         & v_fn(i1) = dble(atoms(1)%atm_z * f_chrg(i1)) / frm_cnf(iw)%d_fn(1,i1)%m
         if(n_at.gt.1_int32) then
            do i2 = 2_int32, n_at
               if (atoms(i2)%atm_z.ne.0_int32) &
               & v_fn(i1) = v_fn(i1) + dble(atoms(i2)%atm_z * f_chrg(i1)) / frm_cnf(iw)%d_fn(i2,i1)%m
            enddo 
         endif
      enddo 
   end subroutine cmp_v_fn
   subroutine cmp_v_ee( iw, v_ee )
      integer(int32),                         intent(in)    :: iw
      real(dp), dimension(n_el*(n_el-1)/2), intent(inout) :: v_ee
      integer(int32) :: i1, i2, i3
      v_ee = 0.0_dp
      i3 = 0_int32
      do i1 = 1_int32, n_el ; do i2 = i1 + 1_int32, n_el
         i3 = i3 + 1_int32
         v_ee(i3) = 1.0_dp / frm_cnf(iw)%d_ee(i3)%m
      enddo ; enddo 
   end subroutine cmp_v_ee
   subroutine cmp_v_pp( iw, v_pp )
      integer(int32),                       intent(in)    :: iw
      real(dp), dimension(n_po*(n_po-1)/2), intent(inout) :: v_pp
      integer(int32) :: i1, i2, i3
      v_pp = 0.0_dp
      i3 = 0_int32
      do i1 = 1_int32, n_po ; do i2 = i1 + 1_int32, n_po
         i3 = i3 + 1_int32
         v_pp(i3) = 1.0_dp / frm_cnf(iw)%d_pp(i3)%m
      enddo ; enddo 
   end subroutine cmp_v_pp
   subroutine cmp_v_ep( iw, v_ep )
      integer(int32),                 intent(in)    :: iw
      real(dp), dimension(n_po,n_el), intent(inout) :: v_ep
      integer(int32) :: i1, i2
      v_ep = 0.0_dp
      do i1 = 1_int32, n_el ; do i2 = 1_int32, n_po
         v_ep(i2,i1) = - 1.0_dp / frm_cnf(iw)%d_ep(i2,i1)%m
      enddo ; enddo 
   end subroutine cmp_v_ep
   subroutine upd_v_fn( iw, v_fn )
      integer(int32),              intent(in)    :: iw
      real(dp), dimension(n_fe), intent(inout) :: v_fn
      integer(int32) :: i1
      if (atoms(1)%atm_z.ne.0_int32) &
      & v_fn(frm_cnf(iw)%i_fe) = dble(atoms(1)%atm_z * f_chrg(frm_cnf(iw)%i_fe)) &
      & / frm_cnf(iw)%d_fn(1,frm_cnf(iw)%i_fe)%m
      if(n_at.ge.2_int32) then
         do i1 = 2_int32, n_at
            if (atoms(i1)%atm_z.ne.0_int32) &
            & v_fn(frm_cnf(iw)%i_fe) = v_fn(frm_cnf(iw)%i_fe) + dble(atoms(i1)%atm_z) *&
            & dble(f_chrg(frm_cnf(iw)%i_fe)) / frm_cnf(iw)%d_fn(i1,frm_cnf(iw)%i_fe)%m
         enddo 
      endif
   end subroutine upd_v_fn
   subroutine upd_v_ee( iw, v_ee )
      integer(int32),                       intent(in)    :: iw
      real(dp), dimension(n_el*(n_el-1)/2), intent(inout) :: v_ee
      integer(int32) :: i1, i2
      if ( frm_cnf(iw)%i_fe.le.n_el ) then
         do i1 = 1, frm_cnf(iw)%i_fe-1
            i2 = n_el * (i1-1)-i1 * (i1+1) / 2 + frm_cnf(iw)%i_fe
            v_ee(i2) = 1.0_dp / frm_cnf(iw)%d_ee(i2)%m
         enddo 
         i2 = n_el * (frm_cnf(iw)%i_fe-1)-frm_cnf(iw)%i_fe* (frm_cnf(iw)%i_fe+1) / 2 + frm_cnf(iw)%i_fe
         do i1 = frm_cnf(iw)%i_fe+1, n_el
            i2 = i2 + 1_int32
            v_ee(i2) = 1.0_dp / frm_cnf(iw)%d_ee(i2)%m
         enddo 
      endif
   end subroutine upd_v_ee
   subroutine upd_v_pp( iw, v_pp )
      integer(int32),                       intent(in)    :: iw
      real(dp), dimension(n_po*(n_po-1)/2), intent(inout) :: v_pp
      integer(int32) :: i1, i2, i3
      if ( frm_cnf(iw)%i_fe.gt.n_el ) then
         i3 = frm_cnf(iw)%i_fe - n_el
         do i1 = 1, i3-1
            i2 = n_po*(i1-1)-i1*(i1+1)/2 + i3
            v_pp(i2) = 1.0_dp / frm_cnf(iw)%d_pp(i2)%m
         enddo 
         i2 = n_po*(i3-1)-i3*(i3+1)/2 + i3
         do i1 = i3+1, n_po
            i2 = i2 + 1_int32
            v_pp(i2) = 1.0_dp / frm_cnf(iw)%d_pp(i2)%m
         enddo 
      endif
   end subroutine upd_v_pp
   subroutine upd_v_ep( iw, v_ep )
      integer(int32),                 intent(in)    :: iw
      real(dp), dimension(n_po,n_el), intent(inout) :: v_ep
      integer(int32) :: i1
      if (frm_cnf(iw)%i_fe.le.n_el) then
         do i1 = 1_int32, n_po
            v_ep(i1,frm_cnf(iw)%i_fe) = - 1.0_dp / frm_cnf(iw)%d_ep(i1,frm_cnf(iw)%i_fe)%m
         enddo 
      else
         do i1 = 1_int32, n_el
            v_ep(frm_cnf(iw)%i_fe-n_el,i1) = - 1.0_dp / frm_cnf(iw)%d_ep(frm_cnf(iw)%i_fe-n_el,i1)%m
         enddo 
      endif
   end subroutine upd_v_ep
end module fermionic_potentials_m
