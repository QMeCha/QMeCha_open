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
module external_field_m
   use fortran_kinds_v,    only: dp, int32
   use molecular_system_v, only: n_at, n_el, n_po, n_fe, f_chrg, r_at, atoms, nuc_chrg
   use fermionic_config_c, only: frm_cnf
   use qdo_system_v,       only: n_qdo, n_pchrg, qdos, r_qdo, pchrgs, r_pchrg
   use drudonic_config_c,  only: drd_cnf
   implicit none
   real(dp), public, save, dimension(3) :: ext_E_field
   real(dp), public, save, dimension(3) :: R_c
   logical,  public, save               :: floating_r_c 
   real(dp), public, save               :: V_E_QDOc, V_E_Zac, V_E_Pc
   public :: init_extfld, cmp_extfld_frm, cmp_extfld_drd
contains
   subroutine init_extfld( v_ext )
      real(dp), intent(inout) :: v_ext
      real(dp) :: tot_fixed_charge
      integer(int32) :: i1 
      V_E_Zac  = 0.0_dp 
      V_E_QDOc = 0.0_dp
      V_E_Pc   = 0.0_dp
      R_c(1:3) = 0.0_dp
      if (n_at.gt.0_int32) then 
         do i1 = 1_int32, n_at 
            V_E_Zac = -atoms(i1)%atm_z * dot_product ( ext_E_field(1:3), r_at(1:3,i1) ) 
         enddo 
      endif
      if (n_qdo.gt.0_int32) then 
         do i1 = 1_int32, n_qdo 
            V_E_QDOc = -qdos(i1)%qdo_q * dot_product ( ext_E_field(1:3), r_qdo(1:3,i1) ) 
         enddo 
      endif 
      if (n_pchrg.gt.0_int32) then 
         do i1 = 1_int32, n_pchrg 
            V_E_Pc = -pchrgs(i1)%pchrg_q * dot_product ( ext_E_field(1:3), r_pchrg(1:3,i1) ) 
         enddo 
      endif 
      v_ext = v_ext + V_E_Zac + V_E_QDOc + V_E_Pc
      floating_r_c = .false. 
      if ( (n_qdo + n_at + n_pchrg).eq.0_int32) then
         if (n_fe.gt.0_int32) floating_r_c = .true.
      else 
         if ((nuc_chrg-n_el+n_po).ne.0_int32) then 
            tot_fixed_charge = 0.0_dp
            if (n_at.gt.0_int32) then 
               do i1 = 1_int32, n_at 
                  R_c(1:3) = R_c(1:3) + atoms(i1)%atm_z * r_at(1:3,i1) 
                  tot_fixed_charge = tot_fixed_charge + atoms(i1)%atm_z
               enddo ! (n_at)
            endif
            if (n_qdo.gt.0_int32) then 
               do i1 = 1_int32, n_qdo 
                  R_c(1:3) = R_c(1:3) + qdos(i1)%qdo_q * r_qdo(1:3,i1)  
                  tot_fixed_charge = tot_fixed_charge + qdos(i1)%qdo_q
               enddo ! (n_at)
            endif 
            if (n_pchrg.gt.0_int32) then 
               do i1 = 1_int32, n_pchrg 
                  R_c(1:3) = R_c(1:3) + pchrgs(i1)%pchrg_q * r_pchrg(1:3,i1)  
                  tot_fixed_charge = tot_fixed_charge + pchrgs(i1)%pchrg_q
               enddo ! (n_at)
            endif 
            v_ext = v_ext + tot_fixed_charge * dot_product ( ext_E_field(1:3), R_c(1:3) ) 
         endif
      endif
   end subroutine init_extfld
   subroutine cmp_extfld_frm( iw, e_pp_loc, v_loc )
      integer(int32),            intent(in)    :: iw
      real(dp), dimension(n_fe), intent(inout) :: e_pp_loc
      real(dp),                  intent(inout) :: v_loc
      real(dp)       :: tmp_cmp
      real(dp)       :: r_cc(3)
      integer(int32) :: i1
      if ( floating_r_c ) then
         do i1 = 1_int32, n_fe
            r_cc(1:3) = r_cc(1:3) + frm_cnf(iw)%r_fe(1:3,i1) 
         end do
         r_cc = r_cc / dble(n_fe)
         do i1 = 1_int32, n_fe
            tmp_cmp = -f_chrg(i1) * dot_product ( ext_E_field(1:3), frm_cnf(iw)%r_fe(1:3,i1)-r_cc(1:3) )
            e_pp_loc(i1) = e_pp_loc(i1) + tmp_cmp
            v_loc = v_loc + tmp_cmp
         end do
      else 
         do i1 = 1_int32, n_fe
            tmp_cmp = -f_chrg(i1) * dot_product ( ext_E_field(1:3), frm_cnf(iw)%r_fe(1:3,i1))
            e_pp_loc(i1) = e_pp_loc(i1) + tmp_cmp
            v_loc = v_loc + tmp_cmp
         end do
      endif
   end subroutine cmp_extfld_frm
   subroutine cmp_extfld_drd( iw, e_pp_loc, v_loc )
      integer(int32),             intent(in)    :: iw
      real(dp), dimension(n_qdo), intent(inout) :: e_pp_loc
      real(dp),                   intent(inout) :: v_loc
      real(dp) :: tmp_cmp
      integer(int32) :: i1
      do i1 = 1_int32, n_qdo
         tmp_cmp = tmp_cmp + qdos(i1)%qdo_q * dot_product( ext_E_field(1:3), drd_cnf(iw)%r_drd(1:3,i1))
         e_pp_loc(i1) = e_pp_loc(i1) + tmp_cmp 
         v_loc = v_loc + tmp_cmp
      enddo
   end subroutine cmp_extfld_drd
end module external_field_m
