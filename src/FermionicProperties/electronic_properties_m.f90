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
module electronic_properties_m
   use fortran_kinds_v,       only: int32, dp
   use quantum_monte_carlo_v, only: n_wlk_max
   use molecular_system_v,    only: n_at, atoms, r_at, n_fe, f_chrg
   use fermionic_config_c,    only: frm_cnf
   implicit none
   logical, public ,save :: elec_prp, dip_mom, quad_mom, expec_dist, expec_pote
   real(dp), private, save, allocatable, dimension(:) :: dip_atom
   real(dp), private, save, allocatable, dimension(:) :: quad_atom
   real(dp), allocatable, dimension(:,:) :: dipole
   real(dp), allocatable, dimension(:,:) :: quadrupole
   public  :: init_mlt_mmn, cmp_mlt_mmn
contains
   subroutine init_mlt_mmn
      real(dp)        :: r2
      integer(int32)  :: i1
      if (dip_mom.or.quad_mom.or.expec_dist ) then
         elec_prp = .true.
      else
         elec_prp = .false.
      endif
      if (dip_mom) then
         allocate( dipole(1:3,1:n_wlk_max) ) ; dipole = 0.0_dp
      endif
      if (quad_mom) then
         allocate( quadrupole(1:6,1:n_wlk_max) ) ; quadrupole = 0.0_dp
      endif
      if (dip_mom) then
         allocate( dip_atom(1:3) ) ; dip_atom(:) = 0.0_dp
         do i1 = 1_int32, n_at
            dip_atom(1) = dip_atom(1) + dble(atoms(i1)%atm_z) * r_at(1,i1)
            dip_atom(2) = dip_atom(2) + dble(atoms(i1)%atm_z) * r_at(2,i1)
            dip_atom(3) = dip_atom(3) + dble(atoms(i1)%atm_z) * r_at(3,i1)
         end do
      endif
      if (quad_mom) then
         allocate( quad_atom(1:6) ) ; quad_atom(:) = 0.0_dp
         do i1 = 1_int32, n_at
            r2 = sum(r_at(1:3,i1)**2)
            quad_atom(1) = quad_atom(1)+dble(atoms(i1)%atm_z)*(3.0_dp * r_at(1,i1)**2 - r2)
            quad_atom(2) = quad_atom(2)+dble(atoms(i1)%atm_z)*(3.0_dp * r_at(2,i1)**2 - r2)
            quad_atom(3) = quad_atom(3)+dble(atoms(i1)%atm_z)*(3.0_dp * r_at(3,i1)**2 - r2)
            quad_atom(4) = quad_atom(4)+dble(atoms(i1)%atm_z)* 3.0_dp * r_at(1,i1) * r_at(2,i1)
            quad_atom(5) = quad_atom(5)+dble(atoms(i1)%atm_z)* 3.0_dp * r_at(1,i1) * r_at(3,i1)
            quad_atom(6) = quad_atom(6)+dble(atoms(i1)%atm_z)* 3.0_dp * r_at(2,i1) * r_at(3,i1)
         end do
      endif
   end subroutine init_mlt_mmn
   subroutine cmp_mlt_mmn( iw )
      integer(int32), intent(in) :: iw
      real(dp)       :: r_cc(3)
      integer(int32) :: i1
      if (dip_mom) then
         dipole(:,iw) = dip_atom
         do i1 = 1_int32, n_fe
            dipole(1,iw) = dipole(1,iw) + dble(f_chrg(i1)) * frm_cnf(iw)%r_fe(1,i1)
            dipole(2,iw) = dipole(2,iw) + dble(f_chrg(i1)) * frm_cnf(iw)%r_fe(2,i1)
            dipole(3,iw) = dipole(3,iw) + dble(f_chrg(i1)) * frm_cnf(iw)%r_fe(3,i1)
         enddo  
      endif
      if (quad_mom) then
         quadrupole(:,iw) = quad_atom
         r_cc = 0.0_dp
         if ( n_at .eq. 0_int32 ) then
            do i1 = 1_int32, n_fe
               r_cc(1:3) = r_cc(1:3) + frm_cnf(iw)%r_fe(1:3,i1)
            end do
            r_cc = r_cc / dble(n_fe)
         end if
         do i1 = 1_int32, n_fe
            !r2 = sum( (frm_cnf(iw)%r_fe(1:3,i1) - r_cc(1:3) )**2 )
            !quadrupole(1) = quadrupole(1) + dble(f_chrg(i1)) * (3.0_dp*(frm_cnf(iw)%r_fe(1,i1) - r_cc(1)  )**2 - r2)
            !quadrupole(2) = quadrupole(2) + dble(f_chrg(i1)) * (3.0_dp*(frm_cnf(iw)%r_fe(2,i1) - r_cc(2)  )**2 - r2)
            !quadrupole(3) = quadrupole(3) + dble(f_chrg(i1)) * (3.0_dp*(frm_cnf(iw)%r_fe(3,i1) - r_cc(3)  )**2 - r2)
            !quadrupole(4) = quadrupole(4) + dble(f_chrg(i1)) * 3.0_dp*( frm_cnf(iw)%r_fe(1,i1) - r_cc(1) )*(frm_cnf(iw)%r_fe(2,i1) - r_cc(2) )
            !quadrupole(5) = quadrupole(5) + dble(f_chrg(i1)) * 3.0_dp*( frm_cnf(iw)%r_fe(1,i1) - r_cc(1) )*(frm_cnf(iw)%r_fe(3,i1) - r_cc(3) )
            !quadrupole(6) = quadrupole(6) + dble(f_chrg(i1)) * 3.0_dp*( frm_cnf(iw)%r_fe(2,i1) - r_cc(2) )*(frm_cnf(iw)%r_fe(3,i1) - r_cc(3) )
            quadrupole(1,iw) = quadrupole(1,iw) + dble(f_chrg(i1)) * ( (frm_cnf(iw)%r_fe(1,i1) - r_cc(1)  )**2 )
            quadrupole(2,iw) = quadrupole(2,iw) + dble(f_chrg(i1)) * ( (frm_cnf(iw)%r_fe(2,i1) - r_cc(2)  )**2 )
            quadrupole(3,iw) = quadrupole(3,iw) + dble(f_chrg(i1)) * ( (frm_cnf(iw)%r_fe(3,i1) - r_cc(3)  )**2 )
            quadrupole(4,iw) = quadrupole(4,iw) + dble(f_chrg(i1)) * ( frm_cnf(iw)%r_fe(1,i1) - r_cc(1) )*(frm_cnf(iw)%r_fe(2,i1) - r_cc(2) )
            quadrupole(5,iw) = quadrupole(5,iw) + dble(f_chrg(i1)) * ( frm_cnf(iw)%r_fe(1,i1) - r_cc(1) )*(frm_cnf(iw)%r_fe(3,i1) - r_cc(3) )
            quadrupole(6,iw) = quadrupole(6,iw) + dble(f_chrg(i1)) * ( frm_cnf(iw)%r_fe(2,i1) - r_cc(2) )*(frm_cnf(iw)%r_fe(3,i1) - r_cc(3) )
         enddo  
      endif
   end subroutine cmp_mlt_mmn
end module electronic_properties_m
