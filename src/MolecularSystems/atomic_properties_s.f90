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
function comp_atomic_charge( atm_name ) result(atm_chrg)
   use fortran_kinds_v, only: int32
   implicit none
   character(2),   intent(in) :: atm_name
   integer(int32)             :: atm_chrg
   character(2)               :: red_atm_name
   logical,          external :: is_numeric
   if ( is_numeric( atm_name(2:2) ) ) then
      red_atm_name = atm_name(1:1)
   else
      red_atm_name = atm_name(1:2)
   endif
   select case ( red_atm_name )
    case('H')  ; atm_chrg =  1_int32
    case('He') ; atm_chrg =  2_int32
    case('Li') ; atm_chrg =  3_int32
    case('Be') ; atm_chrg =  4_int32
    case('B')  ; atm_chrg =  5_int32
    case('C')  ; atm_chrg =  6_int32
    case('N')  ; atm_chrg =  7_int32
    case('O')  ; atm_chrg =  8_int32
    case('F')  ; atm_chrg =  9_int32
    case('Ne') ; atm_chrg = 10_int32
    case('Na') ; atm_chrg = 11_int32
    case('Mg') ; atm_chrg = 12_int32
    case('Al') ; atm_chrg = 13_int32
    case('Si') ; atm_chrg = 14_int32
    case('P')  ; atm_chrg = 15_int32
    case('S')  ; atm_chrg = 16_int32
    case('Cl') ; atm_chrg = 17_int32
    case('Ar') ; atm_chrg = 18_int32
    case('K')  ; atm_chrg = 19_int32
    case('Ca') ; atm_chrg = 20_int32
    case('Sc') ; atm_chrg = 21_int32
    case('Ti') ; atm_chrg = 22_int32
    case('V')  ; atm_chrg = 23_int32
    case('Cr') ; atm_chrg = 24_int32
    case('Mn') ; atm_chrg = 25_int32
    case('Fe') ; atm_chrg = 26_int32
    case('Co') ; atm_chrg = 27_int32
    case('Ni') ; atm_chrg = 28_int32
    case('Cu') ; atm_chrg = 29_int32
    case('Zn') ; atm_chrg = 30_int32
    case('Ga') ; atm_chrg = 31_int32
    case('Ge') ; atm_chrg = 32_int32
    case('As') ; atm_chrg = 33_int32
    case('Se') ; atm_chrg = 34_int32
    case('Br') ; atm_chrg = 35_int32
    case('Kr') ; atm_chrg = 36_int32
    case('Xe') ; atm_chrg = 54_int32
    case('Gh') ; atm_chrg = 0_int32
   end select
end function comp_atomic_charge
subroutine init_atomic_properties()
   use fortran_kinds_v,    only: int32
   use molecular_system_v, only: n_at, atoms
   implicit none
   character(2)                  :: atom_name
   logical,             external :: is_numeric
   integer(int32),      external :: comp_atomic_charge
   integer(int32) :: i1
   do i1 = 1_int32, n_at
      if ( atoms(i1)%atm_name(1:1) .eq. '*') then
         atoms(i1)%i_pp = -1
         atom_name = atoms(i1)%atm_name(2:3)
      else
         atoms(i1)%i_pp = 0
         atom_name = atoms(i1)%atm_name(1:2)
      endif
      atoms(i1)%atm_z = comp_atomic_charge( atom_name )
   enddo 
end subroutine init_atomic_properties
