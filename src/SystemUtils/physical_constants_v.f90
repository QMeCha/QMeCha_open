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
module physical_constants_v
   use fortran_kinds_v, only : dp
   implicit none
   real(dp), public, parameter :: eul    = 0.5772156649015328606065120900824024310422_dp
   real(dp), public, parameter :: pi     = 3.141592653589793238462643383279502884197_dp
   real(dp), public, parameter :: pi_two = 1.57079632679489661923132169163975144209858_dp
   real(dp), public, parameter :: twopi  = 6.283185307179586476925286766559005768394_dp
   real(dp), public, parameter :: Ha_to_eV  = 27.211386245988_dp    ! (53) eV
   real(dp), public, parameter :: eV_to_Ha  =  0.036749322175655_dp ! (71) Ha
   real(dp), public, parameter :: au_to_ang =  0.529177210903_dp ! (80) Ang
   real(dp), public, parameter :: ang_to_au =  1.88972612462577_dp ! (80) au
   real(dp), public, parameter :: debye_to_au =  0.393430_dp
   real(dp), public, parameter :: au_to_debye =  2.541748_dp
end module physical_constants_v
