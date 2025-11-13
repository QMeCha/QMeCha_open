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
module primorbs_operations_m
   use fortran_kinds_v,      only: int32, dp
   use physical_constants_v, only: pi
   implicit none
   real(dp), private, save, allocatable, dimension(:) :: N_g
   real(dp), private, save, allocatable, dimension(:) :: N_s
   real(dp), private, save, allocatable, dimension(:) :: N_m
   real(dp), private, save, allocatable, dimension(:) :: N_e
   real(dp), private, parameter :: eps = 300_dp
   real(dp), private, parameter :: eta = 1.0_dp
   public  :: cmp_orb_ang, cmp_orb_type, cmp_prim_orb, cmp_GF_fact, cmp_Z_fact,&
   & ini_orb_nrm, prnt_orb_type, prnt_orb_ang
   private :: ini_S_norm, ini_G_norm, ini_M_norm, cmp_S_prim, cmp_G_prim, cmp_M_prim, &
   & ini_E_norm, cmp_E_prim
contains
   subroutine cmp_GF_fact( r, n, z, prm_t, G, F )
      real(dp),       intent(in)  :: r
      integer(int32), intent(in)  :: n
      real(dp),       intent(in)  :: z
      integer(int32), intent(in)  :: prm_t
      real(dp),       intent(out) :: G
      real(dp),       intent(out) :: F
      real(dp)                    :: dr_a, dr2_a
      select case ( prm_t )
       case (0,10)
         dr2_a = -2.0_dp * z
         dr_a  = dr2_a * r
       case (1,11)
         dr2_a = 0.0_dp
         dr_a = -z
       case (2,12)
         dr_a  = - r*(2.0_dp*eta+z*r)*(z/(eta+z*r))**2
         dr2_a = - 2.0_dp *(z*eta)**2 / (eta+z*r)**3
       case (3,13)
         dr_a  = - z**2*r/(1+z*r)
         dr2_a = - z**2/(1+z*r)**2
      end select
      if (n.gt.0_int32) then
         G = dble(n) / r + dr_a
         F = G**2 + dr2_a - dble(n) / r**2
      else
         G = dr_a
         F = G**2 + dr2_a
      endif
   end subroutine cmp_GF_fact
   subroutine cmp_Z_fact( r, n, l, z, prm_t, Zf )
      real(dp),       intent(in)  :: r
      integer(int32), intent(in)  :: n
      integer(int32), intent(in)  :: l
      real(dp),       intent(in)  :: z
      integer(int32), intent(in)  :: prm_t
      real(dp),       intent(out) :: Zf
      select case ( prm_t )
       case (0,10)
         Zf = - r**2
       case (1,11)
         Zf = - r
       case (2,12)
         Zf = - z*(2.0_dp * eta + z*r)*(r/(eta+z*r))**2
       case (3,13)
         Zf = - z*r**2/(1+z*r)
      end select
      select case ( prm_t )
       case (0)
         Zf = Zf + (2.0_dp * dble(n+l) + 3.0_dp) / ( 4.0_dp * z )
       case (1,2,3)
         Zf = Zf + (2.0_dp * dble(n+l) + 3.0_dp) / ( 2.0_dp * z )
      end select
   end subroutine cmp_Z_fact
   function cmp_prim_orb( r, z, l, prm_n, prm_t ) result( p )
      real(dp),       intent(in) :: r
      real(dp),       intent(in) :: z
      integer(int32), intent(in) :: l
      integer(int32), intent(in) :: prm_n
      integer(int32), intent(in) :: prm_t
      real(dp)                   :: p
      select case ( prm_t )
       case (0)
         call cmp_G_prim(r, z, l, prm_n, .true., p )
       case (1)
         call cmp_S_prim(r, z, l, prm_n, .true., p )
       case (2)
         call cmp_M_prim(r, z, l, prm_n, .true., p )
       case (3)
         call cmp_E_prim(r, z, l, prm_n, .true., p )
       case (10)
         call cmp_G_prim(r, z, l, prm_n, .false., p )
       case (11)
         call cmp_S_prim(r, z, l, prm_n, .false., p )
       case (12)
         call cmp_M_prim(r, z, l, prm_n, .false., p )
       case (13)
         call cmp_E_prim(r, z, l, prm_n, .false., p )
      end select
   end function cmp_prim_orb
   subroutine cmp_orb_type( orb_name, norm, orb_t, orb_n )
      character(2),   intent(in)  :: orb_name
      logical,        intent(in)  :: norm
      integer(int32), intent(out) :: orb_t
      integer(int32), intent(out) :: orb_n
      read(orb_name(1:1),*) orb_n
      orb_n = orb_n - 1_int32
      select case (orb_name(2:2))
       case ('G','g') ! Gaussian type primitve e^{-zr^2}
         orb_t =  0_int32
       case ('S','s') ! Slater type primitve e^{-zr}
         orb_t =  1_int32
       case ('M','m') ! Mixed type primitive e^{-(zr)^2/(1+zr)}
         orb_t =  2_int32
       case ('E','e') ! Mixed type primitive (1+zr)e^{-zr}
         orb_t =  3_int32
       case default
         print *,"Error!!! : Orbital type not recognized in basis set"
         stop
      end select
      if (.not.norm ) orb_t = orb_t + 10_int32
   end subroutine cmp_orb_type
   subroutine prnt_orb_type( orb_t, orb_n, orb_name )
      integer(int32), intent(in)  :: orb_t
      integer(int32), intent(in)  :: orb_n
      character(2),   intent(out) :: orb_name
      write(orb_name(1:1),'(I1)') orb_n + 1_int32
      select case (orb_t)
       case (0,10)
         orb_name(2:2) = 'G'
       case (1,11)
         orb_name(2:2) = 'S'
       case (2,12)
         orb_name(2:2) = 'M'
       case (3,13)
         orb_name(2:2) = 'E'
       case default
         print *,"Error!!! : Orbital type not recognized in basis set"
         stop
      end select
   end subroutine prnt_orb_type
   subroutine ini_orb_nrm( orb_t )
      integer(int32), intent(in) :: orb_t
      select case (orb_t)
       case (0)
         call ini_G_norm()
       case (1)
         call ini_S_norm()
       case (2)
         call ini_M_norm()
       case (3)
         call ini_E_norm()
       case default
      end select
   end subroutine ini_orb_nrm
   function cmp_orb_ang( orb_ang_mom ) result( orb_l )
      character(1),  intent(in) :: orb_ang_mom
      integer(int32)             :: orb_l
      select case (trim(orb_ang_mom))
       case ('S','s')
         orb_l =  0_int32
       case ('P','p')
         orb_l =  1_int32
       case ('D','d')
         orb_l =  2_int32
       case ('F','f')
         orb_l =  3_int32
       case ('G','g')
         orb_l =  4_int32
       case default
         print *,"Error!!! : Orbital Angular momentum not recognized"
         stop
      end select
   end function cmp_orb_ang
   function prnt_orb_ang( orb_l ) result( orb_ang_mom )
      integer(int32), intent(in) :: orb_l
      character(1)               :: orb_ang_mom
      select case (orb_l)
       case (0)
         orb_ang_mom = 'S'
       case (1)
         orb_ang_mom = 'P'
       case (2)
         orb_ang_mom = 'D'
       case (3)
         orb_ang_mom = 'F'
       case (4)
         orb_ang_mom = 'G'
       case default
         print *,"Error!!! : Orbital Angular momentum not recognized"
         stop
      end select
   end function prnt_orb_ang
   subroutine ini_S_norm()
      integer(int32) :: i1
      if ( .not.allocated(N_s) ) then
         allocate( N_s(0:14) )
         N_s(0) = 2.0_dp
         do i1 = 1_int32, 14_int32
            N_s(i1) = N_s(i1-1) * 2.0_dp &
            & / dsqrt(dble(( 2_int32 * i1 + 2_int32 ) * ( 2_int32 * i1 + 1_int32 )))
         enddo
      endif
   end subroutine ini_S_norm
   subroutine ini_G_norm()
      integer(int32) :: i1
      if ( .not.allocated(N_g) ) then
         allocate( N_g(0:14) )
         N_g(0) = 2.0_dp * ( 8.0_dp / pi ) ** 0.25_dp
         do i1 = 1_int32, 14_int32
            N_g(i1) = N_g(i1-1) * 2.0_dp / dsqrt(dble( 2_int32 * i1 + 1_int32 ))
         enddo
      endif
   end subroutine ini_G_norm
   subroutine ini_M_norm()
      if ( .not.allocated(N_m) ) then
         allocate( N_m(0:4) )
         if ( abs(1.0_dp-eta) .le. 1.0d-16 ) then
            N_m(0) = 1.126467421_dp
            N_m(1) = 0.576609950_dp
            N_m(2) = 0.196581141_dp
            N_m(3) = 0.050275655_dp
            N_m(4) = 0.010280772_dp
         elseif ( abs(0.5_dp-eta) .le. 1.0d-16 ) then
            N_m(0) = 1.402524252_dp
            N_m(1) = 0.769670200_dp
            N_m(2) = 0.274074953_dp
            N_m(3) = 0.072177594_dp
            N_m(4) = 0.015070484_dp
         elseif ( abs(0.1_dp-eta) .le. 1.0d-16 ) then
            N_m(0) = 1.825237677_dp
            N_m(1) = 1.049727352_dp
            N_m(2) = 0.382737444_dp
            N_m(3) = 0.102211707_dp
            N_m(4) = 0.021537912_dp
         elseif ( abs(0.01_dp-eta) .le. 1.0d-16 ) then
            N_m(0) = 1.980293982_dp
            N_m(1) = 1.143267860_dp
            N_m(2) = 0.417455522_dp
            N_m(3) = 0.111568758_dp
            N_m(4) = 0.023520642_dp
         endif
      endif
   end subroutine ini_M_norm
   subroutine ini_E_norm()
      integer(int32) :: i1
      if ( .not.allocated(N_e) ) then
         allocate( N_e(0:14) )
         N_e(0) = 2.0**2.5_dp / sqrt(24.0)
         do i1 = 1_int32, 14_int32
            N_e(i1) = N_e(i1-1) * 2.0**(i1) / sqrt(dble(2*i1+3)*dble(2*i1+4))
         enddo
         do i1 = 0_int32, 14_int32
            N_e(i1) = N_e(i1) * sqrt(dble(2*i1+3) / dble(2*i1+7))
         enddo
      endif
   end subroutine ini_E_norm
   subroutine cmp_G_prim(r, z, l, n, norm, p)
      real(dp),       intent(in) :: r
      real(dp),       intent(in) :: z
      integer(int32), intent(in) :: l
      integer(int32), intent(in) :: n
      logical,        intent(in) :: norm
      real(dp),       intent(out) :: p
      integer(int32)              :: m
      real(dp)                    :: zr2
      m = n+l
      zr2 = z*r**2
      if ( zr2.le.eps ) then
         p = dexp( -zr2 )
      else
         p = 0.0_dp
         return
      endif
      if ( norm ) p = p * N_g(m) * z ** (0.25_dp * dble(2*m+3) )
      if ( n.gt.0_int32 ) p = p * r ** n
   end subroutine cmp_G_prim
   subroutine cmp_S_prim(r, z, l, n, norm, p)
      real(dp),       intent(in) :: r
      real(dp),       intent(in) :: z
      integer(int32), intent(in) :: l
      integer(int32), intent(in) :: n
      logical,        intent(in) :: norm
      real(dp),       intent(out) :: p
      integer(int32)              :: m
      real(dp)                    :: zr
      m = n+l
      zr = z*r
      if ( zr.le.eps ) then
         p = dexp( -zr )
      else
         p = 0.0_dp
         return
      endif
      if ( norm ) p = p * N_s(m) * z ** ( dble(m) + 1.5_dp )
      if ( n.gt.0_int32 ) p = p * r ** n
   end subroutine cmp_S_prim
   subroutine cmp_M_prim(r, z, l, n, norm, p)
      real(dp),       intent(in) :: r
      real(dp),       intent(in) :: z
      integer(int32),  intent(in) :: l
      integer(int32),  intent(in) :: n
      logical,        intent(in) :: norm
      real(dp),       intent(out) :: p
      real(dp)                    :: zr
      zr = ((z*r) ** 2) / (eta+z*r)
      if ( zr.le.eps ) then
         p = dexp( -zr )
      else
         p = 0.0_dp
         return
      endif
      if ( norm ) p = p * N_m(n+l) * z ** ( dble(n+l) + 1.5_dp )
      if ( n.gt.0_int32 ) p = p * r ** n
   end subroutine cmp_M_prim
   subroutine cmp_E_prim(r, z, l, n, norm, p)
      real(dp),       intent(in) :: r
      real(dp),       intent(in) :: z
      integer(int32), intent(in) :: l
      integer(int32), intent(in) :: n
      logical,        intent(in) :: norm
      real(dp),       intent(out) :: p
      real(dp)                    :: zr
      zr = z*r
      if ( zr.le.eps ) then
         p = (1+zr)*dexp( -zr )
      else
         p = 0.0_dp
         return
      endif
      if ( norm ) p = p * N_e(n+l) * z ** ( dble(n+l) + 1.5_dp )
      if ( n.gt.0_int32 ) p = p * r ** n
   end subroutine cmp_E_prim
end module primorbs_operations_m
