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
module atomic_orbital_c
   use fortran_kinds_v,       only: dp, int32
   use cubic_harmonics_m,     only: cubhrm_t
   use molecular_system_v,    only: dist_t
   use primorbs_operations_m, only: cmp_prim_orb, cmp_GF_fact, cmp_Z_fact
   implicit none
   type, public :: orb_t
      integer(int32) :: i_orb, i_par, n_par
      integer(int32)  :: l, l_z, n
      integer(int32), allocatable, dimension(:) :: prm_n
      integer(int32), allocatable, dimension(:) :: prm_t
      real(dp),       allocatable, dimension(:) :: c
      real(dp),       allocatable, dimension(:) :: z
   contains
      procedure :: cmp      => cmp_orb
      procedure :: cmp_DD2  => cmp_DD2_orb
      procedure :: cmp_DD2s => cmp_DD2_sorb
      procedure :: cmp_der  => cmp_der_orb
      procedure :: cmp_dc   => cmp_dc_orb
   end type orb_t
   private :: cmp_orb, cmp_DD2_orb, cmp_DD2_sorb, cmp_der_orb, cmp_dc_orb
contains
   subroutine cmp_orb( obj, r, y, P_n, orb )
      class(orb_t), intent(inout) :: obj
      real(dp),                     intent(in)    :: r
      type(cubhrm_t),               intent(in)    :: y
      real(dp), dimension(obj%n),   intent(inout) :: P_n
      real(dp),                     intent(inout) :: orb
      integer(int32) :: i1
      do i1 = 1_int32, obj%n
         P_n(i1) = cmp_prim_orb( r, obj%z(i1), obj%l, obj%prm_n(i1), obj%prm_t(i1) )
      enddo
      P_n(:) = P_n(:) * y%y_l(obj%l_z)
      orb = dot_product( obj%c(:), P_n(:) )
   end subroutine cmp_orb
   subroutine cmp_DD2_orb( obj, r, y, n_orbs, P_n, orb, DD2_orb )
      class(orb_t), intent(inout) :: obj
      type(dist_t),                            intent(in)    :: r
      type(cubhrm_t),                          intent(in)    :: y
      integer(int32),                          intent(in)    :: n_orbs
      real(dp),       dimension(obj%n),      intent(in)    :: P_n
      real(dp),                                intent(in)    :: orb
      real(dp),       dimension(n_orbs,4), intent(inout) :: DD2_orb
      real(dp)       :: G_n, F_n, G, F
      integer(int32) :: i1
      G = 0.0_dp
      F = 0.0_dp
      do i1 = 1_int32, obj%n
         call cmp_GF_fact( r%m, obj%prm_n(i1), obj%z(i1), obj%prm_t(i1), G_n, F_n)
         G = G + obj%c(i1) * P_n(i1) * G_n
         F = F + obj%c(i1) * P_n(i1) * F_n
      enddo
      DD2_orb(obj%i_orb,1:3) = G * r%v(1:3) / r%m + orb * y%Dy_l(1:3,obj%l_z) / y%y_l(obj%l_z)
      DD2_orb(obj%i_orb,4) = F + dble(2*(obj%l+1)) * G / r%m
   end subroutine cmp_DD2_orb
   subroutine cmp_der_orb( obj, r, y, n_orbs, n_par_orbs, de_opt, P_n, orb, DD2_orb, da_orb )
      class(orb_t), intent(inout) :: obj
      type(dist_t),                            intent(in)    :: r
      type(cubhrm_t),                          intent(in)    :: y
      integer(int32),                          intent(in)    :: n_orbs
      integer(int32),                          intent(in)    :: n_par_orbs
      logical,                                 intent(in)    :: de_opt
      real(dp),       dimension(obj%n),      intent(in)    :: P_n
      real(dp),                                intent(in)    :: orb
      real(dp),       dimension(n_orbs,4), intent(inout) :: DD2_orb
      real(dp),       dimension(n_par_orbs), intent(inout) :: da_orb
      real(dp)       :: G_n, F_n, Z_n, G, F
      integer(int32) :: i1
      G = 0.0_dp
      F = 0.0_dp
      do i1 = 1_int32, obj%n
         call cmp_GF_fact(r%m, obj%prm_n(i1), obj%z(i1), obj%prm_t(i1), G_n, F_n)
         G = G + obj%c(i1) * P_n(i1) * G_n
         F = F + obj%c(i1) * P_n(i1) * F_n
      enddo
      DD2_orb(obj%i_orb,1:3) = G / r%m * r%v(1:3) &
      & + orb / y%y_l(obj%l_z) * y%Dy_l(1:3,obj%l_z)
      DD2_orb(obj%i_orb,4) = F + dble(2*(obj%l+1))/r%m * G
      if ( obj%n_par.gt.0_int32 ) then
         if ( obj%n.gt.1_int32 ) then
            if( de_opt ) then
               da_orb(obj%i_par:obj%i_par+obj%n-1) = P_n(1:obj%n)
               do i1 = 1_int32, obj%n
                  call cmp_Z_fact( r%m, obj%prm_n(i1), obj%l, obj%z(i1), obj%prm_t(i1), Z_n )
                  da_orb(obj%i_par+obj%n-1+i1) = obj%c(i1) * Z_n * P_n(i1)
               enddo
            else
               da_orb(obj%i_par:obj%i_par+obj%n-1) = P_n(1:obj%n)
            endif
         else
            if ( de_opt ) then
               call cmp_Z_fact( r%m, obj%prm_n(1), obj%l, obj%z(1), obj%prm_t(1), Z_n )
               da_orb(obj%i_par) = obj%c(1) * Z_n * P_n(1)
            endif
         endif 
      endif 
   end subroutine cmp_der_orb
   subroutine cmp_dc_orb( obj, r, n_par_orbs, de_opt, P_n, da_orb )
      class(orb_t), intent(inout) :: obj
      type(dist_t),                            intent(in)    :: r
      integer(int32),                          intent(in)    :: n_par_orbs
      logical,                                 intent(in)    :: de_opt
      real(dp),       dimension(obj%n),      intent(in)    :: P_n
      real(dp),       dimension(n_par_orbs), intent(inout) :: da_orb
      real(dp)       :: Z_n
      integer(int32) :: i1
      if ( obj%n_par.gt.0_int32 ) then
         if ( obj%n.gt.1_int32 ) then
            if( de_opt ) then
               da_orb(obj%i_par:obj%i_par+obj%n-1) = P_n(1:obj%n)
               do i1 = 1_int32, obj%n
                  call cmp_Z_fact( r%m, obj%prm_n(i1), obj%l, obj%z(i1), obj%prm_t(i1), Z_n )
                  da_orb(obj%i_par+obj%n-1+i1) = obj%c(i1) * Z_n * P_n(i1)
               enddo
            else
               da_orb(obj%i_par:obj%i_par+obj%n-1) = P_n(1:obj%n)
            endif
         else
            if ( de_opt ) then
               call cmp_Z_fact( r%m, obj%prm_n(1), obj%l, obj%z(1), obj%prm_t(1), Z_n )
               da_orb(obj%i_par) = obj%c(1) * Z_n * P_n(1)
            endif
         endif 
      endif
   end subroutine cmp_dc_orb
   subroutine cmp_DD2_sorb( obj, r, y, P_n, orb, DD2_orb )
      class(orb_t), intent(inout) :: obj
      type(dist_t),                            intent(in)    :: r
      type(cubhrm_t),                          intent(in)    :: y
      real(dp),       dimension(obj%n),      intent(in)    :: P_n
      real(dp),                                intent(in)    :: orb
      real(dp),       dimension(4),          intent(inout) :: DD2_orb
      real(dp)       :: G_n, F_n, G, F
      integer(int32) :: i1
      G = 0.0_dp
      F = 0.0_dp
      do i1 = 1_int32, obj%n
         call cmp_GF_fact(r%m, obj%prm_n(i1), obj%z(i1), obj%prm_t(i1), G_n, F_n)
         G = G + obj%c(i1) * P_n(i1) * G_n
         F = F + obj%c(i1) * P_n(i1) * F_n
      enddo
      DD2_orb(1:3) = G * r%v(1:3) / r%m + orb * y%Dy_l(1:3,obj%l_z) / y%y_l(obj%l_z)
      DD2_orb(4) = F + dble(2*(obj%l+1)) * G / r%m
   end subroutine cmp_DD2_sorb
end module atomic_orbital_c
