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
module cusp_functions_m
   use fortran_kinds_v,         only: dp, int32
   use cusps_functions_params_m,              only: jcspar_t
   use molecular_system_v,      only: dist_t, atoms
   implicit none
   real(dp), private, parameter :: eps = 600_dp
   public  :: cspfun, DD2_cspfun, db_cspfun
   private :: dtncps, DD2_dtncps, db_dtncps
contains
   pure function cspfun( G, c_indx, pseudo, b, d)
      real(dp),       intent(in) :: G
      integer(int32), intent(in) :: c_indx, pseudo
      type(jcspar_t), intent(in) :: b
      real(dp),       intent(in) :: d
      real(dp)                   :: const, b_tilde, zr
      real(dp)                   :: cspfun
      integer(int32) :: i1
      cspfun = 0.0_dp
      if (abs(G).ge.1.0_dp) then
         const = (2.0_dp * abs(G) )**0.25
      else
         const = 1.0_dp
      endif
      if ( const.lt.10.0d-16 ) then
         return
      endif
      if ( pseudo.eq.0_int32) then
         b_tilde = const * b%b(1,c_indx)
         select case(b%typ)
          case(1)
            cspfun = -G * exp( - b_tilde * d ) / b_tilde
          case(2)
            cspfun = -G / ( b_tilde*(1.0_dp + b_tilde * d) )
          case(3)
            cspfun = dtncps ( b%ord, G, b%b(1:b%ord+1,c_indx), d)
          case default
            cspfun = 0.0_dp
         end select
      else
         cspfun = 0.0_dp
      endif
      if ( b%ord.gt.0_int32  .and. b%typ.ne.3_int32 ) then
         do i1 = 1_int32, b%ord
            zr = b%b(b%ord+1+i1,c_indx) * d**2
            if ( zr.le.eps ) cspfun = cspfun + b%b(i1+1,c_indx) * exp( -zr )
         enddo
      endif 
   end function cspfun
   function DD2_cspfun( G, c_indx, pseudo, b, v, m )
      real(dp),               intent(in) :: G
      integer(int32),         intent(in) :: c_indx, pseudo
      type(jcspar_t),         intent(in) :: b
      real(dp), dimension(3), intent(in) :: v
      real(dp),               intent(in) :: m
      real(dp)                           :: DD2_cspfun(1:4)
      real(dp) :: const, b_tilde, zr
      real(dp) :: P_m, dr_F_r
      integer(int32) :: i1
      DD2_cspfun = 0.0_dp
      if (abs(G).ge.1.0_dp) then
         const = (2.0_dp * abs(G) )**0.25
      else
         const = 1.0_dp
      endif
      if ( const.lt.10.0d-16 ) then
         return
      endif
      if ( pseudo.eq.0_int32) then
         b_tilde = const * b%b(1,c_indx)
         select case(b%typ)
          case(1)
            P_m = -G * exp( - b_tilde * m ) / b_tilde
            dr_F_r  = -b_tilde * P_m
            DD2_cspfun(4) = (2.0_dp / m - b_tilde ) * dr_F_r
          case(2)
            P_m = -G / ( b_tilde*(1.0_dp + b_tilde * m) )
            dr_F_r  = -b_tilde * P_m / ( 1.0_dp + b_tilde * m)
            DD2_cspfun(4) = 2.0_dp*( 1.0_dp/m - b_tilde/( 1.0_dp + b_tilde * m ) ) * dr_F_r
          case(3)
            DD2_cspfun = DD2_dtncps(  b%ord, G, b%b(1:b%ord+1,c_indx), v, m)
          case default
            DD2_cspfun(4) = 0.0_dp
         end select
      else
         dr_F_r  = 0.0_dp
      endif
      if ( b%ord.gt.0_int32  .and. b%typ.ne.3_int32 ) then
         do i1 = 1_int32, b%ord
            zr = b%b(b%ord+1+i1,c_indx) * m**2
            if ( zr.le.eps ) then
               P_m = 2.0_dp*b%b(b%ord+1+i1,c_indx)*b%b(i1+1,c_indx)*exp(-zr)
               dr_F_r = dr_F_r - P_m*m
               DD2_cspfun(4) = DD2_cspfun(4) + P_m*(2.0_dp*b%b(b%ord+1+i1,c_indx)*m**2-3.0_dp)
            endif
         enddo
      endif 
      if (m.gt.0.0_dp .and. b%typ.ne.3_int32) DD2_cspfun(1:3) = dr_F_r/m*v(1:3)
   end function DD2_cspfun
   function db_cspfun( G, c_indx, pseudo, b, d, jce_opt )
      real(dp),         intent(in)    :: G
      integer(int32),   intent(in)    :: c_indx , pseudo
      type(jcspar_t),   intent(in)    :: b
      real(dp),         intent(in)    :: d
      logical,          intent(in)    :: jce_opt
      real(dp), dimension(b%n_opt_par) :: db_cspfun
      real(dp)                      :: const, b_tilde, zr
      integer(int32) :: i1
      db_cspfun = 0.0_dp
      if (abs(G).ge.1.0_dp) then
         const = (2.0_dp * abs(G) )**0.25
      else
         const = 1.0_dp
      endif
      if ( const.lt.10.0d-16 ) then
         return
      endif
      if ( pseudo.eq.0_int32) then
         b_tilde = const * b%b(1,c_indx)
         select case(b%typ)
          case(1)
            zr = b_tilde * d
            if ( zr.le.eps ) then
               db_cspfun(1) = G/b%b(1,c_indx)*(1.0_dp+zr)*exp(-zr)/b_tilde
            else
               db_cspfun(1) = 0.0_dp
            endif
          case(2)
            zr = b_tilde * d
            db_cspfun(1) = G/b%b(1,c_indx)*(1.0_dp+2.0_dp*zr)/(b_tilde*(1.0_dp+zr)**2)
          case(3)
            db_cspfun = db_dtncps(  b%ord, G, b%b(1:b%ord+1,c_indx), d)
          case default
            db_cspfun(1) = 0.0_dp
         end select
      endif
      if (  b%ord.gt.0_int32  ) then
         do i1 = 1_int32, b%ord
            zr = b%b(b%ord+i1+1,c_indx) * d**2
            if ( zr.le.eps ) then
               db_cspfun(1+i1) = exp(-zr)
            else
               db_cspfun(1+i1) = 0.0_dp
            endif
            if (jce_opt) db_cspfun(1+i1+b%ord) = - d**2 * b%b(i1+1,c_indx) * db_cspfun(i1+1)
         enddo
      endif 
   end function db_cspfun
   pure function dtncps( ord, G, b, d )
      integer(int32),           intent(in) :: ord
      real(dp), dimension(0:ord), intent(in) :: b
      real(dp),                 intent(in) :: G
      real(dp),                 intent(in) :: d
      real(dp)                             :: dtncps
      integer(int32) :: i1
      if ( d.ge.b(0) ) then
         dtncps = 0.0_dp
      else
         if (ord.gt.0) then
            dtncps = b(1) + (3.0_dp * b(1) / b(0) - G / b(0) ** 3 ) * d
         else
            dtncps = - G / b(0) ** 3 * d
         endif
         if (ord.gt.1) then
            do i1 = 2_int32, ord
               dtncps = dtncps + b(i1) * d ** i1
            enddo 
         endif
         dtncps = dtncps * ( d - b(0) ) ** 3
      endif
   end function dtncps
   pure function DD2_dtncps( ord, G, b, v, m )
      integer(int32),                 intent(in) :: ord
      real(dp), dimension(0:ord), intent(in) :: b
      real(dp),                       intent(in) :: G
      real(dp), dimension(3),       intent(in) :: v
      real(dp),                       intent(in) :: m
      real(dp), dimension(4)                   :: DD2_dtncps
      real(dp)                                   :: P_n
      real(dp)                                   :: F_n
      real(dp)                                   :: G_n
      integer(int32) :: i1
      if ( m.ge.b(0) ) then
         DD2_dtncps = 0.0_dp
      else
         G_n = 0.0_dp
         if (ord.gt.0) then
            F_n = ( 3.0_dp * b(1) / b(0) - G / b(0) ** 3 )
            P_n = b(1) + F_n * m
         else
            F_n = - G/ b(0) ** 3
            P_n = F_n * m
         endif
         if (ord.gt.1) then
            do i1 = 2_int32, ord
               P_n = P_n + b(i1) * m ** i1
               F_n = F_n + b(i1) * i1 * m ** ( i1 - 1 )
               G_n = G_n + b(i1) * i1 * ( i1 - 1 ) * m ** ( i1 - 2 )
            enddo ! i1 (ord)
         endif
         F_n  = F_n / P_n
         G_n = G_n / P_n - F_n ** 2
         F_n = F_n + 3.0_dp / (m - b(0))
         G_n = G_n - 3.0_dp / (m - b(0))**2
         G_n = G_n + F_n**2 + 2.0_dp / m * F_n
         DD2_dtncps(1:3) =  F_n * P_n * (m - b(0))**3 * v(1:3) / m
         DD2_dtncps(4)   =  G_n * P_n * (m - b(0))**3
      endif
   end function DD2_dtncps
   pure function db_dtncps( ord, G, b, d )
      integer(int32),                 intent(in) :: ord
      real(dp), dimension(0:ord), intent(in) :: b
      real(dp),                       intent(in) :: G
      real(dp),                       intent(in) :: d
      real(dp)                                  :: P_n
      real(dp), dimension(0:ord)             :: db_dtncps
      integer(int32) :: i1
      if ( d.ge.b(0) ) then
         db_dtncps = 0.0_dp
      else
         if ( ord.gt.0 ) then
            db_dtncps(0) = - 3.0_dp / b(0) **2 * (b(1) - G / b(0)**2 ) * d
            db_dtncps(1) = 1.0_dp + 3.0_dp / b(0) * d
            P_n = b(1) + (3.0_dp * b(1) / b(0) - G / b(0) ** 3 ) * d
         else
            db_dtncps(0) = 3.0_dp * G / b(0)**4 * d
            P_n = - G / b(0) ** 3  * d
         endif
         if ( ord.gt.1 ) then
            do i1 = 2_int32, ord
               db_dtncps(i1) = d ** i1
               P_n = P_n + b(i1) * db_dtncps(i1)
            enddo ! i1 (N_dtn)
         endif
         db_dtncps  = db_dtncps / P_n
         db_dtncps(0) = db_dtncps(0) - 3.0_dp / (d - b(0))
         db_dtncps =  db_dtncps * P_n * ( d - b(0) )**3
      endif
   end function db_dtncps
end module cusp_functions_m
