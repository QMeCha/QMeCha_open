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
module cubic_harmonics_m
   use fortran_kinds_v,      only: dp, int32
   use physical_constants_v, only: pi, twopi
   use molecular_system_v,           only: dist_t
   use quantum_monte_carlo_v,           only: n_wlk_max
   implicit none
   integer(int32), save, public :: l_max
   real(dp),       save, public, allocatable, dimension(:) :: c_y_l
   type, public :: cubhrm_t
      real(dp), allocatable, dimension(:)   ::   y_l
      real(dp), allocatable, dimension(:,:) ::  Dy_l
   contains
      procedure :: cmp   => cmp_cubic_hrm
      procedure :: cmp_D => cmp_D_cubic_hrm
   end type cubhrm_t
   type(cubhrm_t), save, public, allocatable, dimension(:), target :: y
   real(dp),       save, public, allocatable, dimension(:,:,:) :: y_l_sym_tab
   public  :: init_cubic_hrm
   private :: cmp_cubic_hrm, cmp_D_cubic_hrm, init_cbhrm_coefs, init_symm_cbhrm
contains
   subroutine init_cubic_hrm( )
      integer(int32) :: iw
      call init_cbhrm_coefs()
      call init_symm_cbhrm()
      allocate( y(1:n_wlk_max) )
      do iw = 1_int32, n_wlk_max
         allocate(y(iw)%y_l(0:l_max*(l_max+2))) ; y(iw)%y_l = 0.0_dp
         allocate(y(iw)%Dy_l(1:3,0:l_max*(l_max+2))) ; y(iw)%Dy_l = 0.0_dp
      enddo
   end subroutine init_cubic_hrm
   subroutine init_cbhrm_coefs( )
      allocate(c_y_l(0:l_max*(l_max+2)))
      c_y_l(0) = 1.0_dp / dsqrt(2.0_dp * twopi)
      if( l_max.gt.0 ) then
         c_y_l(1:3) = dsqrt(3.0_dp / ( 2.0_dp * twopi ))
      endif
      if( l_max.gt.1 ) then
         c_y_l(4) = 0.50_dp * dsqrt( 15.0_dp / pi )
         c_y_l(5) = 0.50_dp * dsqrt( 15.0_dp / pi )
         c_y_l(6) = 0.25_dp * dsqrt(  5.0_dp / pi )
         c_y_l(7) = 0.50_dp * dsqrt( 15.0_dp / pi )
         c_y_l(8) = 0.25_dp * dsqrt( 15.0_dp / pi )
      endif
      if ( l_max.gt.2 ) then
         c_y_l(9)  = 0.25_dp * dsqrt(  35.0_dp / twopi )
         c_y_l(10) = 0.50_dp * dsqrt( 105.0_dp / pi )
         c_y_l(11) = 0.25_dp * dsqrt(  21.0_dp / twopi )
         c_y_l(12) = 0.25_dp * dsqrt(   7.0_dp / pi )
         c_y_l(13) = 0.25_dp * dsqrt(  21.0_dp / twopi )
         c_y_l(14) = 0.25_dp * dsqrt( 105.0_dp / pi )
         c_y_l(15) = 0.25_dp * dsqrt(  35.0_dp / twopi )
      endif
      if ( l_max.gt.3 ) then
         c_y_l(16) = 0.75000_dp * dsqrt( 35.0_dp / pi )
         c_y_l(17) = 0.75000_dp * dsqrt( 35.0_dp / twopi )
         c_y_l(18) = 0.75000_dp * dsqrt(  5.0_dp / pi )
         c_y_l(19) = 0.75000_dp * dsqrt(  5.0_dp / twopi )
         c_y_l(20) = 0.18750_dp * dsqrt(  1.0_dp / pi )
         c_y_l(21) = 0.75000_dp * dsqrt(  5.0_dp / twopi )
         c_y_l(22) = 0.37500_dp * dsqrt(  5.0_dp / pi )
         c_y_l(23) = 0.75000_dp * dsqrt( 35.0_dp / twopi )
         c_y_l(24) = 0.18750_dp * dsqrt( 35.0_dp / pi )
      endif
   end subroutine init_cbhrm_coefs
   subroutine init_symm_cbhrm()
      real(dp)  :: s3, s5, s7
      allocate(y_l_sym_tab(1:4,0:4,0:l_max*(l_max+2))) ; y_l_sym_tab = 0.0_dp
      s3 = dsqrt(3.0_dp) ; s5 = dsqrt(5.0_dp) ; s7 = dsqrt(7.0_dp)
      y_l_sym_tab(:,0, 0) = (/ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)
      if ( l_max.gt.0 ) then
         y_l_sym_tab(:,0, 1) = (/ 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0, 2) = (/ 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0, 3) = (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)
      endif
      if ( l_max.gt.1 ) then
         y_l_sym_tab(:,0, 4) = (/ 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0, 5) = (/ 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0, 6) = (/ 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0, 7) = (/ 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0, 8) = (/ 2.0_dp,-2.0_dp, 0.0_dp, 0.0_dp /)
      endif
      if ( l_max.gt.2 ) then
         y_l_sym_tab(:,0, 9) = (/ 2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,10) = (/ 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,11) = (/ 0.0_dp, 1.0_dp, 2.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,12) = (/ 0.0_dp, 0.0_dp, 3.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,13) = (/ 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,14) = (/ 2.0_dp,-2.0_dp, 1.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,15) = (/ 1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp /)
      endif
      if ( l_max.gt.3 ) then
         y_l_sym_tab(:,0,16) = (/ 3.0_dp,-3.0_dp, 0.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,17) = (/ 0.0_dp, 3.0_dp, 1.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,18) = (/ 1.0_dp, 1.0_dp, 2.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,19) = (/ 0.0_dp, 1.0_dp, 3.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,20) = (/ 0.0_dp, 0.0_dp, 4.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,21) = (/ 1.0_dp, 0.0_dp, 3.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,22) = (/ 2.0_dp,-2.0_dp, 2.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,23) = (/ 3.0_dp, 0.0_dp, 1.0_dp, 0.0_dp /)
         y_l_sym_tab(:,0,24) = (/ 2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp /)
      endif
      y_l_sym_tab(:,1,0) = (/ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)
      if ( l_max.gt.0 ) then
         y_l_sym_tab(:,1,1) = (/ 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp /)
         y_l_sym_tab(:,1,2) = (/ 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp /)
         y_l_sym_tab(:,1,3) = (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)
      endif
      if ( l_max.gt.1 ) then
         y_l_sym_tab(:,1:2, 4) = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp /), (/4,2/))
         y_l_sym_tab(:,1:2, 5) = reshape( (/ 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp /), (/4,2/))
         y_l_sym_tab(:,1:2, 6) = reshape( (/ 0.0_dp, 0.0_dp,    s3, 1.0_dp, &
         & 0.0_dp, 0.0_dp,    s3,-1.0_dp /), (/4,2/))
         y_l_sym_tab(:,1:2, 7) = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp /), (/4,2/))
         y_l_sym_tab(:,1:2, 8) = reshape( (/ 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 1.0_dp,-1.0_dp, 0.0_dp, 0.0_dp /), (/4,2/))
      endif
      if ( l_max.gt.2 ) then
         y_l_sym_tab(:,1:3, 9) = reshape( (/ 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         &    s3, 1.0_dp, 0.0_dp, 0.0_dp, &
         &    s3,-1.0_dp, 0.0_dp, 0.0_dp /), (/4,3/))
         y_l_sym_tab(:,1:3,10) = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp /), (/4,3/))
         y_l_sym_tab(:,1:3,11) = reshape( (/ 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp,    s5, 1.0_dp, &
         & 0.0_dp, 0.0_dp,    s5,-1.0_dp /), (/4,3/))
         y_l_sym_tab(:,1:3,12) = reshape( (/ 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp,    s5,    s3, &
         & 0.0_dp, 0.0_dp,    s5,   -s3 /), (/4,3/))
         y_l_sym_tab(:,1:3,13) = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp,    s5, 1.0_dp, &
         & 0.0_dp, 0.0_dp,    s5,-1.0_dp /), (/4,3/))
         y_l_sym_tab(:,1:3,14) = reshape( (/ 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 1.0_dp,-1.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp /), (/4,3/))
         y_l_sym_tab(:,1:3,15) = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         & 1.0_dp,    s3, 0.0_dp, 0.0_dp, &
         & 1.0_dp,   -s3, 0.0_dp, 0.0_dp /), (/4,3/))
      endif
      if ( l_max.gt.3 ) then
         y_l_sym_tab(:,1:4,16) = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 1.0_dp,-1.0_dp, 0.0_dp, 0.0_dp /), (/4,4/))
         y_l_sym_tab(:,1:4,17) = reshape( (/ 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
         &    s3, 1.0_dp, 0.0_dp, 0.0_dp, &
         &    s3,-1.0_dp, 0.0_dp, 0.0_dp /), (/4,4/))
         y_l_sym_tab(:,1:4,18) = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp,    s7, 1.0_dp, &
         & 0.0_dp, 0.0_dp,    s7,-1.0_dp /), (/4,4/))
         y_l_sym_tab(:,1:4,19) = reshape( (/ 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp,    s7,    s3, &
         & 0.0_dp, 0.0_dp,    s7,   -s3 /), (/4,4/))
         y_l_sym_tab(:,1:4,20) = reshape( (/ 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp,-1.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp,-1.0_dp /), (/4,4/))
         y_l_sym_tab(:,1:4,21) = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp,    s7,    s3, &
         & 0.0_dp, 0.0_dp,    s7,   -s3 /), (/4,4/))
         y_l_sym_tab(:,1:4,22) = reshape( (/ 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 1.0_dp,-1.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp,    s7, 1.0_dp, &
         & 0.0_dp, 0.0_dp,    s7,-1.0_dp /), (/4,4/))
         y_l_sym_tab(:,1:4,23) = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         & 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
         & 1.0_dp,    s3, 0.0_dp, 0.0_dp, &
         & 1.0_dp,   -s3, 0.0_dp, 0.0_dp /), (/4,4/))
         y_l_sym_tab(:,1:4,24) = reshape( (/ 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 1.0_dp,-1.0_dp, 0.0_dp, 0.0_dp, &
         & 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         & 1.0_dp,-1.0_dp, 0.0_dp, 0.0_dp /), (/4,4/))
      endif
   end subroutine init_symm_cbhrm
   subroutine cmp_cubic_hrm( obj, d_fn, l_max_in )
      class(cubhrm_t), intent(inout) :: obj
      type(dist_t),    intent(in)    :: d_fn
      integer(int32),  intent(in)    :: l_max_in
      real(dp) :: x, y, z, x2, y2, z2, r2
      obj%y_l(0) = 1.0_dp
      if( l_max_in.gt.0_int32 ) then
         x = d_fn%v(1) ; y = d_fn%v(2) ; z = d_fn%v(3)
         obj%y_l(1) = y
         obj%y_l(2) = z
         obj%y_l(3) = x
         if( l_max_in.gt.1_int32 ) then
            x2 = x**2
            y2 = y**2
            z2 = z**2
            r2 = d_fn%m**2
            obj%y_l(4) = x * y
            obj%y_l(5) = y * z
            obj%y_l(6) = 3.0_dp * z2 - r2
            obj%y_l(7) = x * z
            obj%y_l(8) = x2 - y2
            if ( l_max_in.gt.2_int32 ) then
               obj%y_l( 9) = y * ( 3.0_dp * x2 - y2 )
               obj%y_l(10) = z * x * y
               obj%y_l(11) = y * ( 5.0_dp * z2 - r2 )
               obj%y_l(12) = z * ( 5.0_dp * z2 - 3.0_dp * r2 )
               obj%y_l(13) = x * ( 5.0_dp * z2 - r2 )
               obj%y_l(14) = z * ( x2 - y2 )
               obj%y_l(15) = x * ( x2 - 3.0_dp * y2 )
               if ( l_max_in.gt.3_int32 ) then
                  obj%y_l(16) = x * y * ( x2 - y2 )
                  obj%y_l(17) = z * y * ( 3.0_dp * x2 - y2 )
                  obj%y_l(18) = x * y * ( 7.0_dp * z2 -r2 )
                  obj%y_l(19) = y * z * ( 7.0_dp * z2 -3.0_dp* r2 )
                  obj%y_l(20) = ( 35.0_dp * z2 **2 - 30.0_dp * z2 * r2 + 3.0_dp * r2**2 )
                  obj%y_l(21) = x * z * ( 7.0_dp * z2 -3.0_dp* r2 )
                  obj%y_l(22) = (x2 - y2) * ( 7.0_dp * z2 - r2 )
                  obj%y_l(23) = z * x * ( x2 - 3.0_dp * y2 )
                  obj%y_l(24) = ( x2 * ( x2 - 3.0_dp * y2 ) - y2 * ( 3.0_dp * x2 - y2 ) )
               endif
            endif
         endif
      endif
      obj%y_l(0:l_max_in*(l_max_in+2)) = obj%y_l(0:l_max_in*(l_max_in+2)) * c_y_l(0:l_max_in*(l_max_in+2))
   end subroutine cmp_cubic_hrm
   subroutine cmp_D_cubic_hrm( obj, d_fn, l_max_in )
      class(cubhrm_t), intent(inout) :: obj
      type(dist_t),    intent(in)    :: d_fn
      integer(int32),  intent(in)    :: l_max_in
      real(dp)        :: x, y, z, x2, y2, z2
      integer(int32)  :: i1
      obj%Dy_l(1:3,0) = 0.0_dp
      if( l_max_in.gt.0_int32 ) then
         obj%Dy_l(1:3,1) = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
         obj%Dy_l(1:3,2) = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
         obj%Dy_l(1:3,3) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
         if( l_max_in.gt.1_int32 ) then
            x = d_fn%v(1) ; y = d_fn%v(2) ; z = d_fn%v(3)
            obj%Dy_l(1:3,4) = (/         y,          x,   0.0_dp /)
            obj%Dy_l(1:3,5) = (/    0.0_dp,          z,        y /)
            obj%Dy_l(1:3,6) = (/ -2.0_dp*x,  -2.0_dp*y, 4.0_dp*z /)
            obj%Dy_l(1:3,7) = (/         z,     0.0_dp,        x /)
            obj%Dy_l(1:3,8) = (/ 2.0_dp*x,  -2.0_dp*y,   0.0_dp /)
            if( l_max_in.gt.2_int32 ) then
               x2 = x**2
               y2 = y**2
               z2 = z**2
               obj%Dy_l(1:3,9)  = (/  6.0_dp*x*y            ,  3.0_dp*(x2-y2)        , 0.0_dp                   /)
               obj%Dy_l(1:3,10) = (/         y*z            ,         x*z            ,        x*y               /)
               obj%Dy_l(1:3,11) = (/ -2.0_dp*x*y            ,  4.0_dp*z2-3.0_dp*y2-x2, 8.0_dp*y*z               /)
               obj%Dy_l(1:3,12) = (/ -6.0_dp*x*z            , -6.0_dp*y*z            , 3.0_dp*(2.0_dp*z2-x2-y2) /)
               obj%Dy_l(1:3,13) = (/  4.0_dp*z2-3.0_dp*x2-y2, -2.0_dp*x*y            , 8.0_dp*x*z               /)
               obj%Dy_l(1:3,14) = (/  2.0_dp*x*z            , -2.0_dp*y*z            , x2-y2                    /)
               obj%Dy_l(1:3,15) = (/  3.0_dp*(x2-y2)        , -6.0_dp*x*y            , 0.0_dp                   /)
               if( l_max_in.gt.3_int32 ) then
                  obj%Dy_l(1:3,16) = (/ y*(3.0_dp*x2-y2) ,x*(x2-3.0_dp*y2)  ,0.0_dp /)
                  obj%Dy_l(1:3,17) = (/  6.0_dp*x*y*z, 3.0_dp*z*(x2-y2),   y*(3.0_dp*x2-y2) /)
                  obj%Dy_l(1:3,18) = (/ y*(6.0_dp*z2-3.0_dp*x2-y2), x*(6.0_dp*z2-3.0_dp*y2-x2), 12.0_dp*x*y*z /)
                  obj%Dy_l(1:3,19) = (/ -6.0_dp*x*y*z, z*(4.0_dp*z2-9.0_dp*y2-3.0_dp*x2), 3.0_dp*y*(4.0_dp*z2-x2-y2) /)
                  obj%Dy_l(1:3,20) = (/ 12.0_dp*x*(x2+y2-4.0_dp*z2), 12.0_dp*y*(x2+y2-4.0_dp*z2), 8.0_dp*(4.0_dp*z2-6.0_dp*x2-6.0_dp*y2) /)
                  obj%Dy_l(1:3,21) = (/ z*(4.0_dp*z2-9.0_dp*x2-3.0_dp*y2), -6.0_dp*x*y*z, 3.0_dp*x*(4.0_dp*z2-x2-y2)  /)
                  obj%Dy_l(1:3,22) = (/ 4.0_dp*x*(3.0_dp*z2-x2), 4.0_dp*y*(y2-3.0_dp*z2), 12.0_dp*z*(x2-y2)  /)
                  obj%Dy_l(1:3,23) = (/ 3.0_dp*z*(x2-y2), -6.0_dp*x*y*z  , x*(x2-3.0_dp*y2)  /)
                  obj%Dy_l(1:3,24) = (/ 4.0_dp*x*(x2-3.0_dp*y2), 4.0_dp*y*(y2-3.0_dp*x2), 0.0_dp /)
               endif
            endif
         endif
      endif
      do i1 = 0 , l_max_in * ( l_max_in + 2_int32 )
         obj%Dy_l(1:3,i1) = c_y_l(i1) * obj%Dy_l(1:3,i1)
      enddo
   end subroutine cmp_D_cubic_hrm
end module cubic_harmonics_m
