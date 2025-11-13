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
module nl_pseudo_grid_m
   use fortran_kinds_v,         only: dp, int32, stdout
   use physical_constants_v,    only: pi, twopi
   use openmp_mpi_m,            only: mpi_rank
   use write_lines_m,           only: write_variable_line, write_simple_line, write_empty_line
   use quantum_monte_carlo_v,              only: n_wlk_max, n_wlk
   use mersenne_twister19937_m, only: rndn
   use molecular_system_v,              only: r_at
   use fermionic_config_c,      only: frm_cnf
   use wavefunction_optimization_v,              only: corr_samp
   implicit none
   integer(int32), save, public                               :: n_psd_grd_pnts
   integer(int32), save, public                               :: l_max_grid
   logical,        save, public                               :: no_rot_grid
   real(dp),      save, public, allocatable, dimension(:,:)   :: r_grid
   real(dp),      save, public, allocatable, dimension(:)     :: w_grid
   real(dp),      save, public, allocatable, dimension(:,:)   :: angles
   real(dp),      save, public, allocatable, dimension(:,:,:) :: nlc_rnd_grid
   real(dp),      save, public, allocatable, dimension(:,:)   :: nlc_cos_theta
   public  :: init_nloc_grid, rndm_nlc_grd, P_l
   private :: octahedron_grid, icosahedron_grid
contains
   subroutine init_nloc_grid()
      real(dp)       :: const
      integer(int32) :: iw
      if ( corr_samp ) then
         no_rot_grid = .true.
         allocate( angles(1:3,1:n_wlk) ) ; angles(1:3,1:n_wlk) = 0.0_dp
         do iw = 1_int32, n_wlk
            angles(1,iw) = twopi*rndn(iw)%rndc()
            angles(2,iw) = twopi*rndn(iw)%rndc()
            angles(3,iw) = 1.0_dp-2.0_dp*rndn(iw)%rndc()
         enddo
      else
         no_rot_grid = .false.
      endif
      select case(n_psd_grd_pnts)
       case(1:3)
         n_psd_grd_pnts = 1_int32
       case(4:5)    ! Four points tetrahedral grid.
         n_psd_grd_pnts = 4_int32
       case(6:11)   ! Six points octahedral grid.
         n_psd_grd_pnts = 6_int32
       case(12:17)  ! 12 points icosahedral grid.
         n_psd_grd_pnts = 12_int32
       case(18:25)  ! 18 points octahedral grid.
         n_psd_grd_pnts = 18_int32
       case(26:31)  ! 26 points octahedral grid.
         n_psd_grd_pnts = 26_int32
       case(32:49)  ! 32 points icosahedral grid.
         n_psd_grd_pnts = 32_int32
       case(50:)    ! 50 points octahedral grid.
         n_psd_grd_pnts = 50_int32
      end select
      call write_variable_line(stdout,0,mpi_rank,2,"Int. points for non-local comp.", n_psd_grd_pnts,var_name="n_psd_grd_pnts")
      allocate( w_grid(1:n_psd_grd_pnts), r_grid(1:3,1:n_psd_grd_pnts) )
      w_grid = 0.0_dp ; r_grid = 0.0_dp
      allocate( nlc_rnd_grid(1:3,1:n_psd_grd_pnts,1:n_wlk_max) ) ; nlc_rnd_grid = 0.0_dp
      allocate( nlc_cos_theta(1:n_psd_grd_pnts,1:n_wlk_max) ) ; nlc_cos_theta = 0.0_dp
      select case(n_psd_grd_pnts)
       case(1)
         l_max_grid = 0_int32
         l_max_grid = 2_int32
         const = 1.0_dp / dsqrt(3.0_dp)
         r_grid(1:3,1) = (/  const,  const,  const /)
         r_grid(1:3,2) = (/  const, -const, -const /)
         r_grid(1:3,3) = (/ -const,  const, -const /)
         r_grid(1:3,4) = (/ -const, -const,  const /)
         w_grid(1:4)   = 0.25_dp
       case(6,18,26,50) ! Octahedral grids.
         call octahedron_grid( )
       case(12,32)      ! Icosahedral grids.
         call icosahedron_grid( )
      end select
      call write_variable_line(stdout,0,mpi_rank,2,"Maximum integrable angular momentum", l_max_grid,var_name="l_max_grid")
   end subroutine init_nloc_grid
   subroutine octahedron_grid()
      real(dp) :: const
      l_max_grid = 3_int32
      r_grid(1:3,1) = (/  1.0_dp,  0.0_dp,  0.0_dp /)
      r_grid(1:3,2) = (/ -1.0_dp,  0.0_dp,  0.0_dp /)
      r_grid(1:3,3) = (/  0.0_dp,  1.0_dp,  0.0_dp /)
      r_grid(1:3,4) = (/  0.0_dp, -1.0_dp,  0.0_dp /)
      r_grid(1:3,5) = (/  0.0_dp,  0.0_dp,  1.0_dp /)
      r_grid(1:3,6) = (/  0.0_dp,  0.0_dp, -1.0_dp /)
      w_grid(1:6)   = 1.0_dp / 6.0_dp
      if ( n_psd_grd_pnts.ge.18_int32 ) then
         l_max_grid = 5_int32
         const = sqrt(0.5_dp)
         r_grid(1:3, 7) = (/  const,  const, 0.0_dp /)
         r_grid(1:3, 8) = (/ -const, -const, 0.0_dp /)
         r_grid(1:3, 9) = (/  const, -const, 0.0_dp /)
         r_grid(1:3,10) = (/ -const,  const, 0.0_dp /)
         r_grid(1:3,11) = (/ 0.0_dp,  const,  const /)
         r_grid(1:3,12) = (/ 0.0_dp, -const, -const /)
         r_grid(1:3,13) = (/ 0.0_dp,  const, -const /)
         r_grid(1:3,14) = (/ 0.0_dp, -const,  const /)
         r_grid(1:3,15) = (/  const, 0.0_dp,  const /)
         r_grid(1:3,16) = (/ -const, 0.0_dp, -const /)
         r_grid(1:3,17) = (/  const, 0.0_dp, -const /)
         r_grid(1:3,18) = (/ -const, 0.0_dp,  const /)
         w_grid(1:6)     = 1.0_dp / 30.0_dp
         w_grid(7:18)    = 1.0_dp / 15.0_dp
      endif
      if ( n_psd_grd_pnts.ge.26_int32 ) then
         l_max_grid = 7_int32
         const = 1.0_dp / sqrt(3.0_dp)
         r_grid(1:3,19) = (/  const,  const,  const /)
         r_grid(1:3,20) = (/ -const, -const, -const /)
         r_grid(1:3,21) = (/ -const,  const,  const /)
         r_grid(1:3,22) = (/  const, -const,  const /)
         r_grid(1:3,23) = (/  const,  const, -const /)
         r_grid(1:3,24) = (/  const, -const, -const /)
         r_grid(1:3,25) = (/ -const,  const, -const /)
         r_grid(1:3,26) = (/ -const, -const,  const /)
         w_grid(1:6)     = 1.0_dp / 21.0_dp
         w_grid(7:18)    = 4.0_dp / 105.0_dp
         w_grid(19:26)   = 27.0_dp / 840.0_dp
      endif
      if ( n_psd_grd_pnts.ge.50_int32 ) then
         l_max_grid = 11_int32
         const = 1.0_dp / sqrt(11.0_dp)
         r_grid(1:3,27) = (/  3.0_dp *const,  const,  const /)
         r_grid(1:3,28) = (/ -3.0_dp *const, -const, -const /)
         r_grid(1:3,29) = (/ -3.0_dp *const,  const,  const /)
         r_grid(1:3,30) = (/  3.0_dp *const, -const,  const /)
         r_grid(1:3,31) = (/  3.0_dp *const,  const, -const /)
         r_grid(1:3,32) = (/  3.0_dp *const, -const, -const /)
         r_grid(1:3,33) = (/ -3.0_dp *const,  const, -const /)
         r_grid(1:3,34) = (/ -3.0_dp *const, -const,  const /)
         r_grid(1:3,35) = (/  const,  3.0_dp *const,  const /)
         r_grid(1:3,36) = (/ -const, -3.0_dp *const, -const /)
         r_grid(1:3,37) = (/ -const,  3.0_dp *const,  const /)
         r_grid(1:3,38) = (/  const, -3.0_dp *const,  const /)
         r_grid(1:3,39) = (/  const,  3.0_dp *const, -const /)
         r_grid(1:3,40) = (/  const, -3.0_dp *const, -const /)
         r_grid(1:3,41) = (/ -const,  3.0_dp *const, -const /)
         r_grid(1:3,42) = (/ -const, -3.0_dp *const,  const /)
         r_grid(1:3,43) = (/  const,  const,  3.0_dp *const /)
         r_grid(1:3,44) = (/ -const, -const, -3.0_dp *const /)
         r_grid(1:3,45) = (/ -const,  const,  3.0_dp *const /)
         r_grid(1:3,46) = (/  const, -const,  3.0_dp *const /)
         r_grid(1:3,47) = (/  const,  const, -3.0_dp *const /)
         r_grid(1:3,48) = (/  const, -const, -3.0_dp *const /)
         r_grid(1:3,49) = (/ -const,  const, -3.0_dp *const /)
         r_grid(1:3,50) = (/ -const, -const,  3.0_dp *const /)
         w_grid(1:6)     = 4.0_dp / 315.0_dp
         w_grid(7:18)    = 64.0_dp / 2835.0_dp
         w_grid(19:26)   = 27.0_dp / 1280.0_dp
         w_grid(27:50)   = 14641.0_dp / 725760.0_dp
      endif
   end subroutine octahedron_grid
   subroutine icosahedron_grid( )
      real(dp) :: theta1, theta2, phi
      integer(int32) :: i1
      l_max_grid = 5_int32
      r_grid(1:3,1) = (/  0.0_dp,  0.0_dp,  1.0_dp /)
      r_grid(1:3,2) = (/  0.0_dp,  0.0_dp, -1.0_dp /)
      theta1 = atan(2.0_dp)
      theta2 = pi - theta1
      do i1 = 0, 4
         phi = twopi / 5.0_dp * dble(i1)
         r_grid(1,3+i1) = sin(theta1) * cos(phi)
         r_grid(2,3+i1) = sin(theta1) * sin(phi)
         r_grid(3,3+i1) = cos(theta1)
         phi = phi + pi / 5.0_dp
         r_grid(1,8+i1) = sin(theta2) * cos(phi)
         r_grid(2,8+i1) = sin(theta2) * sin(phi)
         r_grid(3,8+i1) = cos(theta2)
      enddo
      w_grid(1:12) = 1.0_dp / 12.0_dp
      if ( n_psd_grd_pnts.ge.32_int32 ) then
         l_max_grid = 9_int32
         phi = 1.0_dp / sqrt(15.0_dp+6.0*sqrt(5.0))
         theta1 = acos((2.0_dp+sqrt(5.0)) * phi)
         theta2 = acos(phi)
         do i1 = 0, 4
            phi = ( twopi * dble(i1) + pi )/ 5.0_dp
            r_grid(1,13+i1) = sin(theta1) * cos(phi)
            r_grid(2,13+i1) = sin(theta1) * sin(phi)
            r_grid(3,13+i1) = cos(theta1)
            r_grid(1,18+i1) = sin(theta2) * cos(phi)
            r_grid(2,18+i1) = sin(theta2) * sin(phi)
            r_grid(3,18+i1) = cos(theta2)
         enddo
         theta1 = pi - theta1
         theta2 = pi - theta2
         do i1 = 0, 4
            phi = twopi / 5.0_dp * dble(i1)
            r_grid(1,23+i1) = sin(theta1) * cos(phi)
            r_grid(2,23+i1) = sin(theta1) * sin(phi)
            r_grid(3,23+i1) = cos(theta1)
            r_grid(1,28+i1) = sin(theta2) * cos(phi)
            r_grid(2,28+i1) = sin(theta2) * sin(phi)
            r_grid(3,28+i1) = cos(theta2)
         enddo
         w_grid(1:12)  =  5.0_dp / 168.0_dp
         w_grid(13:32) = 27.0_dp / 840.0_dp
      endif
   end subroutine icosahedron_grid
   subroutine rndm_nlc_grd( iw, i_at, i_fe )
      integer(int32),     intent(in) :: iw
      integer(int32),     intent(in) :: i_at
      integer(int32),     intent(in) :: i_fe
      real(dp), dimension(3,1:3) :: rotmat
      real(dp)                     :: phi,   sin_phi,   cos_phi
      real(dp)                     :: psi,   sin_psi,   cos_psi
      real(dp)                     :: theta, sin_theta, cos_theta
      integer(int32) :: i1
      if ( no_rot_grid ) then
         phi       = angles(1,iw)
         psi       = angles(2,iw)
         cos_theta = angles(3,iw)
      else
         phi   = twopi*rndn(iw)%rndc()
         psi   = twopi*rndn(iw)%rndc()
         cos_theta = 1.0_dp-2.0_dp*rndn(iw)%rndc()
      endif
      theta = acos(cos_theta)
      sin_phi   = sin(phi)
      sin_psi   = sin(psi)
      sin_theta = sin(theta)
      cos_phi   = cos(phi)
      cos_psi   = cos(psi)
      rotmat(1,1) = cos_psi * cos_phi - cos_theta * sin_phi * sin_psi
      rotmat(2,1) = cos_psi * sin_phi + cos_theta * cos_phi * sin_psi
      rotmat(3,1) = sin_theta * sin_psi
      rotmat(1,2) =-sin_psi * cos_phi - cos_theta * sin_phi * cos_psi
      rotmat(2,2) =-sin_psi * sin_phi + cos_theta * cos_phi * cos_psi
      rotmat(3,2) = cos_psi * sin_theta
      rotmat(1,3) = sin_theta * sin_phi
      rotmat(2,3) =-sin_theta * cos_phi
      rotmat(3,3) = cos_theta
      call dgemm('T','N', 3, n_psd_grd_pnts, 3, 1.0_dp, rotmat, 3,&
      & r_grid, 3, 0.0_dp, nlc_rnd_grid(:,:,iw), 3 )
      rotmat(1:3,1) = frm_cnf(iw)%d_fn(i_at,i_fe)%v / frm_cnf(iw)%d_fn(i_at,i_fe)%m
      call dgemv('T', 3, n_psd_grd_pnts, 1.0_dp, nlc_rnd_grid(:,:,iw), 3, &
      & rotmat(1:3,1), 1_int32, 0.0_dp, nlc_cos_theta(:,iw), 1_int32 )
      do i1 = 1_int32, n_psd_grd_pnts
         nlc_rnd_grid(1:3,i1,iw) = frm_cnf(iw)%d_fn(i_at,i_fe)%m * nlc_rnd_grid(1:3,i1,iw) + r_at(1:3,i_at)
      enddo
   end subroutine rndm_nlc_grd
   pure double precision function P_l( l, x )
      integer(int32), intent(in) :: l
      real(dp),       intent(in) :: x
      select case(l)
       case(0) ! S
         P_l = 1.0_dp
       case(1) ! P
         P_l = x
       case(2) ! D
         P_l = 0.5_dp*(3.0_dp*x**2 - 1.0_dp )
       case(3) ! F
         P_l = 0.5_dp*x*(5.0_dp*x**2 - 3.0_dp )
       case(4) ! G
         P_l = 0.125_dp*(35.0_dp*x**4 - 30.0_dp*x**2 + 3.0_dp )
      end select
   end function P_l
end module nl_pseudo_grid_m
