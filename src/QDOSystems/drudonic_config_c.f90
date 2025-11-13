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
module drudonic_config_c
   use fortran_kinds_v,  only: int32, dp
   use mersenne_twister19937_m, only: rndn
   use qdo_system_v, only: n_qdo, r_qdo, n_pchrg, r_pchrg
   use molecular_system_v, only: n_fe, n_at, dist_t
   implicit none
   type, public :: drdcnf_t
      integer(int32)                            :: i_drd
      integer(int32)                            :: i_tmv_q
      integer(int32), allocatable, dimension(:) :: i_drd_tmv
      integer(int32)                            :: n_drd_tmv
      real(dp),     allocatable, dimension(:,:) :: r_drd
      type(dist_t), allocatable, dimension(:)   :: d_dd
      type(dist_t), allocatable, dimension(:,:) :: d_dq
      real(dp),                  dimension(3)   :: r_drd_new
      type(dist_t), allocatable, dimension(:)   :: d_dd_new
      type(dist_t), allocatable, dimension(:)   :: d_dq_new
      integer(int32)                            :: i_type
      type(dist_t), allocatable, dimension(:,:) :: d_df
      type(dist_t), allocatable, dimension(:,:) :: d_fq
      type(dist_t), allocatable, dimension(:,:) :: d_dn
      type(dist_t), allocatable, dimension(:,:) :: d_fpch
      type(dist_t), allocatable, dimension(:,:) :: d_dpch
      type(dist_t), allocatable, dimension(:)   :: d_df_new
      type(dist_t), allocatable, dimension(:)   :: d_fd_new
      type(dist_t), allocatable, dimension(:)   :: d_fq_new
      type(dist_t), allocatable, dimension(:)   :: d_dn_new
      type(dist_t), allocatable, dimension(:)   :: d_fpch_new
      type(dist_t), allocatable, dimension(:)   :: d_dpch_new
   contains
      procedure :: ini => ini_drd_cnf
      procedure :: cmp => cmp_drd_cnf
      procedure :: chs => chs_drd_cnf
      procedure :: new => new_drd_cnf
      procedure :: upd => upd_drd_cnf
   end type drdcnf_t
   type(drdcnf_t), save, public, allocatable, dimension(:), target :: drd_cnf
   private :: ini_drd_cnf, cmp_drd_cnf, chs_drd_cnf, new_drd_cnf, upd_drd_cnf
contains
   subroutine ini_drd_cnf( obj, iw )
      class(drdcnf_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      integer(int32) :: i1
      allocate( obj%r_drd(1:3,1:n_qdo) )
      allocate( obj%d_dd(1:n_qdo*(n_qdo-1)/2) )
      allocate( obj%d_dq(1:n_qdo,1:n_qdo) )
      allocate( obj%d_dd_new(1:n_qdo) )
      allocate( obj%d_dq_new(1:n_qdo) )
      allocate( obj%i_drd_tmv(1:n_qdo) )
      do i1 = 1_int32, n_qdo
         obj%r_drd(1,i1) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_qdo(1,i1)
         obj%r_drd(2,i1) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_qdo(2,i1)
         obj%r_drd(3,i1) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_qdo(3,i1)
      enddo  ! i1 (n_qdo)
      if(n_at.gt.0_int32) then
         allocate( obj%d_dn( 1:n_qdo, 1:n_at  ) )
         allocate( obj%d_dn_new( 1:n_at  ) )
      endif
      if (n_fe.gt.0_int32) then
         allocate( obj%d_df_new( 1:n_fe  ) )
         allocate( obj%d_fd_new(1:n_qdo) )
         allocate( obj%d_fq_new(1:n_qdo) )
         allocate( obj%d_df(1:n_qdo,1:n_fe) )
         allocate( obj%d_fq(1:n_fe,1:n_qdo) )
      endif
      if (n_pchrg.gt.0_int32) then
         allocate( obj%d_dpch( 1:n_qdo, 1:n_pchrg ) )
         allocate( obj%d_dpch_new( 1:n_pchrg ) )
      endif
      if ( (n_fe.gt.0_int32).and.(n_pchrg.gt.0_int32) ) then
         allocate( obj%d_fpch( 1:n_fe, 1:n_pchrg ) )
         allocate( obj%d_fpch_new( 1:n_pchrg ) )
      endif
      obj%n_drd_tmv = n_qdo
      do i1 = 1_int32, n_qdo
         obj%i_drd_tmv(i1) = i1
      enddo ! i1 ( n_qdo )
      call obj%cmp( iw )
   end subroutine ini_drd_cnf
   subroutine cmp_drd_cnf( obj, iw )
      use qdo_system_v, only: n_qdo, r_qdo, n_pchrg, r_pchrg
      use molecular_system_v, only: n_fe, n_at, r_at
      use fermionic_config_c, only: frm_cnf
      class(drdcnf_t), intent(inout) :: obj
      integer(int32),  intent(in)    :: iw
      integer(int32) :: i1, i2, i3
      if( n_qdo.ge.2_int32 ) then
         i3 = 0_int32
         do i2 = 1_int32, n_qdo
            do i1 = i2 + 1_int32, n_qdo
               i3 = i3 + 1_int32
               obj%d_dd(i3)%v(:) = obj%r_drd(:,i2) - obj%r_drd(:,i1)
               obj%d_dd(i3)%m    = dsqrt(sum(obj%d_dd(i3)%v(:)**2))
            enddo ; enddo 
      endif
      do i1 = 1_int32, n_qdo ; do i2 = 1_int32, n_qdo
            obj%d_dq(i2,i1)%v(:) = obj%r_drd(:,i1) - r_qdo(:,i2)
            obj%d_dq(i2,i1)%m    = dsqrt(sum(obj%d_dq(i2,i1)%v(:)**2))
         enddo ; enddo 
      if (n_at.gt.0_int32) then
         do i1 = 1_int32, n_qdo
            do i2 = 1_int32, n_at
               obj%d_dn(i1,i2)%v(:) = obj%r_drd(:,i1) - r_at(:,i2)
               obj%d_dn(i1,i2)%m    = dsqrt( sum( obj%d_dn( i1, i2 )%v(:)**2 ) )
            enddo 
         enddo
      endif
      if (n_fe.gt.0_int32) then
         do i1 = 1_int32, n_qdo
            do i2 = 1_int32, n_fe
               obj%d_df(i1,i2)%v(:) = obj%r_drd(:,i1) - frm_cnf(iw)%r_fe(:,i2)
               obj%d_df(i1,i2)%m    = dsqrt( sum( obj%d_df( i1, i2 )%v(:)**2 ) )
               obj%d_fq(i2,i1)%v(:) = frm_cnf(iw)%r_fe(:,i2) - r_qdo(:,i1)
               obj%d_fq(i2,i1)%m    = dsqrt( sum( obj%d_fq( i2, i1 )%v(:)**2 ) )
            enddo
         enddo 
      endif
      if( n_pchrg.gt.0_int32 ) then
         do i1 = 1_int32, n_qdo
            do i2 = 1_int32, n_pchrg
               obj%d_dpch(i1,i2)%v(:)  = obj%r_drd(:,i1) - r_pchrg(:,i2)
               obj%d_dpch(i1,i2)%m     = dsqrt( sum( obj%d_dpch( i1, i2 )%v(:)**2 ) )
            enddo
         enddo
      endif
      if ( ( n_pchrg.gt.0_int32 ).and.(n_fe.gt.0_int32) ) then
         do i1 = 1_int32, n_fe
            do i2 = 1_int32, n_pchrg
               obj%d_fpch(i1,i2)%v(:)  = frm_cnf(iw)%r_fe(:,i1) - r_pchrg(:,i2)
               obj%d_fpch(i1,i2)%m     = dsqrt( sum( obj%d_fpch( i1, i2 )%v(:)**2 ) )
            enddo
         enddo
      endif
   end subroutine cmp_drd_cnf
   subroutine chs_drd_cnf( obj, iw )
      class(drdcnf_t), intent(inout) :: obj
      integer(int32), intent(in) :: iw
      if( obj%n_drd_tmv.eq.1_int32) then
         obj%i_tmv_q = 1_int32
      else
         obj%i_tmv_q = int( dble(obj%n_drd_tmv) * rndn(iw)%rndo() + 1.0_dp )
      endif
      obj%i_drd = obj%i_drd_tmv(obj%i_tmv_q)
      if( obj%n_drd_tmv.eq.1_int32) then
         obj%n_drd_tmv = n_qdo
      else
         if(obj%i_tmv_q.ne.obj%n_drd_tmv) then
            obj%i_drd_tmv(obj%i_tmv_q)      = obj%i_drd_tmv(obj%n_drd_tmv)
            obj%i_drd_tmv(obj%n_drd_tmv)   = obj%i_drd
         endif
         obj%n_drd_tmv = obj%n_drd_tmv - 1_int32
      endif
   end subroutine chs_drd_cnf
   subroutine new_drd_cnf( obj, iw, r_new )
      use qdo_system_v, only: n_qdo, r_qdo, n_pchrg, r_pchrg
      use molecular_system_v, only: n_fe, n_at, r_at
      use fermionic_config_c, only: frm_cnf
      class(drdcnf_t), intent(inout) :: obj
      integer(int32),           intent(in)    :: iw
      real(dp),       optional, intent(in)    :: r_new(1:3)
      integer(int32) :: i1
      if ( present(r_new) ) then
         obj%i_type = 1_int32
      else
         obj%i_type = 0_int32
      endif
      if ( obj%i_type.eq.1_int32 ) then
         obj%r_drd_new(:) = r_new(:)
         do i1 = 1_int32, n_qdo
            obj%d_dq_new(i1)%v(:) = obj%r_drd_new(:) - r_qdo(:,i1)
            obj%d_dq_new(i1)%m    = dsqrt(sum(obj%d_dq_new(i1)%v(:)**2))
         enddo
         if( n_qdo.ge.2_int32 ) then
            do i1 = 1_int32, n_qdo
               if( obj%i_drd.ne.i1 ) then
                  obj%d_dd_new(i1)%v(:) = obj%r_drd_new(:) - obj%r_drd(:,i1)
                  obj%d_dd_new(i1)%m    = dsqrt(sum(obj%d_dd_new(i1)%v(:)**2))
               else
                  obj%d_dd_new(i1)%v(:) = 0.0_dp
                  obj%d_dd_new(i1)%m    = 0.0_dp
               endif
            enddo
         endif
         if (n_fe.gt.0_int32) then
            do i1 = 1_int32, n_fe
               obj%d_df_new( i1 )%v(:) = obj%r_drd_new(:) - frm_cnf(iw)%r_fe(:,i1)
               obj%d_df_new( i1 )%m    = dsqrt( sum( obj%d_df_new( i1 )%v(:)**2 ) )
            enddo 
         endif
         if (n_at.gt.0_int32) then
            do i1 = 1_int32, n_at
               obj%d_dn_new( i1 )%v(:) = obj%r_drd_new(:) - r_at(:,i1)
               obj%d_dn_new( i1 )%m    = dsqrt( sum( obj%d_dn_new( i1 )%v(:)**2 ) )
            enddo
         endif
         if (n_pchrg.gt.0_int32) then
            do i1 = 1_int32, n_pchrg
               obj%d_dpch_new( i1 )%v(:) = obj%r_drd_new(:) - r_pchrg(:,i1)
               obj%d_dpch_new( i1 )%m    = dsqrt( sum( obj%d_dpch_new( i1 )%v(:)**2 ) )
            enddo
         endif
         do i1 = 1_int32, n_qdo
            obj%d_fd_new( i1 )%v(:) = obj%r_drd(:,i1) - frm_cnf(iw)%r_fe_new(:)
            obj%d_fd_new( i1 )%m    = dsqrt( sum( obj%d_fd_new( i1 )%v(:)**2 ) )
            obj%d_fq_new( i1 )%v(:) = frm_cnf(iw)%r_fe_new(:) - r_qdo(:,i1)
            obj%d_fq_new( i1 )%m    = dsqrt( sum( obj%d_fq_new( i1 )%v(:)**2 ) )
         enddo
         if (n_pchrg.gt.0_int32) then
            do i1 = 1_int32, n_pchrg
               obj%d_fpch_new( i1 )%v(:) = frm_cnf(iw)%r_fe_new(:) - r_pchrg(:,i1)
               obj%d_fpch_new( i1 )%m    = dsqrt( sum( obj%d_fpch_new( i1 )%v(:)**2 ) )
            enddo
         endif
      endif
   end subroutine new_drd_cnf
   subroutine upd_drd_cnf( obj, iw )
      use qdo_system_v, only: n_qdo, n_pchrg
      use molecular_system_v, only: n_fe, n_at
      use fermionic_config_c, only: frm_cnf
      class(drdcnf_t),          intent(inout) :: obj
      integer(int32), optional, intent(in)    :: iw
      integer(int32) :: i1, i2
      if ( obj%i_type.eq.1_int32 ) then
         obj%r_drd(:,obj%i_drd) = obj%r_drd_new(:)
         obj%d_dq(:,obj%i_drd) = obj%d_dq_new(:)
         if( n_qdo.ge.2_int32 ) then
            do i1 = 1, obj%i_drd-1
               i2 = n_qdo * (i1-1)-i1 * (i1+1) / 2 + obj%i_drd
               obj%d_dd(i2)%v(:) = - obj%d_dd_new(i1)%v(:)
               obj%d_dd(i2)%m    =   obj%d_dd_new(i1)%m
            enddo
            i2 = n_qdo * (obj%i_drd-1)-obj%i_drd* (obj%i_drd+1) / 2 + obj%i_drd
            do i1 = obj%i_drd+1, n_qdo
               i2 = i2 + 1_int32
               obj%d_dd(i2)%v(:) = obj%d_dd_new(i1)%v(:)
               obj%d_dd(i2)%m    = obj%d_dd_new(i1)%m
            enddo
         endif
         if (n_fe.gt.0_int32)    obj%d_df(obj%i_drd,:)   = obj%d_df_new(:)
         if (n_at.gt.0_int32)    obj%d_dn(obj%i_drd,:)   = obj%d_dn_new(:)
         if (n_pchrg.gt.0_int32) obj%d_dpch(obj%i_drd,:) = obj%d_dpch_new(:)
      else ! obj%i_type.eq.1_int32 moved a fermion
         obj%d_df(:,frm_cnf(iw)%i_fe) = obj%d_fd_new( : )
         obj%d_fq(frm_cnf(iw)%i_fe,:) = obj%d_fq_new( : )
         if (n_pchrg.gt.0_int32) obj%d_fpch(frm_cnf(iw)%i_fe,:) = obj%d_fpch_new(:)
      endif
   end subroutine upd_drd_cnf
end module drudonic_config_c
