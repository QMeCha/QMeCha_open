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
module pfaffian_c
   use fortran_kinds_v,          only: dp, int32
   use molecular_system_m,       only: n_at
   use fermionic_orbitals_m,     only: aoc_opt, n_forbs, n_par_forbs, forbs
   use quantum_monte_carlo_v,    only: grad_updat, n_wlk_max
   use fermionic_wavefunction_v, only: spin
   use pfaffian_params_m,        only: pffcfs_t
   use pfaffian_operats_m,       only: cmp_L_mat, cmp_P_mat, upd_L_mat, upd_P_inv, & 
   & cmp_P_vec, cmp_D_mat, cmp_DD2_vec, upd_D_mat, cmp_DD2_vec, cmp_H_mat, cmp_Lp_mat
   implicit none
   type, public :: pfffun_t
      integer(int32)                                :: d_P
      real(dp),       allocatable, dimension(:)     :: pff_P, g
      integer(int32), allocatable, dimension(:)     :: n_mat_chngs
      real(dp),       allocatable, dimension(:,:,:) :: mat_P
      real(dp),       allocatable, dimension(:,:,:) :: inv_P
      real(dp),       allocatable, dimension(:,:,:) :: mat_L_u
      real(dp),       allocatable, dimension(:,:,:) :: mat_L_d
      real(dp),       allocatable, dimension(:,:)   :: vec_P
      real(dp),       allocatable, dimension(:,:)   :: vec_Pi
      real(dp),       allocatable, dimension(:,:)   :: vec_Pip
      real(dp),       allocatable, dimension(:,:,:) :: mat_D
      real(dp),       allocatable, dimension(:,:)   :: DD2_ln_pff_P
      real(dp),       allocatable, dimension(:,:)   :: DD2_ln_pff_P_new
      real(dp),       allocatable, dimension(:,:)   :: vec_LPi
      real(dp),       allocatable, dimension(:,:)   :: vec_LPip
      real(dp),       allocatable, dimension(:,:,:) :: mat_OPi
      real(dp),       allocatable, dimension(:,:,:) :: mat_Lp
      real(dp),       allocatable, dimension(:,:)   :: vec_lu, vec_ld
      real(dp),       allocatable, dimension(:,:,:) :: mat_H
      real(dp),       allocatable, dimension(:,:)   :: da_ln_pff_P
      real(dp),       allocatable, dimension(:,:)   :: work
      integer(int32), allocatable, dimension(:,:)   :: ipiv
   contains
      procedure :: ini    => ini_pff_fun
      procedure :: cmp    => cmp_pff_fun
      procedure :: ratio  => ratio_pff_fun
      procedure :: upd    => upd_pff_fun
      procedure :: cmp_D  => cmp_DD2_pff_fun
      procedure :: upd_D  => upd_DD2_pff_fun
      procedure :: new_D  => new_DD2_pff_fun
      procedure :: cmp_dl => cmp_dl_pff_fun
      procedure :: cmp_da => cmp_da_pff_fun
      procedure :: cpy_dl => cpy_dl_pff_fun
   end type pfffun_t
   private :: ini_pff_fun, cmp_pff_fun, ratio_pff_fun, upd_pff_fun, cmp_DD2_pff_fun, &
   & upd_DD2_pff_fun, cmp_dl_pff_fun, new_DD2_pff_fun, cmp_da_pff_fun, &
   & cpy_dl_pff_fun
contains
   subroutine ini_pff_fun( obj, n_fe, n_fe_u, n_fe_d, dl_opt )
      class(pfffun_t), intent(inout) :: obj
      logical,         intent(in)    :: dl_opt
      integer(int32),  intent(in)    :: n_fe
      integer(int32),  intent(in)    :: n_fe_u
      integer(int32),  intent(in)    :: n_fe_d
      if( mod(n_fe,2_int32).eq.0_int32 ) then
         obj%d_P = n_fe
      else
         obj%d_P = n_fe + 1_int32
      endif
      allocate( obj%pff_P(1:n_wlk_max) ) ; obj%pff_P = 1.0_dp
      allocate( obj%g(1:n_wlk_max) ) ; obj%g = 1.0_dp
      allocate( obj%n_mat_chngs(1:n_wlk_max) ) ; obj%n_mat_chngs = 0_int32
      allocate( obj%mat_L_u(1:n_forbs,1:obj%d_P,1:n_wlk_max) ) ; obj%mat_L_u = 0.0_dp
      allocate( obj%mat_L_d(1:n_forbs,1:obj%d_P,1:n_wlk_max) ) ; obj%mat_L_d = 0.0_dp
      allocate( obj%mat_P(1:obj%d_P,1:obj%d_P,1:n_wlk_max) )   ; obj%mat_P = 0.0_dp
      allocate( obj%inv_P(1:obj%d_P,1:obj%d_P,1:n_wlk_max) )   ; obj%inv_P = 0.0_dp
      allocate( obj%vec_P(1:obj%d_P,1:n_wlk_max) )             ; obj%vec_P = 0.0_dp
      allocate( obj%vec_Pi(1:obj%d_P,1:n_wlk_max) )            ; obj%vec_Pi = 0.0_dp
      allocate( obj%vec_Pip(1:obj%d_P,1:n_wlk_max) )           ; obj%vec_Pip = 0.0_dp
      allocate( obj%mat_D(1:n_forbs,1:n_fe,1:n_wlk_max) )      ; obj%mat_D = 0.0_dp
      allocate( obj%vec_LPi(1:2*n_forbs,1:n_wlk_max) )         ; obj%vec_LPi = 0.0_dp
      allocate( obj%vec_LPip(1:2*n_forbs,1:n_wlk_max) )        ; obj%vec_LPip = 0.0_dp
      allocate ( obj%DD2_ln_pff_P(1:4*n_fe,1:n_wlk_max) )      ; obj%DD2_ln_pff_P = 0.0_dp
      if ( grad_updat ) then
         allocate( obj%DD2_ln_pff_P_new(1:4,1:n_wlk_max) )      ; obj%DD2_ln_pff_P_new = 0.0_dp
      endif
      if ( dl_opt ) then
         allocate( obj%mat_OPi(1:n_forbs,1:n_fe_u,1:n_wlk_max) ) ; obj%mat_OPi = 0.0_dp
         allocate( obj%mat_Lp(1:n_forbs,1:n_forbs,1:n_wlk_max) ) ; obj%mat_Lp = 0.0_dp
         if(n_fe.ne.obj%d_P) then
            allocate( obj%vec_Lu(1:n_forbs,1:n_wlk_max) ) ; obj%vec_Lu = 0.0_dp
            allocate( obj%vec_Ld(1:n_forbs,1:n_wlk_max) ) ; obj%vec_Ld = 0.0_dp
         endif
         if ( (n_fe_d.gt.1_int32) .and. (spin.eq.'U') ) then
            allocate( obj%mat_H(1:n_forbs,1:2*n_forbs,1:n_wlk_max) ) ; obj%mat_H = 0.0_dp
         else
            allocate( obj%mat_H(1:n_forbs,1:n_forbs,1:n_wlk_max) ) ; obj%mat_H = 0.0_dp
         endif
      endif
      if ( aoc_opt ) then
         allocate( obj%da_ln_pff_P(1:n_par_forbs,1:n_wlk_max) ) ; obj%da_ln_pff_P = 0.0_dp
      endif
      allocate( obj%work(1:2*obj%d_P,1:n_wlk_max) )
      allocate( obj%ipiv(1:obj%d_P,1:n_wlk_max) )
   end subroutine ini_pff_fun
   subroutine cmp_pff_fun( obj, iw, n_fe, n_fe_u, n_fe_d, pff_c, mat_o )
      class(pfffun_t),                   intent(inout) :: obj
      integer(int32),                    intent(in) :: iw
      integer(int32),                    intent(in) :: n_fe
      integer(int32),                    intent(in) :: n_fe_u
      integer(int32),                    intent(in) :: n_fe_d
      type(pffcfs_t),                    intent(in) :: pff_c
      real(dp), dimension(n_forbs,n_fe), intent(in) :: mat_o
      call cmp_L_mat ( pff_c, n_fe, n_fe_u, n_fe_d, obj%d_P, mat_o, obj%mat_L_u(:,:,iw), obj%mat_L_d(:,:,iw) )
      call cmp_P_mat ( n_fe, n_fe_u, n_fe_d, obj%d_P, mat_o, obj%mat_L_u(:,:,iw), obj%mat_L_d(:,:,iw), obj%mat_P(:,:,iw) )
      call hh_trn_pff_inv( obj%d_P, obj%work(:,iw), obj%ipiv(:,iw), obj%mat_P(:,:,iw), obj%pff_P(iw), obj%inv_P(:,:,iw) )
      obj%n_mat_chngs(iw) = 0_int32
   end subroutine cmp_pff_fun
   subroutine ratio_pff_fun( obj, iw, i_fe, n_fe_u, vec_o_new, g )
      class(pfffun_t), intent(inout) :: obj
      integer(int32),               intent(in) :: iw
      integer(int32),               intent(in)    :: i_fe
      integer(int32),               intent(in)    :: n_fe_u
      real(dp), dimension(n_forbs), intent(in)    :: vec_o_new
      real(dp),                     intent(inout) :: g
      if ( i_fe.le.n_fe_u ) then
         call cmp_P_vec( i_fe, obj%d_P, vec_o_new, obj%mat_L_u(:,:,iw), obj%vec_P(:,iw) )
      else
         call cmp_P_vec( i_fe, obj%d_P, vec_o_new, obj%mat_L_d(:,:,iw), obj%vec_P(:,iw) )
      endif
      obj%vec_Pi(:,iw) = - obj%inv_P(:,i_fe,iw)
      obj%g(iw) = dot_product( obj%vec_Pi(:,iw), obj%vec_P(:,iw) )
      g = obj%g(iw)
   end subroutine ratio_pff_fun
   subroutine upd_pff_fun( obj, iw, i_fe, n_fe, n_fe_u, n_fe_d, pff_c, mat_o )
      class(pfffun_t), intent(inout) :: obj
      integer(int32),                    intent(in) :: iw
      integer(int32),                    intent(in) :: i_fe
      integer(int32),                    intent(in) :: n_fe
      integer(int32),                    intent(in) :: n_fe_u
      integer(int32),                    intent(in) :: n_fe_d
      type(pffcfs_t),                    intent(in) :: pff_c
      real(dp), dimension(n_forbs,n_fe), intent(in) :: mat_o
      call upd_L_mat( pff_c, i_fe, n_fe, n_fe_u, n_fe_d, obj%d_P, mat_o, obj%vec_LPip(:,iw),&
      & obj%mat_L_u(:,:,iw), obj%mat_L_d(:,:,iw) )
      obj%mat_P(:,i_fe,iw) = obj%vec_P(:,iw)
      obj%mat_P(i_fe,:,iw) = - obj%vec_P(:,iw)
      obj%pff_P(iw) = obj%pff_P(iw) * obj%g(iw)
      if( obj%n_mat_chngs(iw).eq.n_fe .or. abs(obj%g(iw)).lt.10d-14) then
         call hh_trn_pff_inv( obj%d_P, obj%work(:,iw), obj%ipiv(:,iw), obj%mat_P(:,:,iw), obj%pff_P(iw), obj%inv_P(:,:,iw) )
         obj%n_mat_chngs(iw) = 0_int32
      else
         call upd_P_inv( i_fe, obj%d_P, obj%g(iw), obj%vec_P(:,iw), obj%vec_Pi(:,iw), &
         & obj%vec_Pip(:,iw), obj%inv_P(:,:,iw) )
         obj%n_mat_chngs(iw) = obj%n_mat_chngs(iw) + 1_int32
      endif
   end subroutine upd_pff_fun
   subroutine cmp_DD2_pff_fun( obj, iw, n_fe, n_fe_u, n_fe_d, DD2_o )
      class(pfffun_t), intent(inout) :: obj
      integer(int32),                      intent(in) :: iw
      integer(int32),                      intent(in) :: n_fe
      integer(int32),                      intent(in) :: n_fe_u
      integer(int32),                      intent(in) :: n_fe_d
      real(dp), dimension(n_forbs,4*n_fe), intent(in) :: DD2_o
      call cmp_D_mat( n_fe, n_fe_u, n_fe_d, obj%d_P, obj%mat_L_u(:,:,iw), &
      & obj%mat_L_d(:,:,iw), obj%inv_P(:,:,iw), obj%mat_D(:,:,iw) )
      call cmp_DD2_vec( n_fe, obj%mat_D(:,:,iw), DD2_o, obj%DD2_ln_pff_P(:,iw) )
   end subroutine cmp_DD2_pff_fun
   subroutine new_DD2_pff_fun( obj, iw, i_fe, DD2_o_new )
      class(pfffun_t), intent(inout) :: obj
      integer(int32),                 intent(in) :: iw
      integer(int32),                 intent(in) :: i_fe
      real(dp), dimension(n_forbs,4), intent(in) :: DD2_o_new
      call dgemv('T', n_forbs, 4, 1.0_dp/obj%g(iw), DD2_o_new(:,1:4), n_forbs, &
      & obj%mat_D(:,i_fe,iw), 1_int32, 0.0_dp, obj%DD2_ln_pff_P_new(:,iw), 1_int32 )
      obj%DD2_ln_pff_P_new(4,iw) = obj%DD2_ln_pff_P_new(4,iw) - sum(obj%DD2_ln_pff_P_new(1:3,iw)**2)
   end subroutine new_DD2_pff_fun
   subroutine upd_DD2_pff_fun( obj, iw, n_fe, n_fe_u, n_fe_d, DD2_o )
      class(pfffun_t), intent(inout) :: obj      
      integer(int32),                      intent(in) :: iw
      integer(int32),                      intent(in) :: n_fe
      integer(int32),                      intent(in) :: n_fe_u
      integer(int32),                      intent(in) :: n_fe_d
      real(dp), dimension(n_forbs,4*n_fe), intent(in) :: DD2_o
      if( obj%n_mat_chngs(iw).eq.0_int32 ) then
         call cmp_D_mat( n_fe, n_fe_u, n_fe_d, obj%d_P, obj%mat_L_u(:,:,iw), &
         & obj%mat_L_d(:,:,iw), obj%inv_P(:,:,iw), obj%mat_D(:,:,iw) )
      else
         call upd_D_mat( n_fe, n_fe_u, n_fe_d, obj%d_P, obj%mat_L_u(:,:,iw), obj%mat_L_d(:,:,iw), obj%g(iw), &
         & obj%vec_Pi(:,iw), obj%vec_Pip(:,iw), obj%vec_LPi(:,iw), obj%vec_LPip(:,iw), obj%mat_D(:,:,iw) )
      endif
      call cmp_DD2_vec( n_fe, obj%mat_D(:,:,iw), DD2_o, obj%DD2_ln_pff_P(:,iw)  )
   end subroutine upd_DD2_pff_fun
   subroutine cmp_dl_pff_fun( obj, iw, n_fe, n_fe_u, n_fe_d, pff_c, mat_o, n_par_det, dp_ln_pff_P )
      class(pfffun_t), intent(inout) :: obj      
      integer(int32),                    intent(in)    :: iw
      integer(int32),                    intent(in)    :: n_fe
      integer(int32),                    intent(in)    :: n_fe_u
      integer(int32),                    intent(in)    :: n_fe_d
      type(pffcfs_t),                    intent(in)    :: pff_c
      real(dp), dimension(n_forbs,n_fe), intent(in)    :: mat_o
      integer(int32),                    intent(in)    :: n_par_det
      real(dp), dimension(n_par_det),    intent(inout) :: dp_ln_pff_P
      call cmp_H_mat( n_fe, n_fe_u, n_fe_d, obj%d_P, mat_o, obj%inv_P(:,:,iw), obj%mat_OPi(:,:,iw), obj%mat_H(:,:,iw) )
      call cmp_Lp_mat( n_fe, n_fe_u, n_fe_d, obj%d_P, mat_o, obj%inv_P(:,:,iw), obj%mat_OPi(:,:,iw), obj%mat_Lp(:,:,iw), obj%vec_lu(:,iw), obj%vec_ld(:,iw) )
      call obj%cpy_dl( iw, n_fe, n_fe_d, pff_c, n_par_det, dp_ln_pff_P )
   end subroutine cmp_dl_pff_fun
   subroutine cpy_dl_pff_fun( obj, iw, n_fe, n_fe_d, pff_c, n_par_det, dp_ln_pff_P )
      class(pfffun_t), intent(inout) :: obj      
      integer(int32),                 intent(in)    :: iw
      integer(int32),                 intent(in)    :: n_fe
      integer(int32),                 intent(in)    :: n_fe_d
      type(pffcfs_t),                 intent(in)    :: pff_c
      integer(int32),                 intent(in)    :: n_par_det
      real(dp), dimension(n_par_det), intent(inout) :: dp_ln_pff_P
      integer(int32) :: i1, i2
      i2 = 1_int32
      do i1 = 1_int32, pff_c%n_nz_Z
         dp_ln_pff_P(i2) = obj%mat_H(pff_c%nz_Z(i1)%o1,pff_c%nz_Z(i1)%o2,iw)
         i2 = i2 + 1_int32
      enddo
      if ( (n_fe_d.gt.1_int32) .and. (spin.eq.'U') ) then
         do i1 = 1_int32, pff_c%n_nz_Z
            dp_ln_pff_P(i2) = obj%mat_H(pff_c%nz_Z(i1)%o1,pff_c%nz_Z(i1)%o2+n_forbs,iw)
            i2 = i2 + 1_int32
         enddo
      endif
      do i1 = 1_int32, pff_c%n_nz_Lmbd
         dp_ln_pff_P(i2) = obj%mat_Lp(pff_c%nz_Lmbd(i1)%o1,pff_c%nz_Lmbd(i1)%o2,iw)
         i2 = i2 + 1_int32
      enddo
      if(obj%d_P.gt.n_fe) then
         if (spin.eq.'U') then
            do i1 = 1_int32, pff_c%n_nz_u
               dp_ln_pff_P(i2) = obj%vec_Lu(pff_c%nz_u(i1)%o1,iw)
               i2 = i2 + 1_int32
            enddo
            do i1 = 1_int32, pff_c%n_nz_d
               dp_ln_pff_P(i2) = obj%vec_Ld(pff_c%nz_d(i1)%o1,iw)
               i2 = i2 + 1_int32
            enddo
         else
            do i1 = 1_int32, pff_c%n_nz_u
               dp_ln_pff_P(i2) = obj%vec_Lu(pff_c%nz_u(i1)%o1,iw) + obj%vec_Ld(pff_c%nz_u(i1)%o1,iw)
               i2 = i2 + 1_int32
            enddo
         endif
      endif
   end subroutine cpy_dl_pff_fun
   subroutine cmp_da_pff_fun( obj, iw, n_fe, mat_da_o, da_ln_pff_P )
      class(pfffun_t), intent(inout) :: obj      
      integer(int32),                        intent(in)    :: iw
      integer(int32),                        intent(in)    :: n_fe
      real(dp), dimension(n_par_forbs,n_fe), intent(in)    :: mat_da_o
      real(dp), dimension(n_par_forbs),      intent(inout) :: da_ln_pff_P
      integer(int32) :: i1, i2, i3
      integer(int32) :: io, pi, pf
      obj%da_ln_pff_P(:,iw) = 0.0_dp
      do i1 = 1_int32, n_fe ; do i2 = 1_int32, n_at ; do i3 = 1_int32, forbs(i2)%n_orbs
         io = forbs(i2)%orb(i3)%i_orb
         if( forbs(i2)%orb(i3)%n_par.gt.0_int32 ) then
            pi = forbs(i2)%orb(i3)%i_par
            pf = pi + forbs(i2)%orb(i3)%n_par - 1_int32
            obj%da_ln_pff_P(pi:pf,iw) = obj%da_ln_pff_P(pi:pf,iw) + obj%mat_D(io,i1,iw) * mat_da_o(pi:pf,i1)
         endif
      enddo ; enddo ; enddo 
      da_ln_pff_P = da_ln_pff_P + obj%da_ln_pff_P(:,iw)
   end subroutine cmp_da_pff_fun
end module pfaffian_c
