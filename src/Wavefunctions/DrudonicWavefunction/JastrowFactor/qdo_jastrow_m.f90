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
module qdo_jastrow_m
   use fortran_kinds_v,   only: dp, int32
   use openmp_mpi_m,      only: mpi_rank
   use quantum_monte_carlo_v,        only: n_wlk_max, restart
   use molecular_system_v,        only: dist_t
   use qdo_system_v,        only: n_qdo, qdos
   use drudonic_config_c, only: drd_cnf
   use cusps_functions_params_m,        only : jcspar_t, save_cusps_functions_params, read_cusps_functions_params
   use jastrow_factors_v,        only: jc1_opt, jc2_opt, jce_opt
   use cusp_functions_m,  only: cspfun, DD2_cspfun, db_cspfun
   implicit none
   type, public :: qdocsp_t
      integer(int32), pointer :: i_drd
      type(dist_t), dimension(:,:), pointer :: d_dq
      type(dist_t), dimension(:),   pointer :: d_dd, d_dq_new, d_dd_new
      real(dp) :: jst, g_1b, g_2b, g, new_1b
      real(dp), allocatable, dimension(:) :: vec_1b, vec_2b, vec_2b_new
      real(dp), allocatable, dimension(:) :: DD2_dq, DD2_dq_new, DD2_dd_new
      real(dp), allocatable, dimension(:,:) :: DD2_dd
      real(dp), allocatable, dimension(:) :: db_1b, db_2b
   contains
      procedure :: ini  => init_cusps
      procedure :: cmp  => comp_cusps
      procedure :: r1b  => variation_1body_cusps
      procedure :: r2b  => variation_2body_cusps
      procedure :: upd  => updt_cusps
      procedure :: cDD2 => DD2_cusps
      procedure :: nDD2 => new_DD2_cusps
      procedure :: uDD2 => updt_DD2_cusps
      procedure :: dpar => dp_cusps
   end type qdocsp_t
   logical,        save, public :: jdd_prs, jqd_prs, qdo_cusp_prs, qdo_cusp_opt
   integer(int32), save, public :: n_dd_jst_par, n_qd_jst_par, n_qdo_jst_par
   type(jcspar_t), save, public :: b_dd, b_qd
   type(qdocsp_t), save, public, allocatable, dimension(:) :: qdo_jastrow
   public :: init_qdocusps_par, save_jstcsp, upd_par_jstcsp
   private :: init_cusps, comp_cusps, variation_1body_cusps, variation_2body_cusps, updt_cusps, &
   & DD2_cusps, new_DD2_cusps, updt_DD2_cusps, dp_cusps
contains
   subroutine init_qdocusps_par
      integer(int32) :: i1
      n_dd_jst_par = 0_int32
      n_qd_jst_par = 0_int32
      n_qdo_jst_par = 0_int32
      qdo_cusp_prs = .false.
      qdo_cusp_opt = .false.
      jdd_prs = .false.
      jqd_prs = .false.
      if (n_qdo.gt.1_int32) then
         if ( b_qd%typ.ne.0_int32 ) jqd_prs = .true.
         if ( b_dd%typ.ne.0_int32 ) jdd_prs = .true.
      endif
      if (jqd_prs) then
         b_qd%n_cmp = 1_int32
         b_qd%n_par = 2_int32*b_qd%ord + 1_int32
         if ( jc1_opt ) then
            b_qd%n_opt_par = b_qd%ord + 1_int32
            if ( jce_opt ) b_qd%n_opt_par = b_qd%n_opt_par + b_qd%ord
         else
            b_qd%n_opt_par = 0_int32
         endif
         allocate(b_qd%b(1:b_qd%n_par,1) ) ; b_qd%b = 0.0_dp
         b_qd%b(1,:) = 1.0_dp
         if ( b_qd%ord.gt.0 ) then
            b_qd%b(b_qd%ord+2,:) = 1.0_dp
            if ( b_qd%ord.gt.1 ) then
               do i1 = 2_int32, b_qd%ord
                  b_qd%b(b_qd%ord+i1+1,:) = b_qd%b(b_qd%ord+i1,:) / 2.0_dp
               enddo
            endif
         endif
         n_qd_jst_par = b_qd%n_opt_par * b_qd%n_cmp
      endif
      if (jdd_prs) then
         b_dd%n_cmp = 1_int32
         b_dd%n_par = 2_int32*b_dd%ord + 1_int32
         if ( jc2_opt ) then
            b_dd%n_opt_par = b_dd%ord + 1_int32
            if ( jce_opt ) b_dd%n_opt_par = b_dd%n_opt_par + b_dd%ord
         else
            b_dd%n_opt_par = 0_int32
         endif
         allocate(b_dd%b(1:b_dd%n_par,1) ) ; b_dd%b = 0.0_dp
         b_dd%b(1,:) = 1.0_dp
         if ( b_dd%ord.gt.0 ) then
            b_dd%b(b_dd%ord+2,:) = 1.0_dp
            if ( b_dd%ord.gt.1 ) then
               do i1 = 2_int32, b_dd%ord
                  b_dd%b(b_dd%ord+i1+1,:) = b_dd%b(b_dd%ord+i1,:) / 2.0_dp
               enddo
            endif
         endif
         n_dd_jst_par = b_dd%n_opt_par * b_dd%n_cmp
      endif
      if ( jdd_prs.or.jqd_prs ) qdo_cusp_prs = .true.
      if (restart.and.qdo_cusp_prs) then
         if ( mpi_rank.eq.0_int32 ) open(unit=1,file='wvfn.save/jst_qdo.sav',action='read',form='formatted',status='old')
         if (jqd_prs) then
            if ( mpi_rank.eq.0_int32 ) read(1,*) ! '("# One body Cusp parameters")'
            call read_cusps_functions_params( 1, b_qd )
         endif
         if (jdd_prs) then
            if ( mpi_rank.eq.0_int32 ) read(1,*) ! '("# Two body Cusp parameters")'
            call read_cusps_functions_params( 1, b_dd )
         endif
         if (mpi_rank.eq.0_int32) close(1)
      endif
      n_qdo_jst_par = n_dd_jst_par + n_qd_jst_par
      if ( (jc1_opt.or.jc2_opt).and.qdo_cusp_prs) qdo_cusp_opt =.true.
      if( mpi_rank.eq.0_int32 ) then
         write(*,*)
         write(*,'(3X,"               ________________________________                 ")')
         write(*,'(3X,"QDO-drudon cusp function                  (jqd_prs) : ",L10)') jdd_prs
         write(*,'(3X,"Drudon-drudon cusp function               (jdd_prs) : ",L10)') jqd_prs
         write(*,'(3X,"QDO-QDO cusps functions              (qdo_cusp_prs) : ",L10)') qdo_cusp_prs
         write(*,'(3X,"QDO-QDO cusps optimization           (qdo_cusp_opt) : ",L10)') qdo_cusp_opt
         write(*,*)
         write(*,'(3X,"Total number of cusp parameters     (n_qdo_jst_par) : ",I10)') n_qdo_jst_par
      endif
   end subroutine init_qdocusps_par
   subroutine init_cusps( obj, iw )
      class(qdocsp_t), intent(inout) :: obj
      integer(int32),   intent(in) :: iw
      obj%i_drd    => drd_cnf(iw)%i_drd
      obj%d_dq     => drd_cnf(iw)%d_dq
      obj%d_dd     => drd_cnf(iw)%d_dd
      obj%d_dq_new => drd_cnf(iw)%d_dq_new
      obj%d_dd_new => drd_cnf(iw)%d_dd_new
      allocate( obj%vec_1b(1:n_qdo) ) ; obj%vec_1b = 0.0_dp
      allocate( obj%vec_2b(1:n_qdo*(n_qdo-1)/2 ) ) ; obj%vec_2b = 0.0_dp
      allocate( obj%vec_2b_new(1:n_qdo) ) ; obj%vec_2b_new = 0.0_dp
      allocate( obj%DD2_dq(1:4*n_qdo) ) ; obj%DD2_dq = 0.0_dp
      allocate( obj%DD2_dd(1:4*n_qdo,1:n_qdo) ) ; obj%DD2_dd = 0.0_dp
      allocate( obj%DD2_dd_new(1:4*n_qdo) ) ; obj%DD2_dd_new = 0.0_dp
      allocate( obj%DD2_dq_new(1:4) ) ; obj%DD2_dq_new = 0.0_dp
   end subroutine init_cusps
   subroutine comp_cusps( obj )
      class(qdocsp_t), intent(inout) :: obj
      integer(int32) :: iq1, iq2, ic
      real(dp) :: G
      obj%jst = 0.0_dp
      if (jqd_prs) then
         obj%vec_1b = 0.0_dp
         do iq1 = 1_int32, n_qdo ; do iq2 = 1_int32, n_qdo
               if (iq1 .ne. iq2 ) then
                  G = -qdos(iq1)%qdo_q * qdos(iq2)%qdo_q * qdos(iq2)%qdo_m
                  obj%vec_1b(iq1) = obj%vec_1b(iq1) + cspfun( G, 1, 0, b_qd, obj%d_dq(iq2,iq1)%m )
               endif
            enddo ; enddo
         obj%jst = obj%jst + sum(obj%vec_1b)
      endif
      if (jdd_prs) then
         obj%vec_2b = 0.0_dp
         ic = 0_int32
         do iq1 = 1_int32, n_qdo ; do iq2 = iq1 + 1_int32, n_qdo
               ic = ic + 1_int32
               G = qdos(iq1)%qdo_q * qdos(iq2)%qdo_q
               G = G * qdos(iq1)%qdo_m * qdos(iq2)%qdo_m / (qdos(iq1)%qdo_m + qdos(iq2)%qdo_m )
               obj%vec_2b(ic) = cspfun( G, 1, 0, b_dd, obj%d_dd(ic)%m )
            enddo ;enddo
         obj%jst = obj%jst + sum(obj%vec_2b)
      endif
   end subroutine comp_cusps
   subroutine variation_1body_cusps( obj, g_1b )
      class(qdocsp_t), intent(inout) :: obj
      real(dp),        intent(out)   :: g_1b
      integer(int32) :: iq1
      real(dp) :: G
      if (jqd_prs) then
         obj%new_1b = 0.0_dp
         do iq1 = 1, n_qdo
            if (iq1.ne.obj%i_drd) then
               G = -qdos(iq1)%qdo_q * qdos(obj%i_drd)%qdo_q* qdos(obj%i_drd)%qdo_m
               obj%new_1b = obj%new_1b  + cspfun( G, 1, 0, b_qd, obj%d_dq_new(iq1)%m )
            endif
         enddo
         g_1b = obj%new_1b - obj%vec_1b(obj%i_drd)
         g_1b = exp(g_1b)
      else
         g_1b = 1.0_dp
      endif
   end subroutine variation_1body_cusps
   subroutine variation_2body_cusps( obj, g_2b )
      class(qdocsp_t), intent(inout) :: obj
      real(dp),        intent(out)   :: g_2b
      integer(int32) :: iq1, iq2
      real(dp) :: G
      if (jdd_prs) then
         obj%vec_2b_new = 0.0_dp
         do iq1 = 1, n_qdo
            if (iq1.ne.obj%i_drd) then
               G = qdos(iq1)%qdo_q * qdos(obj%i_drd)%qdo_q
               G = G * qdos(iq1)%qdo_m * qdos(obj%i_drd)%qdo_m / (qdos(iq1)%qdo_m + qdos(obj%i_drd)%qdo_m )
               obj%vec_2b_new(iq1) = cspfun( G, 1, 0, b_dd, obj%d_dd_new(iq1)%m )
            endif
         enddo
         g_2b = sum(obj%vec_2b_new)
         do iq1 = 1, obj%i_drd - 1
            iq2 = n_qdo*(iq1-1) - iq1*(iq1+1)/2 + obj%i_drd
            g_2b = g_2b - obj%vec_2b(iq2)
         enddo
         iq2 = n_qdo*(obj%i_drd-1)-obj%i_drd*(obj%i_drd+1)/2 + obj%i_drd
         do iq1 = obj%i_drd + 1, n_qdo
            iq2 = iq2 + 1_int32
            g_2b = g_2b - obj%vec_2b(iq2)
         enddo
         g_2b = exp(g_2b)
      else
         g_2b = 1.0_dp
      endif
   end subroutine variation_2body_cusps
   subroutine updt_cusps( obj )
      class(qdocsp_t), intent(inout) :: obj
      integer(int32) :: iq1, iq2
      obj%jst = 0.0_dp
      if (jqd_prs) then
         obj%vec_1b(obj%i_drd)  = obj%new_1b
         obj%jst = obj%jst + sum(obj%vec_1b)
      endif
      if (jdd_prs) then
         do iq1 = 1, obj%i_drd - 1
            iq2 = n_qdo*(iq1-1) - iq1*(iq1+1)/2 + obj%i_drd
            obj%vec_2b(iq2) = obj%vec_2b_new(iq1)
         enddo
         iq2 = n_qdo*(obj%i_drd-1)-obj%i_drd*(obj%i_drd+1)/2 + obj%i_drd
         do iq1 = obj%i_drd + 1, n_qdo
            iq2 = iq2 + 1_int32
            obj%vec_2b(iq2) = obj%vec_2b_new(iq1)
         enddo
         obj%jst = obj%jst + sum(obj%vec_2b)
      endif
   end subroutine updt_cusps
   subroutine DD2_cusps( obj, DD2_vec )
      class(qdocsp_t), intent(inout) :: obj
      real(dp), dimension(4*n_qdo), intent(inout) :: DD2_vec
      integer(int32) :: iq1, iq2, ic
      real(dp) :: G
      if (jqd_prs) then
         obj%DD2_dq(:) = 0.0_dp
         do iq1 = 1_int32, n_qdo ; do iq2 = 1_int32, n_qdo
               if (iq1 .ne. iq2 ) then
                  G = -qdos(iq1)%qdo_q * qdos(iq2)%qdo_q*qdos(iq1)%qdo_m
                  obj%DD2_dq(4*iq1-3:4*iq1) = obj%DD2_dq(4*iq1-3:4*iq1) + &
                  & DD2_cspfun( G, 1, 0, b_qd, obj%d_dq(iq2,iq1)%v, obj%d_dq(iq2,iq1)%m )
               endif
            enddo ; enddo
         DD2_vec(:) = DD2_vec(:) + obj%DD2_dq(:)
      endif
      if (jdd_prs) then
         obj%DD2_dd = 0.0_dp
         ic = 0_int32
         do iq2 = 1_int32, n_qdo ; do iq1 = iq2 + 1_int32, n_qdo
               ic = ic + 1_int32
               G = qdos(iq1)%qdo_q * qdos(iq2)%qdo_q
               G = G * qdos(iq1)%qdo_m * qdos(iq2)%qdo_m / (qdos(iq1)%qdo_m + qdos(iq2)%qdo_m )
               obj%DD2_dd(4*iq2-3:4*iq2,iq1) = DD2_cspfun( G, 1, 0, b_dd, obj%d_dd(ic)%v(:), obj%d_dd(ic)%m )
               obj%DD2_dd(4*iq1-3:4*iq1-1,iq2) = - obj%DD2_dd(4*iq2-3:4*iq2-1,iq1)
               obj%DD2_dd(4*iq1,iq2)           =   obj%DD2_dd(4*iq2,iq1)
            enddo ; enddo
         do iq1 = 1_int32, n_qdo
            DD2_vec(:) = DD2_vec(:) + obj%DD2_dd(:,iq1)
         enddo
      endif
   end subroutine DD2_cusps
   subroutine new_DD2_cusps( obj, DD2_new )
      class(qdocsp_t),          intent(inout) :: obj
      real(dp), dimension(4), intent(inout) :: DD2_new
      integer(int32) :: iq1
      real(dp) :: G
      if (jqd_prs) then
         obj%DD2_dq_new = 0.0_dp
         do iq1 = 1_int32, n_qdo
            if ( iq1 .ne. obj%i_drd ) then
               G = -qdos(iq1)%qdo_q * qdos(obj%i_drd)%qdo_q* qdos(obj%i_drd)%qdo_m
               obj%DD2_dq_new(1:4) = obj%DD2_dq_new(1:4) + DD2_cspfun( G, 1, 0, b_qd, obj%d_dq_new(iq1)%v, obj%d_dq_new(iq1)%m )
            endif
         enddo 
         DD2_new(1:4) = DD2_new(1:4) +  obj%DD2_dq_new(1:4)
      endif
      if (jdd_prs) then
         obj%DD2_dd_new = 0.0_dp
         do iq1 = 1_int32, n_qdo
            if ( iq1 .ne. obj%i_drd ) then
               G = qdos(iq1)%qdo_q * qdos(obj%i_drd)%qdo_q
               G = G * qdos(iq1)%qdo_m * qdos(obj%i_drd)%qdo_m / (qdos(iq1)%qdo_m +  qdos(obj%i_drd)%qdo_m )
               obj%DD2_dd_new(4*iq1-3:4*iq1) = DD2_cspfun( G, 1, 0, b_dd, obj%d_dd_new(iq1)%v, obj%d_dd_new(iq1)%m )
               DD2_new(1:4) = DD2_new(1:4) + obj%DD2_dd_new(4*iq1-3:4*iq1)
            endif
         enddo
      endif
   end subroutine new_DD2_cusps
   subroutine updt_DD2_cusps( obj, DD2_vec )
      class(qdocsp_t), intent(inout) :: obj
      real(dp), dimension(4*n_qdo), intent(inout) :: DD2_vec
      integer(int32) :: iq1
      if (jqd_prs) then
         obj%DD2_dq(4*obj%i_drd-3:4*obj%i_drd) = obj%DD2_dq_new(1:4)
         DD2_vec(1:4*n_qdo) = DD2_vec(1:4*n_qdo) + obj%DD2_dq(1:4*n_qdo)
      endif
      if (jdd_prs) then
         do iq1 = 1_int32, n_qdo
            obj%DD2_dd(4*obj%i_drd-3:4*obj%i_drd,iq1) = obj%DD2_dd_new(4*iq1-3:4*iq1)
            obj%DD2_dd(4*iq1-3:4*iq1-1,obj%i_drd) = -obj%DD2_dd_new(4*iq1-3:4*iq1-1)
            obj%DD2_dd(4*iq1,obj%i_drd)           = obj%DD2_dd_new(4*iq1)
         enddo
         do iq1 = 1_int32, n_qdo
            DD2_vec(:) = DD2_vec(:) + obj%DD2_dd(:,iq1)
         enddo
      endif
   end subroutine updt_DD2_cusps
   subroutine dp_cusps( obj, vec_par_var )
      class(qdocsp_t), intent(inout) :: obj
      real(dp), dimension(n_qdo_jst_par), intent(inout) :: vec_par_var
      integer(int32) :: iq1, iq2, ic, ip
      real(dp) :: G
      vec_par_var = 0.0_dp
      ip = 0_int32
      if (jc1_opt.and.jqd_prs) then
         do iq1 = 1_int32, n_qdo ; do iq2 = 1_int32, n_qdo
               if (iq1 .ne. iq2 ) then
                  G = -qdos(iq1)%qdo_q * qdos(iq2)%qdo_q * qdos(iq1)%qdo_m
                  vec_par_var(ip+1:ip+n_qd_jst_par) = vec_par_var(ip+1:ip+n_qd_jst_par) + db_cspfun( G, 1, 0, b_qd, obj%d_dq(iq2,iq1)%m, jce_opt )
               endif
            enddo ; enddo
         ip = ip + n_qd_jst_par
      endif
      if (jc2_opt.and.jdd_prs) then
         ic = 0_int32
         do iq1 = 1_int32, n_qdo ; do iq2 = iq1 + 1_int32, n_qdo
               ic = ic + 1_int32
               G = qdos(iq1)%qdo_q * qdos(iq2)%qdo_q
               G = G * qdos(iq1)%qdo_m * qdos(iq2)%qdo_m / (qdos(iq1)%qdo_m + qdos(iq2)%qdo_m )
               vec_par_var(ip+1:ip+n_dd_jst_par) = vec_par_var(ip+1:ip+n_dd_jst_par) + db_cspfun( G, 1, 0,  b_dd, obj%d_dd(ic)%m, jce_opt )
            enddo ; enddo
         ip = ip + n_dd_jst_par
      endif
   end subroutine dp_cusps
   subroutine upd_par_jstcsp( vec_par_var )
      real(dp), dimension(n_qdo_jst_par), intent(in) :: vec_par_var
      integer(int32) :: ip
      ip = 0_int32
      if (jc1_opt.and.jqd_prs) then
         call update_cusp ( b_qd, n_qd_jst_par, vec_par_var(ip+1:ip+n_qd_jst_par) )
         ip = ip + n_qd_jst_par
      endif
      if (jc2_opt.and.jdd_prs) then
         call update_cusp ( b_dd, n_dd_jst_par, vec_par_var(ip+1:ip+n_dd_jst_par) )
         ip = ip + n_dd_jst_par
      endif
   contains
      subroutine update_cusp ( b, n_par, vec_var )
         type(jcspar_t), intent(inout) :: b
         integer(int32) :: n_par
         real(dp), dimension(n_par), intent(in) :: vec_var
         real(dp) :: z_tmp
         integer(int32) :: i1, i2, i3
         i3 = 0_int32
         do i1 = 1_int32, b%n_cmp
            z_tmp = b%b(1,i1) + vec_var(i3+1)
            if( z_tmp.gt.0.0_dp ) b%b(1,i1) = z_tmp
            i3 = i3 + 1_int32
            if ( b%ord.gt.0_int32 ) then
               b%b(2:b%ord+1,i1) = b%b(2:b%ord+1,i1) + vec_var(i3+1:i3+b%ord)
               i3 = i3 + b%ord
               if ( jce_opt ) then
                  do i2 = 1_int32, b%ord
                     z_tmp = b%b(b%ord+1+i2,i1) + vec_var(i3+i2)
                     if( z_tmp.gt.0.0_dp ) then
                        b%b(b%ord+1+i2,i1) = z_tmp
                     endif
                  enddo
                  i3 = i3 + b%ord
               endif
            endif
         enddo
      end subroutine update_cusp
   end subroutine upd_par_jstcsp
   subroutine save_jstcsp( )
      if ( mpi_rank.eq.0_int32 ) open(unit=1,file='wvfn.save/jst_qdo.sav',action='write',form='formatted',status='unknown')
      if (jqd_prs) then
         if ( mpi_rank.eq.0_int32 ) write(1,'("# One body Cusp parameters")')
         call save_cusps_functions_params( 1, b_qd )
      endif
      if ( jdd_prs ) then
         if ( mpi_rank.eq.0_int32 ) write(1,'("# Two body Cusp parameters")')
         call save_cusps_functions_params( 1, b_dd )
      endif
      if (mpi_rank.eq.0_int32) close(1)
   end subroutine save_jstcsp
end module qdo_jastrow_m
