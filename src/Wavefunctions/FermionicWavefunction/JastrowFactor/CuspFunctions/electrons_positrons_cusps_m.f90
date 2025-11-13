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
module electrons_positrons_cusps_m
   use fortran_kinds_v, only: dp, int32
   use openmp_mpi_m, only: mpi_rank
   use cusps_functions_params_m, only: jcspar_t
   use molecular_system_v, only: dist_t, n_at, atoms, f_spin, mu_ep
   use cusp_functions_m,   only: cspfun, DD2_cspfun, db_cspfun
   implicit none
   public :: init_electrons_positrons_cusps, comp_electrons_positrons_cusps, variation_electrons_positrons_cusps, &
   & DD2_electrons_positrons_cusps, new_DD2_electrons_positrons_cusps, updt_DD2_electrons_positrons_cusps, &
   & db_electrons_positrons_cusps, updt_electrons_positrons_cusps_params
contains
   subroutine init_electrons_positrons_cusps( n_el, n_po, b_ep, jcep_prs, n_par_opt, jc_opt, jce_opt)
      integer(int32), intent(in)    :: n_el, n_po
      type(jcspar_t), intent(inout) :: b_ep
      logical,        intent(inout) :: jcep_prs
      integer(int32), intent(inout) :: n_par_opt
      logical,        intent(in)    :: jc_opt, jce_opt
      integer(int32)                :: i1, i2
      if ( n_el.eq.0_int32 .or. n_po.eq.0_int32 ) b_ep%typ = 0_int32
      if ( b_ep%typ.eq.0_int32 ) then
         b_ep%ord       = 0_int32
         b_ep%n_par     = 0_int32
         b_ep%n_opt_par = 0_int32
         jcep_prs       = .false.
      else 
         jcep_prs       = .true.
         if ( b_ep%typ.eq.3 ) then
            b_ep%n_par = b_ep%ord + 1_int32
         else if ( b_ep%typ.eq.4 ) then
            b_ep%n_par = 2_int32*b_ep%ord + 2_int32
         else
            b_ep%n_par = 2_int32*b_ep%ord + 1_int32
         endif
         b_ep%n_cmp = 1_int32
         allocate( b_ep%b(1:b_ep%n_par,b_ep%n_cmp )) ; b_ep%b = 0.0_dp
         if ( b_ep%typ.eq.4 ) then
            b_ep%b(1,:) = 0.0_dp
            b_ep%b(2,:) = 0.0_dp
            i2 = 2_int32
         else
            b_ep%b(1,:) = 1.0_dp
            i2 = 1_int32
         endif
         if ( b_ep%ord.gt.0 .and. b_ep%typ .ne. 3) then
            b_ep%b(b_ep%ord+i2+1,:) = 0.5_dp
            do i1 = 2_int32, b_ep%ord
               b_ep%b(b_ep%ord+i1+i2,:) = b_ep%b(b_ep%ord+i1+i2-1,:) / 2.0_dp
            enddo ! i1 (b_ep%ord)
         endif
         if ( jc_opt ) then
            b_ep%n_opt_par = b_ep%ord + 1_int32
            if ( b_ep%typ.eq.4 ) b_ep%n_opt_par = b_ep%n_opt_par + 1_int32
            if ( b_ep%typ.ne.3_int32 .and. jce_opt) b_ep%n_opt_par = b_ep%n_opt_par + b_ep%ord
         else
            b_ep%n_opt_par = 0_int32
         endif
      endif
      if ( jc_opt ) then
         n_par_opt = b_ep%n_opt_par
      else
         n_par_opt = 0_int32
      endif
#ifdef _PRNTALL
      if ( mpi_rank.eq.0_int32 ) then
         write(*,'(3X,"  Type of function                                  : ",I10)') b_ep%typ
         write(*,'(3X,"  Order of the polynomial                           : ",I10)') b_ep%ord
         write(*,'(3X,"  Check presence of cusp parameter                  : ",L10)') jcep_prs
         write(*,'(3X,"  Total Number of parameters                        : ",I10)') b_ep%n_par
         write(*,'(3X,"  Total Number of parameters to be optimized        : ",I10)') n_par_opt
         write(*,'(3X,"                                                                ")')
      endif
#endif
   end subroutine init_electrons_positrons_cusps
   subroutine comp_electrons_positrons_cusps( n_el, n_po, d_ep, b_ep, ep_c )
      integer(int32),                      intent(in)  :: n_el
      integer(int32),                      intent(in)  :: n_po
      type(dist_t),  dimension(n_po,n_el), intent(in)  :: d_ep
      type(jcspar_t),                      intent(in)  :: b_ep
      real(dp),   dimension(n_po,n_el),   intent(out) :: ep_c
      integer(int32) :: i1, i2
      do i1 = 1_int32, n_el ; do i2 = 1_int32, n_po
            ep_c(i2,i1) = cspfun( -mu_ep, 1, 0, b_ep, d_ep(i2,i1)%m )
         enddo; enddo
   end subroutine comp_electrons_positrons_cusps
   subroutine variation_electrons_positrons_cusps( n_fe, d_ep_new, b_ep, ep_c_new )
      integer(int32),                   intent(in) :: n_fe
      type(dist_t),  dimension(n_fe), intent(in)    :: d_ep_new
      type(jcspar_t),                   intent(in)    :: b_ep
      real(dp),      dimension(n_fe), intent(inout) :: ep_c_new
      integer(int32) :: i1
      do i1 = 1_int32, n_fe
         ep_c_new(i1) = cspfun( -mu_ep, 1, 0, b_ep, d_ep_new(i1)%m )
      enddo 
   end subroutine variation_electrons_positrons_cusps
   subroutine DD2_electrons_positrons_cusps( n_el, n_po, d_ep, b_ep, DD2_ep_c )
      integer(int32),                          intent(in)  :: n_el, n_po
      type(dist_t),  dimension(n_po,n_el), intent(in)  :: d_ep
      type(jcspar_t),                          intent(in)  :: b_ep
      real(dp),     dimension(4*n_el,n_po),       intent(inout) :: DD2_ep_c
      integer(int32) :: i1, i2
      do i2 = 1_int32, n_po
         do i1 = 1_int32, n_el
            DD2_ep_c(4*i1-3:4*i1,i2) = DD2_cspfun( -mu_ep, 1, 0, b_ep, d_ep(i2,i1)%v, d_ep(i2,i1)%m )
         enddo ; enddo 
   end subroutine DD2_electrons_positrons_cusps
   subroutine new_DD2_electrons_positrons_cusps( n_fe, d_ep_new, b_ep, DD2_ep_c_new )
      integer(int32),                   intent(in)    :: n_fe
      type(dist_t),  dimension(n_fe), intent(in)    :: d_ep_new
      type(jcspar_t),                   intent(in)    :: b_ep
      real(dp), dimension(4*n_fe),    intent(inout) :: DD2_ep_c_new
      integer(int32) :: i1
      do i1 = 1_int32, n_fe
         DD2_ep_c_new(4*i1-3:4*i1) = DD2_cspfun( -mu_ep, 1, 0, b_ep, d_ep_new(i1)%v, d_ep_new(i1)%m )
      enddo 
   end subroutine new_DD2_electrons_positrons_cusps
   subroutine updt_DD2_electrons_positrons_cusps( i_fe, n_fe, n_el, n_po, DD2_ep_c_new, DD2_ep_c )
      integer(int32),                       intent(in)    :: i_fe
      integer(int32),                       intent(in)    :: n_fe
      integer(int32),                       intent(in)    :: n_el, n_po
      real(dp), dimension(4*n_fe),        intent(in)    :: DD2_ep_c_new
      real(dp), dimension(4*n_el,n_po), intent(inout) :: DD2_ep_c
      integer(int32) :: i1
      if (i_fe.le.n_el) then
         do i1 = 1_int32, n_po
            DD2_ep_c(4*i_fe-3:4*i_fe,i1) = DD2_ep_c_new(4*i1-3:4*i1)
         enddo
      else
         do i1 = 1_int32, n_el
            DD2_ep_c(4*i1-3:4*i1-1,i_fe-n_el) = -DD2_ep_c_new(4*i1-3:4*i1-1)
            DD2_ep_c(4*i1,i_fe-n_el) = DD2_ep_c_new(4*i1)
         enddo
      endif
   end subroutine updt_DD2_electrons_positrons_cusps
   subroutine db_electrons_positrons_cusps( n_el, n_po, d_ep, b_ep, jce_opt, n_par_jepc, db_ep_c, db_jst )
      integer(int32),                              intent(in)  :: n_el, n_po
      type(dist_t),  dimension(n_po,n_el),     intent(in)  :: d_ep
      type(jcspar_t),                              intent(in)  :: b_ep
      logical,        intent(in) :: jce_opt
      integer(int32), intent(in) :: n_par_jepc
      real(dp),      dimension(n_par_jepc,n_el*n_po), intent(inout) :: db_ep_c
      real(dp),      dimension(n_par_jepc),            intent(inout) :: db_jst
      integer(int32) :: i1, i2
      do i1 = 1_int32, n_el ; do i2 = 1_int32, n_po
            db_ep_c(:,(i1-1)*n_po+i2) = db_cspfun( -mu_ep, 1, 0, b_ep, d_ep(i2,i1)%m, jce_opt )
            db_jst(:) = db_jst(:) + db_ep_c(:,(i1-1)*n_po+i2)
         enddo ; enddo 
   end subroutine db_electrons_positrons_cusps
   subroutine updt_electrons_positrons_cusps_params( b_ep, jce_opt, n_par_jepc, vec_jst_var )
      type(jcspar_t), intent(inout) :: b_ep
      logical,        intent(in)    :: jce_opt
      integer(int32), intent(in)    :: n_par_jepc
      real(dp),       intent(in)    :: vec_jst_var(1:n_par_jepc)
      real(dp) :: z_tmp
      integer(int32) :: ip, i1
      ip = 0_int32
      z_tmp = b_ep%b(1,1) + vec_jst_var(ip+1_int32)
      if ( b_ep%typ.eq.4_int32 ) then
         if( z_tmp.lt.0.0_dp ) b_ep%b(1,1) = z_tmp
      else
         if( z_tmp.gt.0.0_dp ) b_ep%b(1,1) = z_tmp
      endif
      ip = ip + 1_int32
      if ( b_ep%typ.eq.4_int32 ) then
         z_tmp = b_ep%b(2,1) + vec_jst_var(ip+1_int32)
         if( z_tmp.gt.0.0_dp ) b_ep%b(2,1) = z_tmp
         ip = ip + 1_int32
      endif
      if ( b_ep%ord.gt.0_int32 ) then
         if ( b_ep%typ.eq.4_int32 ) then
            b_ep%b(3:b_ep%ord+2,1) = b_ep%b(3:b_ep%ord+2,1) &
            & + vec_jst_var(ip+1:ip+b_ep%ord)
            ip = ip + b_ep%ord
            if ( jce_opt ) then
               do i1 = 1_int32, b_ep%ord
                  z_tmp = b_ep%b(b_ep%ord+2+i1,1) + vec_jst_var(ip+i1)
                  if( z_tmp.gt.0.0_dp ) b_ep%b(b_ep%ord+2+i1,1) = z_tmp
               enddo
            endif
         else
            b_ep%b(2:b_ep%ord+1,1) = b_ep%b(2:b_ep%ord+1,1) &
            & + vec_jst_var(ip+1:ip+b_ep%ord)
            ip = ip + b_ep%ord
            if ( b_ep%typ.ne.3_int32.and.jce_opt ) then
               do i1 = 1_int32, b_ep%ord
                  z_tmp = b_ep%b(b_ep%ord+1+i1,1) + vec_jst_var(ip+i1)
                  if( z_tmp.gt.0.0_dp ) b_ep%b(b_ep%ord+1+i1,1) = z_tmp
               enddo
            endif
         endif
      endif
   end subroutine updt_electrons_positrons_cusps_params
end module electrons_positrons_cusps_m
