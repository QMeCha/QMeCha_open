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
module fermions_fermions_cusps_m
   use fortran_kinds_v, only: sp, dp, int32
   use openmp_mpi_m, only: mpi_rank
   use cusps_functions_params_m, only: jcspar_t
   use molecular_system_v, only: dist_t
   use cusp_functions_m,   only: cspfun, DD2_cspfun, db_cspfun
   implicit none
   public :: init_fermions_fermions_cusps, comp_fermions_fermions_cusps, &
   & variation_fermions_fermions_cusps, DD2_fermions_fermions_cusps, &
   & new_DD2_fermions_fermions_cusps, updt_DD2_fermions_fermions_cusps, &
   & db_fermions_fermions_cusps
contains
   subroutine init_fermions_fermions_cusps( n_fe, n_fe_u, n_fe_d, b_ff, jcff_prs, n_par_opt, jc_opt, jce_opt )
      integer(int32), intent(in)    :: n_fe, n_fe_u, n_fe_d
      type(jcspar_t), intent(inout) :: b_ff
      logical,        intent(inout) :: jcff_prs
      integer(int32), intent(inout) :: n_par_opt
      logical,        intent(in)    :: jc_opt, jce_opt
      logical                       :: b_ff_a_prs, b_ff_p_prs
      integer(int32)                :: i1
      if ( n_fe.lt.2_int32 ) b_ff%typ = 0_int32
      if ( b_ff%typ.eq.0_int32 ) then
         b_ff%ord       = 0_int32
         b_ff%n_par     = 0_int32
         b_ff%n_opt_par = 0_int32
         jcff_prs       = .false.
         b_ff_a_prs =.false.
         b_ff_p_prs =.false.
      else
         jcff_prs = .true.
         if (n_fe_u.ge.2_int32 .or. n_fe_d.ge.2_int32) then
            b_ff_p_prs = .true.
         else
            b_ff_p_prs = .false.
         endif
         if (n_fe_u.gt.0_int32 .and. n_fe_d.gt.0_int32) then
            b_ff_a_prs = .true.
         else
            b_ff_a_prs = .false.
         endif
         if ( b_ff%typ.ne.3 ) then
            b_ff%n_par = 2_int32*b_ff%ord + 1_int32
         else
            b_ff%n_par = b_ff%ord + 1_int32
         endif
         b_ff%n_cmp = 0_int32
         if ( b_ff_a_prs ) b_ff%n_cmp = b_ff%n_cmp + 1_int32
         if ( b_ff_p_prs ) b_ff%n_cmp = b_ff%n_cmp + 1_int32
         allocate(b_ff%b(1:b_ff%n_par,1:b_ff%n_cmp)) ; b_ff%b = 0.0_dp
         b_ff%b(1,:) = 1.0_dp
         if ( b_ff%ord.gt.0 .and. b_ff%typ.ne.3 ) then
            b_ff%b(b_ff%ord+2,:) = 1.0_dp
            if ( b_ff%ord.gt.1 ) then
               do i1 = 2_int32, b_ff%ord
                  b_ff%b(b_ff%ord+i1+1,:) = b_ff%b(b_ff%ord+i1,:) / 2.0_dp
               enddo 
            endif
         endif
         if ( jc_opt ) then
            if ( b_ff%ord.gt.0_int32 ) then
               b_ff%n_opt_par = b_ff%ord + 1_int32
               if ( jce_opt.and.b_ff%typ.ne.3 ) b_ff%n_opt_par = b_ff%n_opt_par + b_ff%ord
            else
               b_ff%n_opt_par = b_ff%ord + 1_int32
            endif
         else
            b_ff%n_opt_par = 0_int32
         endif
      endif
      if (jc_opt.and.jcff_prs) then
         if (b_ff_a_prs.and.b_ff_p_prs) then
            n_par_opt = 2 * b_ff%n_opt_par
         else
            n_par_opt = b_ff%n_opt_par
         endif
      else
         n_par_opt = 0_int32
      endif
#ifdef _PRNTALL
      if ( mpi_rank.eq.0_int32 ) then
         write(*,'(3X,"  Type of function                                  : ",I10)') b_ff%typ
         write(*,'(3X,"  Order of the polynomial                           : ",I10)') b_ff%ord
         write(*,'(3X,"  Check presence of antiparallel/parallel coupling  : ",2L5)') b_ff_a_prs, b_ff_p_prs
         write(*,'(3X,"  Total Number of parameters                        : ",I10)') b_ff%n_par
         write(*,'(3X,"  Total Number of parameters to be optimized        : ",I10)') n_par_opt
         write(*,'(3X,"                                                                ")')
      endif
#endif
   end subroutine init_fermions_fermions_cusps
   subroutine comp_fermions_fermions_cusps( n_fe, d_ff, spin, b_ff, ff_c )
      integer(int32),                            intent(in)  :: n_fe
      type(dist_t),  dimension(n_fe*(n_fe-1)/2), intent(in)  :: d_ff
      integer(int32),      dimension(n_fe),            intent(in)  :: spin
      type(jcspar_t),                              intent(in)  :: b_ff
      real(dp),      dimension(n_fe*(n_fe-1)/2), intent(inout) :: ff_c
      integer(int32) :: i1, i2, i3, ia
      if (b_ff%n_cmp.eq.2_int32) then
         ia = 2
      else
         ia = 1
      endif
      i3 = 1_int32
      do i1 = 1_int32, n_fe ; do i2 = i1 + 1_int32, n_fe
         if ( spin(i1).eq.spin(i2) ) then
            ff_c(i3) = cspfun( 0.25_dp, 1, 0, b_ff, d_ff(i3)%m )
         else
            ff_c(i3) = cspfun( 0.50_dp, ia, 0, b_ff, d_ff(i3)%m )
         endif
         i3 = i3 + 1_int32
      enddo; enddo
   end subroutine comp_fermions_fermions_cusps
   subroutine variation_fermions_fermions_cusps( i_fe, n_fe, d_ff_new, spin, b_ff, ff_c_new )
      integer(int32),                    intent(in)  :: i_fe
      integer(int32),                    intent(in)    :: n_fe
      type(dist_t),   dimension(n_fe), intent(in)    :: d_ff_new
      integer(int32), dimension(n_fe), intent(in)    :: spin
      type(jcspar_t),                    intent(in)    :: b_ff
      real(dp),      dimension(n_fe),  intent(inout) :: ff_c_new
      integer(int32) :: i1, ia
      if (b_ff%n_cmp.eq.2_int32) then
         ia = 2
      else
         ia = 1
      endif
      do i1 = 1_int32, n_fe
         if ( i1.ne.i_fe ) then
            if ( spin(i1).eq.spin(i_fe) ) then
               ff_c_new(i1)= cspfun( 0.25_dp, 1, 0, b_ff, d_ff_new(i1)%m )
            else
               ff_c_new(i1) = cspfun( 0.50_dp, ia, 0, b_ff, d_ff_new(i1)%m )
            endif
         else
            ff_c_new(i1) = 0.0_dp
         endif
      enddo 
   end subroutine variation_fermions_fermions_cusps
   subroutine DD2_fermions_fermions_cusps( n_fe, d_ff, spin, b_ff, DD2_ff )
      integer(int32),                             intent(in)    :: n_fe
      type(dist_t), dimension(n_fe*(n_fe-1)/2), intent(in)    :: d_ff
      integer(int32), dimension(n_fe),          intent(in)    :: spin
      type(jcspar_t),                             intent(in)    :: b_ff
      real(dp),     dimension(4*n_fe,n_fe),   intent(inout) :: DD2_ff
      integer(int32) :: i1, i2, i3, ia
      if (b_ff%n_cmp.eq.2_int32) then
         ia = 2
      else
         ia = 1
      endif
      i3 = 1_int32
      do i2 = 1_int32, n_fe ; do i1 = i2 + 1_int32, n_fe
         if ( spin(i1).eq.spin(i2) ) then
            DD2_ff(4*i1-3:4*i1,i2) = DD2_cspfun( 0.25_dp, 1, 0, b_ff, -d_ff(i3)%v(:), d_ff(i3)%m )
         else
            DD2_ff(4*i1-3:4*i1,i2) = DD2_cspfun( 0.50_dp, ia, 0, b_ff, -d_ff(i3)%v(:), d_ff(i3)%m )
         endif
         DD2_ff(4*i2-3:4*i2-1,i1) = - DD2_ff(4*i1-3:4*i1-1,i2)
         DD2_ff(4*i2,i1)          =   DD2_ff(4*i1,i2)
         i3 = i3 + 1
      enddo ; enddo
   end subroutine DD2_fermions_fermions_cusps
   subroutine new_DD2_fermions_fermions_cusps( i_fe, n_fe, d_ff_new, spin, b_ff, DD2_ff_c_new )
      integer(int32),                              intent(in)  :: i_fe
      integer(int32),                              intent(in)    :: n_fe
      type(dist_t),  dimension(n_fe),            intent(in)    :: d_ff_new
      integer(int32),      dimension(n_fe),      intent(in)    :: spin
      type(jcspar_t),                              intent(in)    :: b_ff
      real(dp),      dimension(4*n_fe),          intent(inout) :: DD2_ff_c_new
      integer(int32) :: i1, ia
      if (b_ff%n_cmp.eq.2_int32) then
         ia = 2
      else
         ia = 1
      endif
      do i1 = 1_int32, n_fe
         if (i1.ne.i_fe) then
            if ( spin(i1).eq.spin(i_fe) ) then
               DD2_ff_c_new(4*i1-3:4*i1) = DD2_cspfun( 0.25_dp, 1, 0, b_ff, d_ff_new(i1)%v(:), d_ff_new(i1)%m )
            else
               DD2_ff_c_new(4*i1-3:4*i1) = DD2_cspfun( 0.50_dp, ia, 0, b_ff, d_ff_new(i1)%v(:), d_ff_new(i1)%m )
            endif
         else
            DD2_ff_c_new(4*i1-3:4*i1) = 0.0_dp
         endif
      enddo
   end subroutine new_DD2_fermions_fermions_cusps
   subroutine updt_DD2_fermions_fermions_cusps( i_fe, n_fe, DD2_ff_c_new, DD2_ff_c )
      integer(int32),                            intent(in)    :: i_fe
      integer(int32),                            intent(in)    :: n_fe
      real(dp),     dimension(4*n_fe),         intent(in)    :: DD2_ff_c_new
      real(dp),     dimension(4*n_fe,n_fe),  intent(inout) :: DD2_ff_c
      integer(int32) :: i1
      do i1 = 1_int32, n_fe
         DD2_ff_c(4*i_fe-3:4*i_fe,i1) = DD2_ff_c_new(4*i1-3:4*i1)
         DD2_ff_c(4*i1-3:4*i1-1,i_fe) = -DD2_ff_c_new(4*i1-3:4*i1-1)
         DD2_ff_c(4*i1,i_fe) = DD2_ff_c_new(4*i1)
      enddo
   end subroutine updt_DD2_fermions_fermions_cusps
   subroutine db_fermions_fermions_cusps( n_fe, spin, d_ff, b_ff, jce_opt, n_par_jffc, db_ff_c, db_jst )
      integer(int32),                               intent(in)  :: n_fe
      type(dist_t),   dimension(n_fe*(n_fe-1)/2), intent(in)  :: d_ff
      integer(int32), dimension(n_fe),            intent(in)  :: spin
      type(jcspar_t),                               intent(in)  :: b_ff
      integer(int32), intent(in) :: n_par_jffc
      logical,        intent(in) :: jce_opt
      real(dp), dimension(n_par_jffc,n_fe*(n_fe-1)/2), intent(inout) :: db_ff_c
      real(dp), dimension(n_par_jffc),                   intent(inout) :: db_jst
      integer(int32) :: i1, i2, i3, ia1, ia2
      if (b_ff%n_cmp.eq.2_int32) then
         ia1 = b_ff%n_opt_par ; ia2 = 2
      else
         ia1 = 0 ; ia2 = 1
      endif
      i3 = 1_int32
      do i2 = 1_int32, n_fe ; do i1 = i2 + 1_int32, n_fe
            if ( spin(i1).eq.spin(i2) ) then
               db_ff_c(1:b_ff%n_opt_par,i3) = db_cspfun( 0.25_dp, 1, 0, b_ff, d_ff(i3)%m, jce_opt )
               db_jst(1:b_ff%n_opt_par) = db_jst(1:b_ff%n_opt_par) + db_ff_c(1:b_ff%n_opt_par,i3)
            else
               db_ff_c(ia1+1:ia1+b_ff%n_opt_par,i3) = db_cspfun( 0.50_dp, ia2, 0, b_ff, d_ff(i3)%m, jce_opt )
               db_jst(ia1+1:ia1+b_ff%n_opt_par) = db_jst(ia1+1:ia1+b_ff%n_opt_par) + db_ff_c(ia1+1:ia1+b_ff%n_opt_par,i3)
            endif
            i3 = i3 + 1_int32
         enddo ; enddo 
   end subroutine db_fermions_fermions_cusps
   subroutine updt_fermions_fermions_cusps_params( b_ff, jce_opt, n_par_jffc, vec_jst_var )
      type(jcspar_t), intent(inout) :: b_ff
      logical,        intent(in)    :: jce_opt
      integer(int32), intent(in)    :: n_par_jffc
      real(dp),       intent(in)    :: vec_jst_var(1:n_par_jffc)
      real(dp) :: z_tmp
      integer(int32) :: i1, i2, ip
      ip = 0_int32
      do i1 = 1_int32, b_ff%n_cmp
         z_tmp = b_ff%b(1,i1) + vec_jst_var(ip+1)
         if( z_tmp.gt.0.0_dp ) then
            b_ff%b(1,i1) = z_tmp
         endif
         ip = ip + 1_int32
         if ( b_ff%ord.gt.0_int32 ) then
            b_ff%b(2:b_ff%ord+1,i1) = b_ff%b(2:b_ff%ord+1,i1) &
            & + vec_jst_var(ip+1:ip+b_ff%ord)
            ip = ip + b_ff%ord
            if ( jce_opt.and.b_ff%ord.ne.3 ) then
               do i2 = 1_int32, b_ff%ord
                  z_tmp = b_ff%b(b_ff%ord+1+i2,i1) + vec_jst_var(ip+i2)
                  if( z_tmp.gt.0.0_dp ) then
                     b_ff%b(b_ff%ord+1+i2,i1) = z_tmp
                  endif
               enddo
               ip = ip + b_ff%ord
            endif
         endif
      enddo 
   end subroutine updt_fermions_fermions_cusps_params
end module fermions_fermions_cusps_m
