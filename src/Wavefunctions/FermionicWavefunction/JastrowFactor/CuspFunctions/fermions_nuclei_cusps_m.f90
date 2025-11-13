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
module fermions_nuclei_cusps_m
   use fortran_kinds_v,    only: dp, int32
   use openmp_mpi_m,       only: mpi_rank
   use cusps_functions_params_m,         only: jcspar_t
   use molecular_system_v, only: dist_t, n_at, atoms, f_spin, n_indp_at, atm_sym
   use cusp_functions_m,   only: cspfun, DD2_cspfun, db_cspfun
   implicit none
   public :: init_fermions_nuclei_cusps, comp_fermions_nuclei_cusps, &
   & variation_fermions_nuclei_cusps, DD2_fermions_nuclei_cusps, &
   & new_DD2_fermions_nuclei_cusps, updt_DD2_fermions_nuclei_cusps, &
   & db_fermions_nuclei_cusps, updt_fermions_nuclei_cusps_params
contains
   subroutine init_fermions_nuclei_cusps( b_fn, jcfn_prs, n_par_opt, jc_opt, jce_opt )
      type(jcspar_t), intent(inout) :: b_fn
      logical,        intent(inout) :: jcfn_prs
      integer(int32), intent(inout) :: n_par_opt
      logical,        intent(in)    :: jc_opt, jce_opt
      integer(int32)                :: i1, i2
      if ( n_at.eq.0_int32 ) b_fn%typ = 0_int32
      if ( b_fn%typ.eq.0_int32 ) then
         b_fn%ord       = 0_int32
         b_fn%n_par     = 0_int32
         b_fn%n_opt_par = 0_int32
         jcfn_prs       = .false.
      else 
         jcfn_prs       = .true.
         if ( b_fn%typ.ne.3 ) then
            b_fn%n_par = 2_int32 * b_fn%ord + 1_int32
         else
            b_fn%n_par = b_fn%ord + 1_int32
         endif
         b_fn%n_cmp = n_at
         allocate( b_fn%b(1:b_fn%n_par,1:b_fn%n_cmp)) ; b_fn%b = 0.0_dp
         if ( b_fn%typ.ne.3 .and. b_fn%ord.gt.0_int32 ) then
            do i2 = 1_int32, b_fn%n_cmp
               if ( atoms(i2)%i_pp.eq.0_int32 ) then
                  b_fn%b(1,i2) = 1.0_dp
                  b_fn%b(b_fn%ord+2,i2) = dble(atoms(i2)%atm_z)
                  do i1 = 2_int32, b_fn%ord
                     b_fn%b(b_fn%ord+i1+1,i2) = b_fn%b(b_fn%ord+i1,i2) / 4.0_dp
                  enddo 
               else
                  b_fn%b(1,i2) = 0.0_dp
                  b_fn%b(b_fn%ord+2,i2) = sqrt(dble(atoms(i2)%atm_z) /2.0_dp)
                  do i1 = 2_int32, b_fn%ord
                     b_fn%b(b_fn%ord+i1+1,i2) = b_fn%b(b_fn%ord+i1,i2) / 2.0_dp
                  enddo 
               endif
            enddo
         else
            b_fn%b(1,:) = 1.0_dp
         endif
         if ( jc_opt ) then
            b_fn%n_opt_par = b_fn%ord + 1_int32
            if ( b_fn%ord.gt.0_int32 .and. b_fn%typ.ne.3 .and. jce_opt ) then
               b_fn%n_opt_par = b_fn%n_opt_par + b_fn%ord
            endif
         else
            b_fn%n_opt_par = 0_int32
         endif
      endif
      if ( jc_opt ) then
         n_par_opt = b_fn%n_opt_par*n_indp_at
      else
         n_par_opt = 0_int32
      endif
#ifdef _PRNTALL
      if ( mpi_rank.eq.0_int32 ) then
         write(*,'(3X,"  Type of function                                  : ",I10)') b_fn%typ
         write(*,'(3X,"  Order of the polynomial                           : ",I10)') b_fn%ord
         write(*,'(3X,"  Check presence of cusp parameter                  : ",L10)') jcfn_prs
         write(*,'(3X,"  Total Number of parameters in cusp function       : ",I10)') b_fn%n_par
         write(*,'(3X,"  Total Number of parameters to be optimized        : ",I10)') n_par_opt
         write(*,'(3X,"                                                                ")')
      endif
#endif
   end subroutine init_fermions_nuclei_cusps
   subroutine comp_fermions_nuclei_cusps( chrg, n_fe, d_fn, b_fn, fn_c )
      real(dp),                           intent(in)  :: chrg
      integer(int32),                     intent(in)  :: n_fe
      type(dist_t), dimension(n_at,n_fe), intent(in)  :: d_fn
      type(jcspar_t),                     intent(in)  :: b_fn
      real(dp),      dimension(n_fe),     intent(out) :: fn_c
      integer(int32) :: i1, i2
      do i1 = 1_int32, n_fe ; do i2 = 1_int32, n_at
         fn_c(i1) = fn_c(i1) + cspfun( chrg * atoms(i2)%atm_z , i2, atoms(i2)%i_pp, b_fn, d_fn(i2,i1)%m )
      enddo; enddo
   end subroutine comp_fermions_nuclei_cusps
   subroutine variation_fermions_nuclei_cusps( chrg, d_fn_new, b_fn, fn_c_new )
      real(dp),                        intent(in)    :: chrg
      type(dist_t), dimension(n_at), intent(in)    :: d_fn_new
      type(jcspar_t),                  intent(in)    :: b_fn
      real(dp),                        intent(inout) :: fn_c_new
      integer(int32) :: i1
      fn_c_new = cspfun( chrg * atoms(1)%atm_z , 1, atoms(1)%i_pp, b_fn, d_fn_new(1)%m )
      do i1 = 2_int32, n_at
         fn_c_new = fn_c_new + cspfun( chrg * atoms(i1)%atm_z , i1, atoms(i1)%i_pp, b_fn, d_fn_new(i1)%m )
      enddo 
   end subroutine variation_fermions_nuclei_cusps
   subroutine DD2_fermions_nuclei_cusps( chrg, n_fe, d_fn, b_fn, DD2_fn_c )
      real(dp),                           intent(in)    :: chrg
      integer(int32),                     intent(in)    :: n_fe
      type(dist_t), dimension(n_at,n_fe), intent(in)    :: d_fn
      type(jcspar_t),                     intent(in)    :: b_fn
      real(dp),     dimension(4*n_fe),    intent(inout) :: DD2_fn_c
      integer(int32) :: i1, i2
      do i1 = 1_int32, n_fe
         DD2_fn_c(4*i1-3:4*i1) = DD2_cspfun( chrg * atoms(1)%atm_z, 1, atoms(1)%i_pp, b_fn, d_fn(1,i1)%v, d_fn(1,i1)%m )
         do i2 = 2_int32, n_at
            DD2_fn_c(4*i1-3:4*i1) = DD2_fn_c(4*i1-3:4*i1) + DD2_cspfun( chrg * atoms(i2)%atm_z, i2, atoms(i2)%i_pp, b_fn, d_fn(i2,i1)%v, d_fn(i2,i1)%m )
      enddo ; enddo 
   end subroutine DD2_fermions_nuclei_cusps
   subroutine new_DD2_fermions_nuclei_cusps( chrg, d_fn_new, b_fn, DD2_fn_c_new )
      real(dp),                        intent(in)    :: chrg
      type(dist_t), dimension(n_at), intent(in)    :: d_fn_new
      type(jcspar_t),                  intent(in)    :: b_fn
      real(dp),     dimension(4),    intent(inout) :: DD2_fn_c_new
      integer(int32) :: i1
      DD2_fn_c_new(1:4) = DD2_cspfun( chrg * atoms(1)%atm_z, 1, atoms(1)%i_pp, b_fn, d_fn_new(1)%v, d_fn_new(1)%m )
      do i1 = 2_int32, n_at
         DD2_fn_c_new(1:4) = DD2_fn_c_new(1:4) + DD2_cspfun( chrg * atoms(i1)%atm_z, i1, atoms(i1)%i_pp, b_fn, d_fn_new(i1)%v, d_fn_new(i1)%m )
      enddo 
   end subroutine new_DD2_fermions_nuclei_cusps
   subroutine updt_DD2_fermions_nuclei_cusps( i_fe, n_fe, DD2_fn_c_new, DD2_fn_c )
      integer(int32),                     intent(in)    :: i_fe
      integer(int32),                     intent(in)    :: n_fe
      real(dp),     dimension(4),       intent(in)    :: DD2_fn_c_new
      real(dp),     dimension(4*n_fe),  intent(inout) :: DD2_fn_c
      DD2_fn_c(4*i_fe-3:4*i_fe) = DD2_fn_c_new(:)
   end subroutine updt_DD2_fermions_nuclei_cusps
   subroutine db_fermions_nuclei_cusps( chrg, n_fe, d_fn, b_fn, jce_opt, n_par_jfnc, db_fn_c, db_jst )
      real(dp),                                 intent(in)    :: chrg
      integer(int32),                           intent(in)    :: n_fe
      type(dist_t),  dimension(n_at,n_fe),  intent(in)    :: d_fn
      type(jcspar_t),                           intent(in)    :: b_fn
      logical,                                  intent(in)    :: jce_opt
      integer(int32),                           intent(in)    :: n_par_jfnc
      real(dp), dimension(n_par_jfnc,n_fe), intent(inout) :: db_fn_c
      real(dp), dimension(n_par_jfnc),        intent(inout) :: db_jst
      integer(int32) :: i1, i2, i3, ia
      do i1 = 1_int32, n_fe
         db_fn_c(:,i1) = 0.0_dp
         do i2 = 1_int32, n_indp_at ; do i3 = 1_int32, atm_sym(i2)%n_cnct_at
               ia = atm_sym(i2)%cnct(i3)
               db_fn_c(b_fn%n_opt_par*(i2-1)+1:b_fn%n_opt_par*i2,i1) = &
               & db_fn_c(b_fn%n_opt_par*(i2-1)+1:b_fn%n_opt_par*i2,i1) + &
               & db_cspfun( chrg * atoms(ia)%atm_z, ia, atoms(ia)%i_pp,  b_fn, d_fn(ia,i1)%m, jce_opt )
            enddo ; enddo 
         db_jst(:) = db_jst(:) + db_fn_c(:,i1)
      enddo 
   end subroutine db_fermions_nuclei_cusps
   subroutine updt_fermions_nuclei_cusps_params( b_fn, jce_opt, n_par_jfnc, vec_jst_var )
      type(jcspar_t), intent(inout) :: b_fn
      logical,        intent(in)    :: jce_opt
      integer(int32), intent(in)    :: n_par_jfnc
      real(dp),       intent(in)    :: vec_jst_var(1:n_par_jfnc)
      real(dp) :: z_tmp
      integer(int32) :: i1, i2, i3, ic, ip
      do i1 = 1_int32, n_indp_at
         do i2 = 1_int32, atm_sym(i1)%n_cnct_at
            ic = atm_sym(i1)%cnct(i2)
            ip = b_fn%n_opt_par*(i1-1)
            z_tmp = b_fn%b(1,ic) + vec_jst_var(ip+1)
            if( z_tmp .gt.0.0_dp ) then
               b_fn%b(1,ic) = z_tmp
            endif
            ip = ip + 1_int32
            if ( b_fn%ord.gt.0_int32 ) then
               b_fn%b(2:b_fn%ord+1,ic) = b_fn%b(2:b_fn%ord+1,ic) &
               & + vec_jst_var(ip+1:ip+b_fn%ord)
               ip = ip + b_fn%ord
               if ( b_fn%typ.ne.3_int32.and.jce_opt ) then
                  do i3 = 1_int32, b_fn%ord
                     z_tmp = b_fn%b(b_fn%ord+1+i3,ic) + vec_jst_var(ip+i3)
                     if( z_tmp.gt.0.0001_dp ) then
                        b_fn%b(b_fn%ord+1+i3,ic) = z_tmp
                     else
                        b_fn%b(b_fn%ord+1+i3,ic) = 0.0001_dp
                     endif
                  enddo
                  ip = ip + b_fn%ord
               endif
            endif
         enddo
      enddo  
   end subroutine updt_fermions_nuclei_cusps_params
end module fermions_nuclei_cusps_m
