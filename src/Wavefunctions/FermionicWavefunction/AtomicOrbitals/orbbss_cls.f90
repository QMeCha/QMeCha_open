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
module orbbss_cls
   use fortran_kinds_v,   only: int32, dp, stdout
   use openmp_mpi_m
   use atomic_orbital_c,        only: orb_t
   use basisset_v,        only: n_fbs, bssset_t
   use molecular_system_v,        only: n_at, atoms
   use write_lines_m,     only: write_simple_line
   use merge_sort_m,      only: mrg_srt
   use molecular_system_m,        only: n_indp_at, n_sys_sym, sys_sym_tab
   use molec_symmetry_m,  only: sym_op
   use cubic_harmonics_m, only: y_l_sym_tab
   implicit none
   type, public :: sym_t
      integer(int32) :: n_cnct
      integer(int32), allocatable, dimension(:) :: cnct
   end type sym_t
   type, public :: sym_tab_t
      integer(int32) :: n_indp
      type(sym_t), allocatable, dimension(:) :: indp
   end type sym_tab_t
   type, public :: orbspa_t
      integer(int32) :: n_orbs
      integer(int32) :: l_max
      type(orb_t), allocatable, dimension(:) :: orb
   end type orbspa_t
   type, public :: sysorb_t
      logical, pointer                          :: oc_opt
      logical, pointer                          :: oe_opt
      integer(int32)                            :: n_orbs
      integer(int32)                            :: n_par_orbs
      integer(int32)                            :: n_par_orbs_s
      type(orbspa_t),    allocatable, dimension(:) :: orbs
      type(sym_tab_t)                           :: orbs_sym_tab
      type(sym_tab_t)                           :: par_sym_tab
   contains
      procedure :: ini     => ini_sys_orbs
      procedure :: sym     => apl_sym_orbs_vec
      procedure :: upd     => upd_sys_orbs
   end type sysorb_t
   private :: ini_sys_orbs, cmp_orbs_sym_tab, cmp_opar_sym_tab, cnt_orbs_par, &
   & upd_sys_orbs, apl_sym_orbs_vec
   public  :: cmp_orbs_sym, cmp_opar_sym, cmp_orbs_trsf
contains
   subroutine ini_sys_orbs( obj, bsst, dc_opt, de_opt )
      class(sysorb_t), intent(inout) :: obj
      type(bssset_t), dimension(n_fbs),         intent(in) :: bsst
      logical,                           target, intent(in) :: dc_opt, de_opt
      integer(int32) :: io, ia
      obj%oc_opt => dc_opt
      obj%oe_opt => de_opt
      if ( obj%oe_opt ) obj%oc_opt = .true.
      allocate( obj%orbs(1:n_at) )
      io = 1_int32
      do ia = 1_int32, n_at
         call init_orbitals_per_atom( obj%orbs(ia), ia, io, bsst )
      enddo ! i1 (n_at)
      obj%n_orbs = sum( obj%orbs(:)%n_orbs )
      obj%n_par_orbs   = 0_int32
      obj%n_par_orbs_s = 0_int32
      if ( obj%n_orbs.gt.0_int32 ) then
         call cmp_orbs_sym_tab( obj )
         if ( obj%oc_opt.or.obj%oe_opt ) then
            call cnt_orbs_par( obj )
            call cmp_opar_sym_tab( obj )
         endif
      endif
   end subroutine ini_sys_orbs
   subroutine init_orbitals_per_atom( orbs_per_atom, i_at, i_orbs, bsst )
      type(orbspa_t),                     intent(inout) :: orbs_per_atom
      integer(int32),                     intent(in)    :: i_at
      integer(int32),                     intent(inout) :: i_orbs
      type(bssset_t), dimension(n_fbs), intent(in)    :: bsst
      integer(int32) :: i1, i2, i3, l, n
      orbs_per_atom%n_orbs = 0_int32
      do i1 = 1_int32, n_fbs
         if(trim(bsst(i1)%atm_name).eq.trim(atoms(i_at)%atm_name)) then
            orbs_per_atom%n_orbs = bsst(i1)%n_orbs
            atoms(i_at)%i_bs   = i1
            allocate( orbs_per_atom%orb(1:orbs_per_atom%n_orbs) )
         endif
      enddo 
      orbs_per_atom%l_max = 0_int32
      i3 = 1_int32
      do i1 = 1_int32, bsst(atoms(i_at)%i_bs)%n_shlls
         l = bsst(atoms(i_at)%i_bs)%shlls(i1)%l
         n = bsst(atoms(i_at)%i_bs)%shlls(i1)%n
         if( orbs_per_atom%l_max.lt.l ) orbs_per_atom%l_max = l
         do i2 = l**2, l*(l+2)
            allocate( orbs_per_atom%orb(i3)%c(1:n), orbs_per_atom%orb(i3)%z(1:n) )
            orbs_per_atom%orb(i3)       = bsst(atoms(i_at)%i_bs)%shlls(i1)
            orbs_per_atom%orb(i3)%l_z   = i2
            orbs_per_atom%orb(i3)%i_orb = i_orbs
            i_orbs = i_orbs + 1_int32
            i3 = i3 + 1_int32
         enddo 
      enddo 
   end subroutine init_orbitals_per_atom
   subroutine cmp_orbs_sym_tab( obj )
      class(sysorb_t), intent(inout) :: obj
      type(sym_tab_t) :: orbs_sym_tab_t
      integer(int32)  :: i1
      orbs_sym_tab_t%n_indp = obj%n_orbs
      allocate( orbs_sym_tab_t%indp(1:obj%n_orbs) )
      do i1 = 1_int32, obj%n_orbs
         orbs_sym_tab_t%indp(i1)%n_cnct = 0_int32
         allocate(orbs_sym_tab_t%indp(i1)%cnct(1:obj%n_orbs)) ; orbs_sym_tab_t%indp(i1)%cnct = 0_int32
      enddo 
      call cmp_orbs_sym( obj%orbs, orbs_sym_tab_t )
      obj%orbs_sym_tab%n_indp = orbs_sym_tab_t%n_indp
      allocate( obj%orbs_sym_tab%indp(1:obj%orbs_sym_tab%n_indp) )
      do i1 = 1_int32, obj%orbs_sym_tab%n_indp
         obj%orbs_sym_tab%indp(i1)%n_cnct = orbs_sym_tab_t%indp(i1)%n_cnct
         allocate( obj%orbs_sym_tab%indp(i1)%cnct(1:obj%orbs_sym_tab%indp(i1)%n_cnct))
         obj%orbs_sym_tab%indp(i1)%cnct(1:obj%orbs_sym_tab%indp(i1)%n_cnct) = orbs_sym_tab_t%indp(i1)%cnct(1:obj%orbs_sym_tab%indp(i1)%n_cnct)
      enddo
      deallocate( orbs_sym_tab_t%indp )
   end subroutine cmp_orbs_sym_tab
   subroutine cnt_orbs_par( obj )
      class(sysorb_t), intent(inout) :: obj
      integer(int32) :: i1, i2, i3
      obj%n_par_orbs = 0_int32
      if ( obj%oc_opt ) then
         i3 = 1_int32
         do i1 = 1_int32, n_at ; do i2 = 1_int32, obj%orbs(i1)%n_orbs
               obj%orbs(i1)%orb(i2)%i_par = i3
               if ( obj%oe_opt ) then
                  if( obj%orbs(i1)%orb(i2)%n.gt.1_int32 ) then
                     obj%orbs(i1)%orb(i2)%n_par = 2 * obj%orbs(i1)%orb(i2)%n
                  else
                     obj%orbs(i1)%orb(i2)%n_par = 1_int32
                  endif
               else
                  if ( obj%orbs(i1)%orb(i2)%n.gt.1_int32 ) then
                     obj%orbs(i1)%orb(i2)%n_par = obj%orbs(i1)%orb(i2)%n
                  else
                     obj%orbs(i1)%orb(i2)%n_par = 0_int32
                  endif
               endif
               i3 = i3 + obj%orbs(i1)%orb(i2)%n_par
               obj%n_par_orbs = obj%n_par_orbs + obj%orbs(i1)%orb(i2)%n_par
            enddo ; enddo 
      endif
      if( obj%n_par_orbs.eq.0_int32 ) then
         obj%oc_opt = .false.
         obj%oe_opt = .false.
      endif
   end subroutine cnt_orbs_par
   subroutine cmp_opar_sym_tab(obj )
      class(sysorb_t), intent(inout) :: obj
      type(sym_tab_t) :: opar_sym_tab_t
      integer(int32)  :: n_cnct_max
      integer(int32)  :: i1
      n_cnct_max = obj%orbs_sym_tab%indp(1)%n_cnct
      do i1 = 2_int32, obj%orbs_sym_tab%n_indp
         if ( n_cnct_max .lt. obj%orbs_sym_tab%indp(i1)%n_cnct ) n_cnct_max = obj%orbs_sym_tab%indp(i1)%n_cnct
      enddo
      opar_sym_tab_t%n_indp = obj%orbs_sym_tab%n_indp
      allocate( opar_sym_tab_t%indp(1:opar_sym_tab_t%n_indp) )
      do i1 = 1_int32, opar_sym_tab_t%n_indp
         opar_sym_tab_t%indp(i1)%n_cnct = 0_int32
         allocate( opar_sym_tab_t%indp(i1)%cnct(1:2*n_cnct_max+2) )
         opar_sym_tab_t%indp(i1)%cnct = 0_int32
      enddo
      call cmp_opar_sym ( obj%orbs, obj%orbs_sym_tab, opar_sym_tab_t, obj%n_par_orbs_s)
      obj%par_sym_tab%n_indp = opar_sym_tab_t%n_indp
      allocate( obj%par_sym_tab%indp(1:obj%par_sym_tab%n_indp) )
      do i1 = 1_int32, obj%par_sym_tab%n_indp
         obj%par_sym_tab%indp(i1)%n_cnct = opar_sym_tab_t%indp(i1)%n_cnct
         allocate( obj%par_sym_tab%indp(i1)%cnct(1:2*obj%par_sym_tab%indp(i1)%n_cnct+2))
         obj%par_sym_tab%indp(i1)%cnct = opar_sym_tab_t%indp(i1)%cnct(1:2*obj%par_sym_tab%indp(i1)%n_cnct+2)
      enddo
      deallocate( opar_sym_tab_t%indp )
   end subroutine cmp_opar_sym_tab
   subroutine upd_sys_orbs( obj, da_s, z_cut )
      class(sysorb_t), intent(inout), target :: obj
      real(dp), dimension(obj%n_par_orbs_s), intent(in) :: da_s
      real(dp)                               , intent(in) :: z_cut
      real(dp), dimension(obj%n_par_orbs)               :: da
      real(dp)                :: z_new
      integer(int32) :: i1, i2, i3, i4
      integer(int32), pointer :: ip, fp, ips, fps
      da = 0.0_dp
      do i1 = 1_int32, obj%par_sym_tab%n_indp
         ips => obj%par_sym_tab%indp(i1)%cnct(1)
         fps => obj%par_sym_tab%indp(i1)%cnct(2)
         do i2 = 1_int32, obj%par_sym_tab%indp(i1)%n_cnct
            ip => obj%par_sym_tab%indp(i1)%cnct(2*i2+1)
            fp => obj%par_sym_tab%indp(i1)%cnct(2*i2+2)
            da(ip:fp) = da_s(ips:fps)
         enddo ; enddo
      i4 = 1_int32
      do i1 = 1_int32, n_at ; do i2 = 1_int32, obj%orbs(i1)%n_orbs
            if ( obj%orbs(i1)%orb(i2)%n_par .gt. 0_int32 ) then
               if ( obj%orbs(i1)%orb(i2)%n .eq. 1_int32 ) then
                  if( obj%oe_opt ) then
                     z_new = obj%orbs(i1)%orb(i2)%z(1) + da(i4)
                     if ( z_new .gt. z_cut ) then
                        obj%orbs(i1)%orb(i2)%z(1) = z_new
                     else
                        obj%orbs(i1)%orb(i2)%z(1) = z_cut
                     endif
                     i4 = i4 + 1_int32
                  endif
               else
                  if( obj%oe_opt ) then
                     obj%orbs(i1)%orb(i2)%c(1:obj%orbs(i1)%orb(i2)%n) &
                     & = obj%orbs(i1)%orb(i2)%c(1:obj%orbs(i1)%orb(i2)%n) &
                     & + da(i4:i4+obj%orbs(i1)%orb(i2)%n-1)
                     i4 = i4 + obj%orbs(i1)%orb(i2)%n
                     do i3 = 1_int32, obj%orbs(i1)%orb(i2)%n
                        z_new = obj%orbs(i1)%orb(i2)%z(i3) + da(i4)
                        if ( z_new .gt. z_cut ) then
                           obj%orbs(i1)%orb(i2)%z(i3) = z_new
                        else
                           obj%orbs(i1)%orb(i2)%z(i3) = z_cut
                        endif
                        i4 = i4 + 1_int32
                     enddo
                  else
                     obj%orbs(i1)%orb(i2)%c(1:obj%orbs(i1)%orb(i2)%n) &
                     & = obj%orbs(i1)%orb(i2)%c(1:obj%orbs(i1)%orb(i2)%n) &
                     & + da(i4:i4+obj%orbs(i1)%orb(i2)%n-1)
                     i4 = i4 + obj%orbs(i1)%orb(i2)%n
                  endif
               endif
            endif
         enddo ; enddo
   end subroutine upd_sys_orbs
   subroutine apl_sym_orbs_vec( obj, da, da_s )
      class(sysorb_t), target, intent(inout) :: obj
      real(dp), dimension(obj%n_par_orbs),   intent(in)    :: da
      real(dp), dimension(obj%n_par_orbs_s), intent(inout) :: da_s
      integer(int32) :: i1, i2
      integer(int32), pointer :: ip, fp, ips, fps
      do i1 = 1_int32, obj%par_sym_tab%n_indp
         ips => obj%par_sym_tab%indp(i1)%cnct(1)
         fps => obj%par_sym_tab%indp(i1)%cnct(2)
         do i2 = 1_int32, obj%par_sym_tab%indp(i1)%n_cnct
            ip => obj%par_sym_tab%indp(i1)%cnct(2*i2+1)
            fp => obj%par_sym_tab%indp(i1)%cnct(2*i2+2)
            da_s(ips:fps) = da_s(ips:fps) + da(ip:fp)
         enddo ; enddo
   end subroutine apl_sym_orbs_vec
   subroutine cmp_orbs_sym( o, o_sym )
      type(orbspa_t), dimension(n_at), intent(in)    :: o
      type(sym_tab_t),                   intent(inout) :: o_sym
      logical        :: acc, orb_already_cnct
      integer(int32) :: sgn, cnct
      integer(int32) :: i_at, i_ac, i_sym, n_cnct
      integer(int32) :: i1, i2, i3, i4, i5
      o_sym%n_indp = 0_int32
      do i1 = 1_int32, n_at
         i_at = sys_sym_tab(1,i1)
         do i2 = 1_int32, o(i_at)%n_orbs
            o_sym%n_indp = o_sym%n_indp + 1_int32
            o_sym%indp(o_sym%n_indp)%n_cnct = 1_int32
            o_sym%indp(o_sym%n_indp)%cnct(1) = o(i_at)%orb(i2)%i_orb
         enddo
      enddo 
      i5 = 0_int32
      do i1 = 1_int32, n_at
         i_at = sys_sym_tab(1,i1)
         do i2 = 1_int32, o(i_at)%n_orbs
            i5 = i5 + 1_int32
            do i3 = 2_int32, n_sys_sym
               i_ac  = sys_sym_tab(i3,i1)
               i_sym = sys_sym_tab(i3,0)
               call cmp_orbs_trsf( i_sym, o(i_ac)%orb(i2), acc, sgn, cnct )
               if ( acc ) then
                  orb_already_cnct = .false.
                  do i4 = 1_int32, o_sym%indp(i5+cnct)%n_cnct
                     if ( o_sym%indp(i5+cnct)%cnct(i4) .eq. o(i_ac)%orb(i2)%i_orb ) then
                        orb_already_cnct = .true.
                     endif
                  enddo
                  if (.not.orb_already_cnct) then
                     o_sym%indp(i5+cnct)%n_cnct = o_sym%indp(i5+cnct)%n_cnct + 1_int32
                     o_sym%indp(i5+cnct)%cnct(o_sym%indp(i5+cnct)%n_cnct) = o(i_ac)%orb(i2)%i_orb
                  endif
               endif
            enddo
         enddo 
      enddo 
      do i1 = 1_int32, o_sym%n_indp
         n_cnct = o_sym%indp(i1)%n_cnct
         call mrg_srt( n_cnct, o_sym%indp(i1)%cnct(1:n_cnct) )
      enddo
      do i1 = 1_int32, o_sym%n_indp
         if ( o_sym%indp(i1)%n_cnct.gt.0_int32 ) then
            do i2 = i1 + 1_int32, o_sym%n_indp
               if ( o_sym%indp(i1)%cnct(1).eq.o_sym%indp(i2)%cnct(1) ) o_sym%indp(i2)%n_cnct = 0_int32
            enddo
         endif
      enddo
      i1 = 1_int32
      do while ( i1.le.o_sym%n_indp )
         if ( o_sym%indp(i1)%n_cnct.eq.0_int32 ) then
            o_sym%indp(i1:o_sym%n_indp-1) = o_sym%indp(i1+1:o_sym%n_indp)
            o_sym%indp(o_sym%n_indp)%n_cnct = 0_int32
            o_sym%n_indp = o_sym%n_indp - 1_int32
         else
            i1 = i1 + 1_int32
         endif
      enddo
      n_cnct = sum( o(1:n_at)%n_orbs )-sum(o_sym%indp(1:o_sym%n_indp)%n_cnct)
      if (n_cnct .ne.0_int32) then
         if (mpi_rank.eq.0_int32) then
            call write_simple_line(stdout,0,mpi_rank,2,'l',"ERROR!!! Table of orbital symmetries is wrong.")
            call write_simple_line(stdout,0,mpi_rank,11,'l',"Check fort.12 and contact matteo.barborini@gmail.com")
            write(12,*) o_sym%n_indp
            do i1 = 1_int32, o_sym%n_indp
               write(12,*) o_sym%indp(i1)%n_cnct, o_sym%indp(i1)%cnct(1:o_sym%indp(i1)%n_cnct)
            enddo
         endif
         call fnlz_ompmpi_env()
         stop
      endif
   end subroutine cmp_orbs_sym
   subroutine cmp_orbs_trsf( i_sym , orb , acc, sgn, cnct )
      integer(int32),  intent(in)  :: i_sym
      type(orb_t),     intent(in)  :: orb
      logical,         intent(out) :: acc
      integer(int32),  intent(out) :: sgn
      integer(int32),  intent(out) :: cnct
      real(dp), dimension(4,1:4) :: y_l_ro
      integer(int32), allocatable, dimension(:) :: acc_sgn
      integer(int32), allocatable, dimension(:) :: chckd_cmp
      integer(int32) :: i1, i2, i3
      if (orb%l.eq.0_int32) then
         sgn  = 1_int32
         cnct = 0_int32
         acc = .true.
         return
      endif
      allocate( acc_sgn(1:orb%l) ) ; acc_sgn = 0_int32
      allocate( chckd_cmp(1:orb%l) ) ; chckd_cmp = 0_int32
      do i1 = 1_int32, orb%l
         call dgemv('T', 3_int32, 3_int32, 1.0_dp, dble(sym_op(:,:,i_sym)), 3_int32, &
         & y_l_sym_tab(1:3,i1,orb%l_z), 1_int32, 0.0_dp, y_l_ro(1:3,i1), 1_int32)
      enddo
      y_l_ro(4,1:orb%l) = y_l_sym_tab(4,1:orb%l,orb%l_z)
      cnct = 0_int32
      acc = .false.
      do i2 = orb%l**2, orb%l*(orb%l+2)
         acc_sgn = 0_int32
         chckd_cmp = 1_int32
         sgn = 1_int32
         do i1 = 1, orb%l
            i3loop: do i3 = 1, orb%l
               if ( chckd_cmp(i3) .eq. 1_int32) then
                  if (sum(abs(y_l_sym_tab(1:3,i1,i2) - y_l_ro(1:3,i3))).lt.10d-5) then
                     acc_sgn(i1) = 1_int32
                     chckd_cmp(i3) = 0_int32
                     exit i3loop
                  endif
                  if ( acc_sgn(i1).ne.1_int32 ) then
                     if (sum(abs(y_l_sym_tab(1:3,i1,i2) + y_l_ro(1:3,i3))).lt.10d-5) then
                        acc_sgn(i1) = -1_int32
                        chckd_cmp(i3) = 0_int32
                        exit i3loop
                     endif
                  endif
               endif
            enddo i3loop
         enddo
         do i1 = 1_int32, orb%l
            sgn = sgn * acc_sgn(i1)
         enddo
         if (sgn.ne.0_int32) acc = .true.
         if ( acc ) then
            cnct = i2 - orb%l_z
            exit
         endif
      enddo 
   end subroutine cmp_orbs_trsf
   subroutine cmp_opar_sym( o, o_sym, o_sym_par, n_indp_par )
      type(orbspa_t),  dimension(n_at), intent(in)    :: o
      type(sym_tab_t),                    intent(in)    :: o_sym
      type(sym_tab_t),                    intent(inout) :: o_sym_par
      integer(int32),                     intent(inout) :: n_indp_par
      integer(int32), allocatable, dimension(:,:) :: par_lst_tmp
      integer(int32) :: n_orbs
      logical :: check_pairs
      integer(int32) :: i1, i2, i3, i4
      n_orbs = sum ( o(:)%n_orbs )
      allocate( par_lst_tmp(1:2,1:n_orbs) )
      i3 = 1_int32
      do i1 = 1_int32, n_at ; do i2 = 1_int32, o(i1)%n_orbs
            par_lst_tmp(1,i3) = o(i1)%orb(i2)%i_par
            par_lst_tmp(2,i3) = o(i1)%orb(i2)%n_par
            i3 = i3 + 1_int32
         enddo ; enddo
      n_indp_par = 0_int32
      o_sym_par%n_indp = o_sym%n_indp
      i4 = 1_int32
      do i1 = 1_int32, o_sym%n_indp
         n_indp_par = n_indp_par + par_lst_tmp(2,o_sym%indp(i1)%cnct(1))
         if ( par_lst_tmp(2,o_sym%indp(i1)%cnct(1)).gt.0_int32) then
            o_sym_par%indp(i4)%n_cnct = o_sym%indp(i1)%n_cnct
            do i2 = 1_int32, o_sym%indp(i1)%n_cnct
               i3 = o_sym%indp(i1)%cnct(i2)
               o_sym_par%indp(i4)%cnct(2*i2+1) = par_lst_tmp(1,i3)
               o_sym_par%indp(i4)%cnct(2*i2+2) = par_lst_tmp(1,i3) + par_lst_tmp(2,i3) - 1_int32
            enddo
            i4 = i4 + 1_int32
         endif
      enddo
      i1 = 1_int32
      do while (i1.lt.o_sym_par%n_indp)
         check_pairs = .false.
         if ( o_sym_par%indp(i1)%n_cnct .eq. o_sym_par%indp(i1+1)%n_cnct ) then
            check_pairs = .true.
            i2loop : do i2 = 1_int32, o_sym_par%indp(i1)%n_cnct
               if (.not.((o_sym_par%indp(i1)%cnct(2*i2+2) + 1_int32) .eq. o_sym_par%indp(i1+1)%cnct(2*i2+1)))then
                  check_pairs = .false.
                  exit i2loop
               endif
            enddo i2loop
         endif
         if ( check_pairs ) then
            do i2 = 1_int32, o_sym_par%indp(i1)%n_cnct
               o_sym_par%indp(i1)%cnct(2*i2+2) = o_sym_par%indp(i1+1)%cnct(2*i2+2)
            enddo
            o_sym_par%indp(i1+1:o_sym_par%n_indp-1) = o_sym_par%indp(i1+2:o_sym_par%n_indp)
            o_sym_par%n_indp = o_sym_par%n_indp - 1_int32
         else
            i1 = i1 + 1_int32
         endif
      enddo
      o_sym_par%indp(1)%cnct(1:2) = o_sym_par%indp(1)%cnct(3:4)
      i1 = 2_int32
      do while (i1.le.o_sym_par%n_indp)
         if (o_sym_par%indp(i1)%n_cnct.gt.0_int32) then
            o_sym_par%indp(i1)%cnct(1) = o_sym_par%indp(i1-1)%cnct(2) + 1_int32
            o_sym_par%indp(i1)%cnct(2) = o_sym_par%indp(i1)%cnct(1) + o_sym_par%indp(i1)%cnct(4) &
            & - o_sym_par%indp(i1)%cnct(3)
            i1 = i1 + 1_int32
         else
            o_sym_par%n_indp = o_sym_par%n_indp - 1
         endif
      enddo
      i2 = par_lst_tmp(1,n_orbs) + par_lst_tmp(2,n_orbs) - 1_int32
      deallocate( par_lst_tmp )
      do i1 = 1_int32, o_sym_par%n_indp
         i2 = i2 - o_sym_par%indp(i1)%n_cnct*(o_sym_par%indp(i1)%cnct(2) - o_sym_par%indp(i1)%cnct(1) + 1 )
      enddo
      if (i2 .ne.0_int32) then
         call write_simple_line(stdout,0,mpi_rank,2,'l',"ERROR!!! Table of orbital parameter symmetries is wrong.")
         call write_simple_line(stdout,0,mpi_rank,11,'l',"Check fort.12 and contact matteo.barborini@gmail.com")
         if (mpi_rank.eq.0_int32) then
            write(12,*) o_sym_par%n_indp
            do i1 = 1_int32, o_sym_par%n_indp
               write(12,*) o_sym_par%indp(i1)%n_cnct, o_sym_par%indp(i1)%cnct(1:2*o_sym_par%indp(i1)%n_cnct+2)
            enddo
         endif
         stop
      endif
   end subroutine cmp_opar_sym
end module orbbss_cls
