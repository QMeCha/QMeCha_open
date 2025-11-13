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
module geminal_symmetries_m
   use fortran_kinds_v, only: int32, dp
   use orbbss_cls, only: orbspa_t, orb_t
   use molecular_system_v, only: n_at, d_nn
   use cubic_harmonics_m, only: y_l_sym_tab
   use merge_sort_m, only: mrg_srt_abs
   use molec_symmetry_m, only: sys_rot
   use molecular_system_m, only: n_sys_sym, sys_sym_tab
   use orbbss_cls, only: cmp_orbs_trsf
   implicit none
   type, public :: lmbd_ele_t
      integer(int32) :: o1
      integer(int32) :: o2
      integer(int32) :: io
   end type lmbd_ele_t
   type, public :: orb_lst_t
      integer(int32) :: at
      integer(int32) :: io
   end type orb_lst_t
   type, public :: lmbd_sym_t
      logical :: is_indp
      integer(int32) :: n_c
      integer(int32),  allocatable, dimension(:)  :: c
   end type lmbd_sym_t
contains
   subroutine cmp_nzr_fee( cpl, n_orbs, o, n_nn_ele, nn_ele )
      character(1),                            intent(in)    :: cpl
      integer(int32),                          intent(in)    :: n_orbs
      type(orbspa_t),  dimension(n_at),      intent(in)    :: o
      integer(int32),                          intent(inout) :: n_nn_ele
      type(lmbd_ele_t), dimension(n_nn_ele), intent(inout) :: nn_ele
      logical        :: ovrlp
      integer(int32) :: ia1, ia2
      integer(int32) :: io1, io2
      integer(int32) :: o1, o2
      n_nn_ele = 0_int32
      io2 = 0_int32
      do ia2 = 1_int32, n_at ; do o2 = 1_int32, o(ia2)%n_orbs
            io2 = io2 + 1_int32
            io1 = 0_int32
            do ia1 = 1_int32, n_at ; do o1 = 1_int32, o(ia1)%n_orbs
                  io1 = io1 + 1_int32
                  ovrlp = .false.
                  if ( cpl.eq.'a') then
                     if ( ia1 .eq. ia2 ) ovrlp = .true.
                  elseif ( cpl.eq.'d') then
                     if ( (io1.eq.io2) ) ovrlp = .true.
                  else
                     ovrlp = .true.
                  endif
                  if ( ovrlp ) then
                     n_nn_ele = n_nn_ele + 1_int32
                     nn_ele(n_nn_ele)%o1 = io1
                     nn_ele(n_nn_ele)%o2 = io2
                     nn_ele(n_nn_ele)%io = io1 + n_orbs * (io2 - 1)
                  endif
               enddo ; enddo
         enddo ; enddo
   end subroutine cmp_nzr_fee
   subroutine cmp_sgm_sym( build_nz, spin, gem_typ, n_orbs, o, n_nn_ele, nn_ele, n_sym_tb, sym_tb  )
      logical,                                 intent(in) :: build_nz
      character(1),                            intent(in) :: spin
      character(1),                            intent(in) :: gem_typ
      integer(int32),                          intent(in) :: n_orbs
      type(orbspa_t), dimension(n_at),       intent(in) :: o
      integer(int32),                          intent(inout) :: n_nn_ele
      type(lmbd_ele_t), dimension(n_nn_ele), intent(inout) :: nn_ele
      integer(int32),                          intent(inout) :: n_sym_tb
      type(lmbd_sym_t), dimension(n_nn_ele), intent(inout) :: sym_tb
      logical        :: acc1, acc2
      integer(int32) :: shft1, shft2
      integer(int32) :: sgn1, sgn2, sgn
      integer(int32) :: i1, i2, i3
      integer(int32) :: ia1, ia2
      integer(int32) :: io1, io2, io
      integer(int32) :: fo1, fo2, fo
      integer(int32) :: is
      if (spin.eq.'R' .or. gem_typ.eq.'t') then
         do i1 = 1_int32, n_nn_ele
            if( nn_ele(i1)%o1.lt.nn_ele(i1)%o2 ) then
               sym_tb(i1)%is_indp = .false.
               do i2 = 1_int32, n_nn_ele
                  if( (nn_ele(i2)%o1.eq.nn_ele(i1)%o2).and.(nn_ele(i2)%o2.eq.nn_ele(i1)%o1) ) then
                     sym_tb(i2)%n_c  = sym_tb(i2)%n_c + 1_int32
                     if ( gem_typ.eq.'t') then
                        sym_tb(i2)%c(2) = - sym_tb(i1)%c(1)
                     else
                        sym_tb(i2)%c(2) =   sym_tb(i1)%c(1)
                     endif
                  endif
               enddo
            endif
         enddo
      endif
      if(sys_rot) then
         do i1 = 2_int32, n_sys_sym
            is = sys_sym_tab(i1,0)
            do ia1 = 1_int32, n_at ; do i2 = 1_int32, o(sys_sym_tab(i1,ia1))%n_orbs
                  io1 = o(sys_sym_tab(i1,ia1))%orb(i2)%i_orb
                  call cmp_orbs_trsf( is, o(sys_sym_tab(i1,ia1))%orb(i2), acc1, sgn1, shft1 )
                  fo1 = o(sys_sym_tab(1,ia1))%orb(i2)%i_orb + shft1
                  if ( acc1 ) then
                     do ia2 = 1_int32, n_at ; do i3 = 1_int32, o(sys_sym_tab(i1,ia2))%n_orbs
                           io2 = o(sys_sym_tab(i1,ia2))%orb(i3)%i_orb
                           call cmp_orbs_trsf( is, o(sys_sym_tab(i1,ia2))%orb(i3), acc2, sgn2, shft2 )
                           fo2 = o(sys_sym_tab(1,ia2))%orb(i3)%i_orb + shft2
                           if ( acc2 ) then
                              io = io2 + n_orbs * (io1 - 1 )
                              fo = fo2 + n_orbs * (fo1 - 1 )
                              call conv_indx( n_nn_ele, nn_ele, io )
                              call conv_indx( n_nn_ele, nn_ele, fo )
                              if ( (io.gt.0).and.(fo.gt.0) ) then
                                 sgn = sgn1 * sgn2
                                 call fll_cpl_pos( gem_typ, n_nn_ele, io, fo, sgn, sym_tb  )
                              endif
                           endif
                        enddo ; enddo
                  endif
               enddo ; enddo
         enddo
      endif
      if ( gem_typ.eq.'s') then
         do i1 = 1_int32, n_nn_ele
            if( sym_tb(i1)%is_indp ) then
               i2_lp1: do i2 = 2_int32, sym_tb(i1)%n_c
                  if ( sym_tb(i1)%c(i2) .eq. -sym_tb(i1)%c(1) ) then
                     sym_tb(i1)%is_indp = .false.
                     exit i2_lp1
                  endif
               enddo i2_lp1
            endif
         enddo
      else if ( gem_typ.eq.'t') then
         do i1 = 1_int32, n_nn_ele
            if( sym_tb(i1)%is_indp ) then
               if ( nn_ele(abs(sym_tb(i1)%c(1)))%o1 .eq. nn_ele(abs(sym_tb(i1)%c(1)))%o2 ) then
                  sym_tb(i1)%is_indp = .false.
               endif
            endif
         enddo
      endif
      n_sym_tb = 0_int32
      i1 = 1_int32
      do while (i1.le.n_nn_ele)
         if(.not.sym_tb(i1)%is_indp) then
            i2_lp2: do i2 = i1+1, n_nn_ele
               if(sym_tb(i2)%is_indp) then
                  sym_tb(i1)%is_indp = sym_tb(i2)%is_indp
                  sym_tb(i1)%n_c = sym_tb(i2)%n_c
                  sym_tb(i1)%c = sym_tb(i2)%c
                  sym_tb(i2)%is_indp = .false.
                  exit i2_lp2
               endif
            enddo i2_lp2
         endif
         if(sym_tb(i1)%is_indp) then
            call mrg_srt_abs ( sym_tb(i1)%n_c, sym_tb(i1)%c(1:sym_tb(i1)%n_c))
            sym_tb(i1)%c(1:sym_tb(i1)%n_c) = sign(1_int32, sym_tb(i1)%c(1) ) * sym_tb(i1)%c(1:sym_tb(i1)%n_c)
            n_sym_tb = n_sym_tb + 1_int32
         endif
         i1 = i1 + 1_int32
      enddo
      if (build_nz) then
         do i1 = 1_int32, n_sym_tb
            do i2 = 1_int32, sym_tb(i1)%n_c
               nn_ele(abs(sym_tb(i1)%c(i2)))%io = 0_int32
            enddo
         enddo
         i2 = 0_int32
         do i1 = 1_int32, n_nn_ele
            if (nn_ele(i1)%io .eq. 0_int32) then
               i2 = i2 + 1_int32
               nn_ele(i1)%io = i2
            else
               nn_ele(i1)%io = 0_int32
            endif
         enddo
         do i1 = 1_int32, n_sym_tb
            do i2 = 1_int32, sym_tb(i1)%n_c
               sgn = sign(1,sym_tb(i1)%c(i2))
               sym_tb(i1)%c(i2) = sgn * nn_ele(abs(sym_tb(i1)%c(i2)))%io
            enddo
         enddo
         i1 = 1_int32
         do while (i1.le.n_nn_ele)
            if( nn_ele(i1)%io.eq.0_int32 ) then
               i2_lp3: do i2 = i1+1, n_nn_ele
                  if(nn_ele(i2)%io.ne.0_int32) then
                     nn_ele(i1)%io = nn_ele(i2)%io
                     nn_ele(i1)%o1 = nn_ele(i2)%o1
                     nn_ele(i1)%o2 = nn_ele(i2)%o2
                     nn_ele(i2)%io = 0_int32
                     exit i2_lp3
                  endif
               enddo i2_lp3
            endif
            i1 = i1 + 1_int32
         enddo
         i2 = 0_int32
         do i1 = 1_int32, n_nn_ele
            if(nn_ele(i1)%io.ne.0_int32) then
               i2 = i2 + 1_int32
               nn_ele(i1)%io = nn_ele(i1)%o1 + n_orbs * (nn_ele(i1)%o2 - 1)
            endif
         enddo
         n_nn_ele = i2
      endif
   end subroutine cmp_sgm_sym
   subroutine conv_indx( n_ele, nn_ele, io )
      integer(int32),                     intent(in)    :: n_ele
      type(lmbd_ele_t), dimension(n_ele), intent(in)    :: nn_ele
      integer(int32),                     intent(inout) :: io
      integer(int32) :: i1
      logical        :: is_found
      is_found = .false.
      do i1 = 1_int32, n_ele
         if ( io.eq.nn_ele(i1)%io) then
            is_found = .true.
            io = i1
            exit
         endif
      enddo
      if (.not.is_found) io = 0_int32
   end subroutine conv_indx
   subroutine fll_cpl_pos( gem_typ, n_ele, io, fo, sgn, sym_tb )
      character(1),                            intent(in) :: gem_typ
      integer(int32),                       intent(in)    :: n_ele
      integer(int32),                       intent(inout) :: io
      integer(int32),                       intent(inout) :: fo
      integer(int32),                       intent(inout) :: sgn
      type(lmbd_sym_t), dimension(n_ele), intent(inout) :: sym_tb
      integer(int32) :: i1, i2
      integer(int32) :: ro
      if ( .not.sym_tb(fo)%is_indp ) then
         i2 = 1_int32
         i1_lp1: do i1 = 1_int32, n_ele
            if ( sym_tb(i1)%is_indp ) then
               do i2 = 2_int32, sym_tb(i1)%n_c
                  ro = abs ( sym_tb(i1)%c(i2) )
                  if ( ro .eq. fo) exit i1_lp1
               enddo
            endif
         enddo i1_lp1
         fo = i1
         sgn = sgn * sign( 1_int32, sym_tb(fo)%c(i2) )
      endif
      if ( (io.ne.fo) .and. (.not.sym_tb(io)%is_indp)  ) then
         i2 = 1_int32
         i1_lp2: do i1 = 1_int32, n_ele
            if ( sym_tb(i1)%is_indp ) then
               do i2 = 2_int32, sym_tb(i1)%n_c
                  ro = abs ( sym_tb(i1)%c(i2) )
                  if ( ro .eq. io) exit i1_lp2
               enddo
            endif
         enddo i1_lp2
         io = i1
         sgn = sgn * sign( 1_int32, sym_tb(io)%c(i2) )
      endif
      if (io.lt.fo) then
         i1 = sym_tb(io)%n_c
         i2 = sym_tb(fo)%n_c
         sym_tb(io)%c(i1+1:i1+i2) = sgn * sym_tb(fo)%c(1:i2)
         sym_tb(io)%n_c = sym_tb(io)%n_c + i2
         sym_tb(fo)%is_indp = .false.
      elseif(fo.lt.io) then
         i1 = sym_tb(fo)%n_c
         i2 = sym_tb(io)%n_c
         sym_tb(fo)%c(i1+1:i1+i2) = sgn * sym_tb(io)%c(1:i2)
         sym_tb(fo)%n_c = sym_tb(fo)%n_c + i2
         sym_tb(io)%is_indp = .false.
      elseif((fo.eq.io).and.(gem_typ.eq.'s')) then
         if (sgn.lt.0) then
            ro = 0
            do i1 = 2_int32, sym_tb(io)%n_c
               if ( sym_tb(io)%c(i1) .eq. -io ) then
                  ro = 1
                  exit
               endif
            enddo
            if ( ro .eq. 0 ) then
               sym_tb(io)%n_c = sym_tb(io)%n_c + 2
               sym_tb(io)%c(sym_tb(io)%n_c-1:sym_tb(io)%n_c) = -sym_tb(io)%c(1:2)
            endif
         endif
      endif
   end subroutine fll_cpl_pos
   subroutine rmv_cpl_pos( n_ele, io, sym_tb )
      integer(int32),                       intent(in)    :: n_ele
      integer(int32),                       intent(inout) :: io
      type(lmbd_sym_t), dimension(n_ele), intent(inout) :: sym_tb
      sym_tb(io)%c(:) = 0
      sym_tb(io)%n_c = 0
      sym_tb(io)%is_indp = .false.
   end subroutine rmv_cpl_pos
   subroutine cmp_orbs_ovrlp( dr, o1, o2, ovr )
      real(dp),       dimension(3),  intent(in)  :: dr
      type(orb_t),                     intent(in)  :: o1, o2
      logical,                         intent(out) :: ovr
      real(dp),       dimension(3)               :: vec
      integer(int32)                               :: i1
      do i1 = 1_int32, 3_int32
         if ( abs(dr(i1)).gt.0.0_dp ) then
            vec(i1) = 0.0_dp
         else
            vec(i1) = abs(y_l_sym_tab(i1,0,o1%l_z)) + abs(y_l_sym_tab(i1,0,o2%l_z))
            vec(i1) = sign( vec(i1), y_l_sym_tab(i1,0,o1%l_z) * y_l_sym_tab(i1,0,o2%l_z) )
         endif
      enddo
      ovr = .false.
      if ( abs( mod(vec(1),2.0_dp)+mod(vec(2),2.0_dp)+mod(vec(3),2.0_dp) ).lt.10d-13 ) then
         if ( abs(vec(1)).lt.10d-13 .or. abs(vec(1)+vec(2)).gt.10d-13 ) ovr = .true.
      endif
   end subroutine cmp_orbs_ovrlp
   subroutine ord_orbs( n_orbs, a1, a2, aout )
      integer(int32), intent(in)    :: n_orbs
      integer(int32), intent(in)    :: a1
      integer(int32), intent(in)    :: a2
      integer(int32), intent(out)   :: aout
      if ( a1 .gt. a2 ) then
         aout = a1 + n_orbs * (a2 -1 )
      else
         aout = a2 + n_orbs * (a1 -1 )
      endif
   end subroutine ord_orbs
   subroutine iord_orbs( a1, a2 )
      integer(int32), intent(inout) :: a1
      integer(int32), intent(inout) :: a2
      integer(int32)                :: a3
      if ( a1 .lt. a2 ) then
         a3 = a1 ; a1 = a2 ; a2 = a3
      endif
   end subroutine iord_orbs
   subroutine ord_cpls( n_orbs, a1, a2, b1, b2, a, b )
      integer(int32), intent(in)    :: n_orbs
      integer(int32), intent(inout) :: a1
      integer(int32), intent(inout) :: a2
      integer(int32), intent(inout) :: b1
      integer(int32), intent(inout) :: b2
      integer(int32), intent(out)   :: a
      integer(int32), intent(out)   :: b
      integer(int32)                 :: c1, c2
      a = a1 + n_orbs * (a2 -1 )
      b = b1 + n_orbs * (b2 -1 )
      if ( a .gt. b ) then
         c1 = a1 ; c2 = a2
         a1 = b1 ; a2 = b2
         b1 = c1 ; b2 = c2
         a = a1 + n_orbs * (a2 -1 )
         b = b1 + n_orbs * (b2 -1 )
      endif
   end subroutine ord_cpls
end module geminal_symmetries_m
