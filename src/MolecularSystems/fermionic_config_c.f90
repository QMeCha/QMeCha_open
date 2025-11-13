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
module fermionic_config_c
   use fortran_kinds_v,         only: int32, dp
   use mersenne_twister19937_m, only: rndn
   use quantum_monte_carlo_v,   only: spc_map, mapping_type
   use molecular_system_v,      only: dist_t, n_fe, n_el, n_po, n_at, r_at, &
   & atoms, nuc_chrg
   use fermionic_positions_m,   only: frmdst
   implicit none
   type, public :: fermionic_config_t
      integer(int32)                            :: i_fe
      integer(int32)                            :: i_tmv
      integer(int32), allocatable, dimension(:) :: i_fe_tmv
      integer(int32)                            :: n_fe_tmv
      real(dp),     allocatable, dimension(:,:) :: r_fe
      type(dist_t), allocatable, dimension(:)   :: d_ee, d_pp
      type(dist_t), allocatable, dimension(:,:) :: d_ep
      type(dist_t), allocatable, dimension(:,:) :: d_fn
      real(dp),                  dimension(3) :: r_fe_new
      type(dist_t), allocatable, dimension(:) :: d_ff_new
      type(dist_t), allocatable, dimension(:) :: d_fn_new
      real(dp),       allocatable, dimension(:) :: sm_el
      real(dp)                                  :: sm_el_new
      integer(int32), allocatable, dimension(:) :: fe_i_at
      integer(int32)                            :: fe_i_at_new
   contains
      procedure :: ini => init_fermionic_conf
      procedure :: cmp => comp_fermionic_conf
      procedure :: chs => chose_fermion_to_move
      procedure :: new => new_fermionic_conf
      procedure :: upd => update_fermionic_conf
   end type fermionic_config_t
   type(fermionic_config_t), save, public, allocatable, dimension(:), target :: frm_cnf
   private :: init_fermionic_conf, new_fermionic_conf, comp_fermionic_conf, &
   & update_fermionic_conf, chose_fermion_to_move
   external :: space_map
contains
   subroutine init_fermionic_conf( frm_cnf_in, iw )
      class(fermionic_config_t), intent(inout) :: frm_cnf_in
      integer(int32),            intent(in)    :: iw
      integer(int32) :: i1, i2, i3
      allocate( frm_cnf_in%r_fe(1:3,1:n_fe) )
      if (spc_map) then
         allocate( frm_cnf_in%sm_el(1:n_fe) ) ; frm_cnf_in%sm_el(1:n_fe) = 1.0_dp
         allocate( frm_cnf_in%fe_i_at(1:n_fe) ) ; frm_cnf_in%fe_i_at(1:n_fe) = 0_int32
      endif
      if ( n_at .gt.0_int32 ) then
         allocate( frm_cnf_in%d_fn(1:n_at,1:n_fe) )
         allocate( frm_cnf_in%d_fn_new(1:n_at) )
      endif
      if ( n_el.ge.2_int32 ) then
         allocate( frm_cnf_in%d_ee(1:n_el*(n_el-1)/2) )
      endif
      if ( n_po.ge.2_int32 ) then
         allocate( frm_cnf_in%d_pp(1:n_po*(n_po-1)/2) )
      endif
      if ( n_po.ge.1_int32 .and. n_el .ge. 1_int32) then
         allocate( frm_cnf_in%d_ep(1:n_po,1:n_el) )
      endif
      allocate( frm_cnf_in%d_ff_new(1:n_fe) )
      allocate( frm_cnf_in%i_fe_tmv(1:n_fe) )
      i3 = 0_int32
      if ( n_at.gt.0_int32 ) then
         if ( sum(frmdst).eq.n_fe) then
            do i1 = 1_int32, n_at
               do i2 = 1_int32, frmdst(1,i1)
                  i3 = i3 + 1_int32
                  frm_cnf_in%r_fe(1,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(1,i1)
                  frm_cnf_in%r_fe(2,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(2,i1)
                  frm_cnf_in%r_fe(3,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(3,i1)
               enddo
            enddo
            do i1 = 1_int32, n_at
               do i2 = 1_int32, frmdst(2,i1)
                  i3 = i3 + 1_int32
                  frm_cnf_in%r_fe(1,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(1,i1)
                  frm_cnf_in%r_fe(2,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(2,i1)
                  frm_cnf_in%r_fe(3,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(3,i1)
               enddo
            enddo
            if (n_po.gt.0_int32) then
               do i1 = 1_int32, n_at
                  do i2 = 1_int32, frmdst(3,i1)
                     i3 = i3 + 1_int32
                     frm_cnf_in%r_fe(1,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(1,i1)
                     frm_cnf_in%r_fe(2,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(2,i1)
                     frm_cnf_in%r_fe(3,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(3,i1)
                  enddo
               enddo
               do i1 = 1_int32, n_at
                  do i2 = 1_int32, frmdst(4,i1)
                     i3 = i3 + 1_int32
                     frm_cnf_in%r_fe(1,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(1,i1)
                     frm_cnf_in%r_fe(2,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(2,i1)
                     frm_cnf_in%r_fe(3,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(3,i1)
                  enddo
               enddo
            endif
         else
            if (n_el.gt.0_int32) then
               do i1 = 1_int32, n_at ; do i2 = 1_int32, abs(atoms(i1)%atm_z)
                     if( i3.lt.n_el ) then
                        i3 = i3 + 1_int32
                        frm_cnf_in%r_fe(1,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(1,i1)
                        frm_cnf_in%r_fe(2,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(2,i1)
                        frm_cnf_in%r_fe(3,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo() + r_at(3,i1)
                     endif
                  enddo ; enddo ! i1 (n_at) i2 (atoms(i1)%atm_z)
               if(i3.lt.n_el) then
                  do while (i3.lt.n_el)
                     i3 = i3 + 1_int32
                     frm_cnf_in%r_fe(1,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo()
                     frm_cnf_in%r_fe(2,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo()
                     frm_cnf_in%r_fe(3,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo()
                  enddo
               endif
            endif
            if (n_po.gt.0_int32) then
               if(i3.lt.n_fe) then
                  do while (i3.lt.n_fe)
                     i3 = i3 + 1_int32
                     frm_cnf_in%r_fe(1,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo()
                     frm_cnf_in%r_fe(2,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo()
                     frm_cnf_in%r_fe(3,i3) = 1.0_dp - 2.0_dp * rndn(iw)%rndo()
                  enddo
               endif
            endif
         endif
      else
         do i2 = 1_int32, n_fe
            frm_cnf_in%r_fe(1,i2) = 1.0_dp - 2.0_dp * rndn(iw)%rndo()
            frm_cnf_in%r_fe(2,i2) = 1.0_dp - 2.0_dp * rndn(iw)%rndo()
            frm_cnf_in%r_fe(3,i2) = 1.0_dp - 2.0_dp * rndn(iw)%rndo()
         enddo
      endif
      frm_cnf_in%n_fe_tmv = n_fe
      do i1 = 1_int32, n_fe
         frm_cnf_in%i_fe_tmv(i1) = i1
      enddo ! i1 ( n_el )
      frm_cnf_in%i_tmv = 1_int32
      call frm_cnf_in%cmp()
   end subroutine init_fermionic_conf
   subroutine comp_fermionic_conf( frm_cnf_in )
      class(fermionic_config_t), intent(inout) :: frm_cnf_in
      integer(int32) :: i1, i2, i3
      if( n_el.ge.2_int32 ) then
         i3 = 0_int32
         do i2 = 1_int32, n_el ; do i1 = i2 + 1_int32, n_el
               i3 = i3 + 1_int32
               frm_cnf_in%d_ee(i3)%v(:) = frm_cnf_in%r_fe(:,i2) - frm_cnf_in%r_fe(:,i1)
               frm_cnf_in%d_ee(i3)%m    = dsqrt(sum(frm_cnf_in%d_ee(i3)%v(:)**2))
            enddo ; enddo ! i1 (n_el) i2 (n_el)
      endif
      if( n_po.ge.2_int32 ) then
         i3 = 0_int32
         do i2 = 1_int32, n_po ; do i1 = i2 + 1_int32, n_po
               i3 = i3 + 1_int32
               frm_cnf_in%d_pp(i3)%v(:) = frm_cnf_in%r_fe(:,i2+n_el) - frm_cnf_in%r_fe(:,i1+n_el)
               frm_cnf_in%d_pp(i3)%m    = dsqrt(sum(frm_cnf_in%d_pp(i3)%v(:)**2))
            enddo ; enddo ! i1 (n_po) i2 (n_po)
      endif
      if( n_po.ge.1_int32 .and.n_el.ge.1_int32) then
         do i2 = 1_int32, n_el ; do i1 = 1_int32, n_po
               frm_cnf_in%d_ep(i1,i2)%v(:) = frm_cnf_in%r_fe(:,i2) - frm_cnf_in%r_fe(:,i1+n_el)
               frm_cnf_in%d_ep(i1,i2)%m    = dsqrt(sum(frm_cnf_in%d_ep(i1,i2)%v(:)**2))
            enddo ; enddo ! i1 (n_el) i2 (n_po)
      endif
      if (n_at.gt.0_int32) then
         do i1 = 1_int32, n_fe ; do i2 = 1_int32, n_at
               frm_cnf_in%d_fn(i2,i1)%v(:) = frm_cnf_in%r_fe(:,i1) - r_at(:,i2)
               frm_cnf_in%d_fn(i2,i1)%m    = dsqrt(sum(frm_cnf_in%d_fn(i2,i1)%v(:)**2))
            enddo ; enddo ! i1 (n_el) i2 (n_at)
      endif
      if (spc_map) then
         do i1 = 1_int32, n_fe
            call space_map(frm_cnf_in%d_fn(1:n_at,i1),mapping_type, frm_cnf_in%fe_i_at(i1), frm_cnf_in%sm_el(i1) )
         enddo ! i1 (n_el)
      endif
   end subroutine comp_fermionic_conf
   subroutine chose_fermion_to_move( frm_cnf_in, iw )
      class(fermionic_config_t), intent(inout) :: frm_cnf_in
      integer(int32),            intent(in)    :: iw
      if( frm_cnf_in%n_fe_tmv.eq.1_int32) then
         frm_cnf_in%i_tmv = 1_int32
      else
         frm_cnf_in%i_tmv = int( dble(frm_cnf_in%n_fe_tmv) * rndn(iw)%rndo() + 1.0_dp )
      endif
      frm_cnf_in%i_fe = frm_cnf_in%i_fe_tmv(frm_cnf_in%i_tmv)
      if( frm_cnf_in%n_fe_tmv.eq.1_int32) then
         frm_cnf_in%n_fe_tmv = n_fe
      else
         if(frm_cnf_in%i_tmv.ne.frm_cnf_in%n_fe_tmv) then
            frm_cnf_in%i_fe_tmv(frm_cnf_in%i_tmv)    = frm_cnf_in%i_fe_tmv(frm_cnf_in%n_fe_tmv)
            frm_cnf_in%i_fe_tmv(frm_cnf_in%n_fe_tmv) = frm_cnf_in%i_fe
         endif
         frm_cnf_in%n_fe_tmv = frm_cnf_in%n_fe_tmv - 1_int32
      endif
   end subroutine chose_fermion_to_move
   subroutine new_fermionic_conf( frm_cnf_in )
      class(fermionic_config_t), intent(inout) :: frm_cnf_in
      integer(int32) :: i1
      if (n_at.gt.0_int32) then
         do i1 = 1_int32, n_at
            frm_cnf_in%d_fn_new(i1)%v(:) = frm_cnf_in%r_fe_new(:) - r_at(:,i1)
            frm_cnf_in%d_fn_new(i1)%m    = dsqrt(sum(frm_cnf_in%d_fn_new(i1)%v(:)**2))
         enddo ! i1 (n_at)
      endif
      if( n_fe.ge.2_int32 ) then
         do i1 = 1_int32, n_fe
            if( frm_cnf_in%i_fe.ne.i1 ) then
               frm_cnf_in%d_ff_new(i1)%v(:) = frm_cnf_in%r_fe_new(:) - frm_cnf_in%r_fe(:,i1)
               frm_cnf_in%d_ff_new(i1)%m    = dsqrt(sum(frm_cnf_in%d_ff_new(i1)%v(:)**2))
            else
               frm_cnf_in%d_ff_new(i1)%v(:) = 0.0_dp
               frm_cnf_in%d_ff_new(i1)%m    = 0.0_dp
            endif
         enddo 
      endif
      if (spc_map) then
         call space_map( frm_cnf_in%d_fn_new(1:n_at), mapping_type, frm_cnf_in%fe_i_at_new, frm_cnf_in%sm_el_new )
         if (frm_cnf_in%i_fe.gt.n_el) frm_cnf_in%sm_el_new  = 1.0_dp
      endif
   end subroutine new_fermionic_conf
   subroutine update_fermionic_conf( frm_cnf_in )
      class(fermionic_config_t),  intent(inout) :: frm_cnf_in
      integer(int32) :: i1, i2, ie
      frm_cnf_in%r_fe(:,frm_cnf_in%i_fe) = frm_cnf_in%r_fe_new(:)
      if (n_at.gt.0_int32) then
         frm_cnf_in%d_fn(:,frm_cnf_in%i_fe) = frm_cnf_in%d_fn_new(:)
      endif
      if( frm_cnf_in%i_fe.le.n_el ) then
         if (n_el.ge.2_int32) then
            do i1 = 1, frm_cnf_in%i_fe-1
               i2 = n_el * (i1-1)-i1 * (i1+1) / 2 + frm_cnf_in%i_fe
               frm_cnf_in%d_ee(i2)%v(:) = - frm_cnf_in%d_ff_new(i1)%v(:)
               frm_cnf_in%d_ee(i2)%m    =   frm_cnf_in%d_ff_new(i1)%m
            enddo 
            i2 = n_el * (frm_cnf_in%i_fe-1)-frm_cnf_in%i_fe* (frm_cnf_in%i_fe+1) / 2 + frm_cnf_in%i_fe
            do i1 = frm_cnf_in%i_fe+1, n_el
               i2 = i2 + 1_int32
               frm_cnf_in%d_ee(i2)%v(:) = frm_cnf_in%d_ff_new(i1)%v(:)
               frm_cnf_in%d_ee(i2)%m    = frm_cnf_in%d_ff_new(i1)%m
            enddo 
         endif
         if ( n_po.ge.1_int32 ) then
            frm_cnf_in%d_ep(:,frm_cnf_in%i_fe) = frm_cnf_in%d_ff_new(n_el+1:n_fe)
         endif
      else
         ie = frm_cnf_in%i_fe - n_el
         if (n_el.ge.1_int32) then
            do i1 = 1_int32, n_el
               frm_cnf_in%d_ep(ie,i1)%v =-frm_cnf_in%d_ff_new(i1)%v
               frm_cnf_in%d_ep(ie,i1)%m = frm_cnf_in%d_ff_new(i1)%m
            enddo
         endif
         if (n_po.ge.2_int32) then
            do i1 = 1, ie-1
               i2 = n_po * (i1-1)-i1 * (i1+1) / 2 + ie
               frm_cnf_in%d_pp(i2)%v(:) = - frm_cnf_in%d_ff_new(i1+n_el)%v(:)
               frm_cnf_in%d_pp(i2)%m    =   frm_cnf_in%d_ff_new(i1+n_el)%m
            enddo ! i2 ( ie )
            i2 = n_po * (ie-1)-ie* (ie+1) / 2 + ie
            do i1 = ie+1, n_po
               i2 = i2 + 1_int32
               frm_cnf_in%d_pp(i2)%v(:) = frm_cnf_in%d_ff_new(i1+n_el)%v(:)
               frm_cnf_in%d_pp(i2)%m    = frm_cnf_in%d_ff_new(i1+n_el)%m
            enddo 
         endif
      endif
      if (spc_map) then
         frm_cnf_in%sm_el(frm_cnf_in%i_fe)   = frm_cnf_in%sm_el_new
         frm_cnf_in%fe_i_at(frm_cnf_in%i_fe) = frm_cnf_in%fe_i_at_new
      endif
   end subroutine update_fermionic_conf
end module fermionic_config_c
