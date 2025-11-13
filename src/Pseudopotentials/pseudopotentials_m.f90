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
module pseudopotentials_m
   use fortran_kinds_v,          only: dp, int32, stdout
   use openmp_mpi_m
   use write_lines_m,            only: write_separator_line, &
   & write_simple_line, write_empty_line, write_variable_line
   use molecular_system_v,       only: n_at, n_fe, atoms, f_chrg
   use fermionic_config_c,       only: frm_cnf
   use fermionic_wavefunction_m, only: wvfn
#ifdef _EQJ
   use drudonic_config_c,        only: drd_cnf
   use qdo_wavefunction_m,       only: wvfnq
   use qdo_system_v,             only: qdosys_prs
#endif
   use pseudopotentials_v
   use nl_pseudo_grid_m
   implicit none
   public  :: init_pseudopotentials, eval_pseudo_cmp, eval_pseudo_ene
   private :: link_psdpot_atoms, eval_pseudo_cutoff
contains
   subroutine init_pseudopotentials( )
      if (psd_prs) then
         call write_separator_line(stdout,0,mpi_rank,2,"_")
         call write_simple_line(stdout,0,mpi_rank,2,"c","INITIALIZING  PSEUDOPOTENTIALS")
         call write_empty_line(stdout,0,mpi_rank)
         call write_variable_line(stdout,0,mpi_rank,2,"Energy cut-off for pseudopotentials", ene_psd_cut,var_name="ene_psd_cut")
         if (psd_int_mthd.eq.0_int32) then
            call write_variable_line(stdout,0,mpi_rank,2,"Locality approximation (LA)", psd_int_mthd,var_name="psd_int_mthd")
         else
            call write_variable_line(stdout,0,mpi_rank,2,"Determinant locality approximation (DLA)", psd_int_mthd,var_name="psd_int_mthd")
         endif
         call link_psdpot_atoms()
         call eval_pseudo_cutoff()
         call init_nloc_grid()
         call prnt_psdpot()
      endif
   end subroutine init_pseudopotentials
   subroutine link_psdpot_atoms()
      character(2)               :: red_atm_name = ""
      integer(int32)             :: i1, i2, ierr
      logical,          external :: is_numeric
      logical                    :: pseudo_is_present
      ierr = 0
      do i1 = 1_int32, n_at
         if ( atoms(i1)%i_pp.lt.0_int32 ) then
            pseudo_is_present = .false.
            if ( is_numeric( atoms(i1)%atm_name(3:3) ) ) then
               red_atm_name(1:1) = atoms(i1)%atm_name(2:2)
            else
               red_atm_name(1:2) = atoms(i1)%atm_name(2:3)
            endif
            do i2 = 1_int32, n_psdpot
               if ( trim(red_atm_name(1:2)).eq.trim(psdpot(i2)%atm_name(1:2)) ) then
                  atoms(i1)%i_pp = i2
                  atoms(i1)%atm_z = atoms(i1)%atm_z - int( psdpot(i2)%z_rep )
                  pseudo_is_present = .true.
                  exit
               endif
            enddo
            if (.not.pseudo_is_present) then
               ierr = -1
               exit
            endif
         endif 
      enddo 
      if ( ierr .ne. 0_int32 ) then
         call write_variable_line(stdout,0,mpi_rank,2,"ERROR!!! Pseudopotential need by pseudo atom is missing", red_atm_name)
         call fnlz_ompmpi_env()
         stop
      endif
   end subroutine link_psdpot_atoms
   subroutine eval_pseudo_cutoff( )
      real(dp)       :: psd_val(0:5), d_af
      integer(int32) :: i1
      do i1 = 1_int32, n_psdpot
         psdpot(i1)%r_cut = 0.0_dp
         d_af = 0.5_dp
         psd_val = 10.0_dp
         do while ( sum(abs(psd_val(0:5))).gt.ene_psd_cut )
            call eval_pseudo_cmp( d_af, psdpot(i1), psd_val(0:5) )
            if ( sum(abs(psd_val(0:5))).le.ene_psd_cut ) then
               psdpot(i1)%r_cut = d_af
            else
               d_af = d_af + 0.01_dp
            endif
         enddo
      enddo ! i1 (n_psdpot)
   end subroutine eval_pseudo_cutoff
   subroutine eval_pseudo_cmp( d_af, pseudo, psd_pot)
      real(dp),                 intent(in)    :: d_af
      type(psdpot_t),           intent(inout) :: pseudo
      real(dp), dimension(0:5), intent(inout) :: psd_pot
      integer(int32)  :: i1, i2, i3
      psd_pot = 0.0_dp
      i3 = 0_int32
      do i1 = 1_int32, pseudo%n_psd_c
         do i2 = 1_int32, pseudo%n_psd_ele(i1)
            i3 = i3 + 1_int32
            psd_pot(i1-1) = psd_pot(i1-1) + pseudo%c(i3) * d_af**(pseudo%n(i3)-2_int32) * dexp(-pseudo%e(i3)*d_af**2 )
         enddo
      enddo
   end subroutine eval_pseudo_cmp
   subroutine eval_local_cmp( d_af, pseudo, psd_loc)
      real(dp),                 intent(in)    :: d_af
      type(psdpot_t),           intent(inout) :: pseudo
      real(dp),                 intent(inout) :: psd_loc
      real(dp)                                :: zr2
      integer(int32) :: i1
      do i1 = 1_int32, pseudo%n_psd_ele(1)
         zr2 =  pseudo%e(i1)* d_af**2
         if ( zr2.lt.300_dp ) then
            psd_loc = psd_loc + pseudo%c(i1) * d_af**(pseudo%n(i1)-2_int32) * dexp( -zr2 )
         endif
      enddo
   end subroutine eval_local_cmp
   subroutine eval_pseudo_ene( iw, v_pp )
      integer(int32),             intent(in)    :: iw
      real(dp),  dimension(n_fe), intent(inout) :: v_pp
      real(dp), dimension(0:5) :: psd_pot
      integer(int32) :: i1, i2
      v_pp = 0.0_dp
      do i1 = 1_int32, n_at
         if ( atoms(i1)%i_pp.gt.0_int32 ) then
            do i2 = 1_int32, n_fe
               if ( frm_cnf(iw)%d_fn(i1,i2)%m.le.psdpot(atoms(i1)%i_pp)%r_cut ) then
                  call eval_pseudo_cmp(frm_cnf(iw)%d_fn(i1,i2)%m, psdpot(atoms(i1)%i_pp), psd_pot(0:5) )
                  if ( f_chrg(i2).lt.0_int32 .and. psdpot(atoms(i1)%i_pp)%n_psd_c.ge.2_int32 &
                  &.and. sum(abs(psd_pot(1:psdpot(atoms(i1)%i_pp)%n_psd_c))).gt.0.0_dp ) then
                     frm_cnf(iw)%i_fe = i2
                     call rndm_nlc_grd( iw, i1, i2 )
                     call eval_nlc_intgr( iw, psdpot(atoms(i1)%i_pp)%n_psd_c-1_int32, psd_pot(0:5) )
                  endif
                  if ( f_chrg(i2).lt.0_int32 ) then
                     v_pp(i2) = v_pp(i2) + sum( psd_pot(0:5) )
                  else
                     v_pp(i2) = v_pp(i2) - psd_pot(0)
                  endif
               endif
            enddo 
         endif
      enddo 
   end subroutine eval_pseudo_ene
   subroutine eval_nlc_intgr( iw, n_nlc_c, psd_pot )
      integer(int32),           intent(in)    :: iw
      integer(int32),           intent(in)    :: n_nlc_c
      real(dp), dimension(0:5), intent(inout) :: psd_pot
      real(dp), dimension(n_nlc_c) :: psd_int
      real(dp)                       :: g
#ifdef _EQJ
      real(dp)                       :: g_qe
#endif
      integer(int32)  :: i1, i2
      psd_int = 0.0_dp
      do i1 = 1_int32, n_psd_grd_pnts
         frm_cnf(iw)%r_fe_new = nlc_rnd_grid(:,i1,iw)
         call frm_cnf(iw)%new(  )
#ifdef _EQJ
         if( qdosys_prs ) call drd_cnf(iw)%new( iw )
#endif
         if ( psd_int_mthd.eq.0_int32 ) then
            call wvfn%ratla( iw, g, pseudo_update=.true. )
#ifdef _EQJ
            if( qdosys_prs ) then
               call wvfnq(iw)%ratio( iw, g_qe )
               g = g * g_qe
            endif
#endif
         else
            call wvfn%ratde( iw, g, pseudo_update=.true. )
         endif
         g = g * w_grid(i1)
         do i2 = 1_int32, n_nlc_c
            psd_int(i2) = psd_int(i2) + g * P_l(i2-1_int32,nlc_cos_theta(i1,iw))
         enddo 
      enddo 
      do i1 = 1_int32, n_nlc_c
         psd_pot(i1) = psd_pot(i1) * psd_int(i1) * dble(2*(i1-1)+1) ! / (2.0_dp * twopi)
      enddo
   end subroutine eval_nlc_intgr
end module pseudopotentials_m
