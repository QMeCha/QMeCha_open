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
module wvfcnv_mod
   use fortran_kinds_v, only: dp, int32, stdout
   use write_lines_m
   use quantum_monte_carlo_v, only: restart
   implicit none
   character(3), save, public :: wf_name_cnv
   real(dp), save, public, allocatable, dimension(:,:) :: ovr_mat, tmp_mat
   integer(int32), save, public, allocatable, dimension(:)   :: cor_tab
   integer(int32), save, public :: wf_type_e_new, wf_type_p_new
   public  :: bcst_input_wvfcnv, eval_wvfcnv
contains
   subroutine bcst_input_wvfcnv( )
      use openmp_mpi_m
      use molecular_system_v, only: fle_coord
      use basisset_v, only: fle_basis
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast(fle_coord, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(fle_basis, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(restart,     1,   MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
      call mpi_bcast(wf_name_cnv, 5, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
#endif
   end subroutine bcst_input_wvfcnv
   subroutine eval_wvfcnv( )
      use fermionic_wavefunction_v, only: wf_type_e
      select case ( abs(wf_type_e) )
       case(1_int32)
         call sd_overlaps()
       case default
      end select
      select case ( abs(wf_type_e) )
       case(1_int32)
         select case ( abs(wf_type_e_new) )
          case(1_int32)
          case(2_int32)
            call sd2sg( )
         case(3_int32)
           call sd2tg( )
         case(4_int32)
           call sd2pf( )
          case default
         end select
       case default
      end select
   end subroutine eval_wvfcnv
   subroutine sd_overlaps( )
      use openmp_mpi_m, only: mpi_rank
      use molecular_system_v, only: n_el, n_el_u, n_el_d
      use slater_determinant_params_m, only: sld_c_e
      use fermionic_orbitals_m, only: n_forbs
      use fermionic_wavefunction_v, only: wf_type_e
      integer(int32) :: i1, i2
      real(dp)       :: norm
      call write_empty_line(stdout,0,mpi_rank)
      call write_separator_line(stdout,0,mpi_rank,2,"=")
      call write_simple_line(stdout,0,mpi_rank,2,"c","Computing Molecular orbital overlap")
      call write_empty_line(stdout,0,mpi_rank)
      allocate( ovr_mat(1:n_el_u,1:n_el_u) ) ; ovr_mat = 0.0_dp
      if (wf_type_e.lt.0_int32) then
         allocate( tmp_mat(1:n_forbs,1:n_el_u) ) ; tmp_mat = 0.0_dp
      endif
      if (n_el.gt.0_int32) then
         if (mpi_rank.eq.0_int32) &
         & write(*,'(3X,"Renormalizing molecular orbitals...                      ")')
         do i1 = 1_int32, n_el_u
            norm = dot_product( sld_c_e%mat_L_u(1:n_forbs,i1), sld_c_e%mat_L_u(1:n_forbs,i1) )
            sld_c_e%mat_L_u(1:n_forbs,i1) = sld_c_e%mat_L_u(1:n_forbs,i1) / sqrt(norm)
         enddo
         do i1 = 1_int32, n_el_d
            norm = dot_product( sld_c_e%mat_L_d(1:n_forbs,i1), sld_c_e%mat_L_d(1:n_forbs,i1) )
            sld_c_e%mat_L_d(1:n_forbs,i1) = sld_c_e%mat_L_d(1:n_forbs,i1) / sqrt(norm)
         enddo
         if (mpi_rank.eq.0_int32) call write_empty_line(stdout,0,mpi_rank)
         if (mpi_rank.eq.0_int32) &
         & write(*,'(3X,"Computing overlap matrix...                              ")')
         call dgemm('T','N', n_el_u, n_el_d, n_forbs, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,1:n_el_u), n_forbs, &
         & sld_c_e%mat_L_d(1:n_forbs,1:n_el_d), n_forbs, 0.0_dp, ovr_mat(1:n_el_u,1:n_el_d), n_el_u )
         ovr_mat = abs(ovr_mat)
         if (wf_type_e.lt.0_int32) then
            do i1 = 1_int32, n_el_d
               if ( sum(ovr_mat(1:n_el_u,i1)).lt.0.01_dp ) then
                  do i2 = 1_int32, n_el_u
                     if ( sum(ovr_mat(i2,1:n_el_d)).lt.0.01_dp ) then
                        ovr_mat(i2,i1) = 1_int32
                        exit
                     endif
                  enddo
               endif
            enddo
            if ( n_el_u.gt.n_el_d ) then
               i1 = n_el_d + 1
               do i2 = 1_int32, n_el_u
                  if ( sum(ovr_mat(i2,1:n_el_d)).lt.0.01_dp ) then
                     ovr_mat(i2,i1) = 1_int32
                     i1 = i1 + 1_int32
                  endif
               enddo
            endif
            call write_empty_line(stdout,0,mpi_rank)
            call dgemm('N','N', n_forbs, n_el_u, n_el_u, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,1:n_el_u), n_forbs, &
            & ovr_mat(1:n_el_u,1:n_el_u), n_el_u, 0.0_dp, tmp_mat(1:n_forbs,1:n_el_u), n_forbs )
            sld_c_e%mat_L_u(1:n_forbs,1:n_el_u) = tmp_mat(1:n_forbs,1:n_el_u)
            call dgemm('T','N', n_el_u, n_el_d, n_forbs, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,1:n_el_u), n_forbs, &
            & sld_c_e%mat_L_d(1:n_forbs,1:n_el_d), n_forbs, 0.0_dp, ovr_mat(1:n_el_u,1:n_el_d), n_el_u )
            ovr_mat = abs(ovr_mat)
         endif
         if (mpi_rank.eq.0_int32) then
            write(*,'(3X,"Found non zero overlaps ...                              ")')
            call write_empty_line(stdout,0,mpi_rank)
            write(*,'(3X,"                    Mol Num.         Overlap                  ")')
            write(*,'(3X,"                 Up.      Down.                                   ")')
            do i1 = 1_int32, n_el_u ; do i2 = 1_int32, n_el_d
                  if ( ovr_mat(i1,i2).gt.0.0_dp ) write(*,'(12X,2I10,F15.5)') i1, i2, ovr_mat(i1,i2)
               enddo ; enddo
            call write_empty_line(stdout,0,mpi_rank)
         endif
      endif
      if (allocated(tmp_mat)) deallocate(tmp_mat)
   end subroutine sd_overlaps
   subroutine sd2tg( )
      use openmp_mpi_m, only: mpi_rank
      use molecular_system_v, only: n_el, n_po, n_el_u, n_el_d, n_po_u, n_po_d
      use slater_determinant_params_m,               only: sld_c_e
      use triplet_geminal_params_m,               only: tgm_c_e
      use fermionic_orbitals_m,     only: n_forbs
      use fermionic_wavefunction_v, only: wf_type_e
      integer(int32) :: i1, i2
      integer(int32) :: mat_dim
      ovr_mat = 0.0_dp 
      do i1 = 1_int32, n_el_u-1
         ovr_mat(i1,i1+1_int32) = 1.0_dp
         ovr_mat(i1+1_int32,i1) = -1.0_dp
      enddo 
      allocate( tmp_mat(1:n_forbs,1:n_el_u) ) ; tmp_mat = 0.0_dp
      if (tgm_c_e%d_Zu.gt.n_forbs) then
         tgm_c_e%mat_Zu(1:n_forbs,tgm_c_e%d_Zu) = sld_c_e%mat_L_u(1:n_forbs,n_el_u)
         mat_dim = n_el_u-1
      else
         mat_dim = n_el_u
      endif
      call dgemm('N','N', n_forbs, mat_dim, mat_dim, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,1:mat_dim), n_forbs, &
      & ovr_mat(1:mat_dim,1:mat_dim), mat_dim, 0.0_dp, tmp_mat(1:n_forbs,1:mat_dim), n_forbs )
      call dgemm('N','T', n_forbs, n_forbs, mat_dim, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,1:mat_dim), n_forbs, &
      & tmp_mat(1:n_forbs,1:mat_dim), n_forbs, 0.0_dp, tgm_c_e%mat_Zu(1:n_forbs,1:n_forbs), n_forbs )
      if (tgm_c_e%d_Zd.gt.n_forbs) then
         mat_dim = n_el_d-1
         tgm_c_e%mat_Zd(1:n_forbs,tgm_c_e%d_Zd) = sld_c_e%mat_L_d(1:n_forbs,n_el_d)
      else
         mat_dim = n_el_d
      endif
      call dgemm('N','N', n_forbs, mat_dim, mat_dim, 1.0_dp, sld_c_e%mat_L_d(1:n_forbs,1:mat_dim), n_forbs, &
      & ovr_mat(1:mat_dim,1:mat_dim), mat_dim, 0.0_dp, tmp_mat(1:n_forbs,1:mat_dim), n_forbs )
      call dgemm('N','T', n_forbs, n_forbs, mat_dim, 1.0_dp, sld_c_e%mat_L_d(1:n_forbs,1:mat_dim), n_forbs, &
      & tmp_mat(1:n_forbs,1:mat_dim), n_forbs, 0.0_dp, tgm_c_e%mat_Zd(1:n_forbs,1:n_forbs), n_forbs )
      deallocate(tmp_mat)
   end subroutine sd2tg
   subroutine sd2pf( )
      use openmp_mpi_m, only: mpi_rank
      use molecular_system_v, only: n_fe, n_el, n_po, n_el_u, n_el_d, n_po_u, n_po_d
      use slater_determinant_params_m,               only: sld_c_e
      use pfaffian_params_m,               only: pff_c_e
      use fermionic_orbitals_m,     only: n_forbs
      use fermionic_wavefunction_v, only: wf_type_e
      integer(int32) :: i1, i2
      integer(int32) :: mat_dim
      ovr_mat = 0.0_dp 
      do i1 = 1_int32, n_el_u-1
         ovr_mat(i1,i1+1_int32) = 0.001_dp
         ovr_mat(i1+1_int32,i1) = -0.001_dp
      enddo 
      allocate( tmp_mat(1:n_forbs,1:n_el_u) ) ; tmp_mat = 0.0_dp
      call dgemm('N','T', n_forbs, n_forbs, n_el_d, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,1:n_el_d), n_forbs, &
      & sld_c_e%mat_L_d(1:n_forbs,1:n_el_d), n_forbs, 0.0_dp, pff_c_e%mat_Lmbd(1:n_forbs,1:n_forbs), n_forbs )
      !pff_c_e%mat_Lmbd(1:n_forbs,1:n_forbs) = 0.0_dp
      if (mod(n_fe,2_int32).gt.0_int32) then
         pff_c_e%vec_u(1:n_forbs) = sld_c_e%mat_L_u(1:n_forbs,n_el_u)
         mat_dim = n_el_u-1
      else
         mat_dim = n_el_u
      endif
      call dgemm('N','N', n_forbs, mat_dim, mat_dim, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,1:mat_dim), n_forbs, &
      & ovr_mat(1:mat_dim,1:mat_dim), mat_dim, 0.0_dp, tmp_mat(1:n_forbs,1:mat_dim), n_forbs )
      call dgemm('N','T', n_forbs, n_forbs, mat_dim, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,1:mat_dim), n_forbs, &
      & tmp_mat(1:n_forbs,1:mat_dim), n_forbs, 0.0_dp, pff_c_e%mat_Z_u(1:n_forbs,1:n_forbs), n_forbs )
      if (mod(n_fe,2_int32).gt.0_int32) then
         mat_dim = n_el_d-1
         pff_c_e%vec_d(1:n_forbs)  = sld_c_e%mat_L_d(1:n_forbs,n_el_d)
      else
         mat_dim = n_el_d
      endif
      call dgemm('N','N', n_forbs, mat_dim, mat_dim, 1.0_dp, sld_c_e%mat_L_d(1:n_forbs,1:mat_dim), n_forbs, &
      & ovr_mat(1:mat_dim,1:mat_dim), mat_dim, 0.0_dp, tmp_mat(1:n_forbs,1:mat_dim), n_forbs )
      call dgemm('N','T', n_forbs, n_forbs, mat_dim, 1.0_dp, sld_c_e%mat_L_d(1:n_forbs,1:mat_dim), n_forbs, &
      & tmp_mat(1:n_forbs,1:mat_dim), n_forbs, 0.0_dp, pff_c_e%mat_Z_d(1:n_forbs,1:n_forbs), n_forbs )
      !pff_c_e%mat_Z_d = 0.0_dp 
      !pff_c_e%mat_Z_u = 0.0_dp 
      deallocate(tmp_mat)
   end subroutine sd2pf
   subroutine sd2sg( )
      use openmp_mpi_m, only: mpi_rank
      use molecular_system_v, only: n_el, n_po, n_el_u, n_el_d, n_po_u, n_po_d
      use slater_determinant_params_m, only: sld_c_e, sld_c_p
      use singlet_geminal_params_m, only: sgm_c_e, sgm_c_p
      use fermionic_orbitals_m, only: n_forbs
      use fermionic_wavefunction_v, only: wf_type_e
      integer(int32) :: i1, i2
      if (n_el.gt.0_int32) then
         if ( wf_type_e.lt.0_int32 .and. wf_type_e_new.gt.0_int32) then
            if (mpi_rank.eq.0_int32) write(*,'(3X,"Converting unrestricted determinant to restricted geminal ...")')
            sgm_c_e%mat_Lmbd = 0.0_dp
            do i1 = 1_int32, n_el_d
               if (sum(ovr_mat(1:n_el_d,i1)).gt.0.0_dp ) then
                  do i2 = 1_int32, n_el_d
                     if (ovr_mat(i2,i1).gt.0.0_dp) &
                     & call dger( n_forbs, n_forbs, ovr_mat(i2,i1), sld_c_e%mat_L_u(1:n_forbs,i2), 1_int32, &
                     &  sld_c_e%mat_L_d(1:n_forbs,i1), 1_int32, sgm_c_e%mat_Lmbd(1:n_forbs,1:n_forbs), n_forbs )
                  enddo
               else
                  !if (wf_type_e_new.lt.0_int32) then
                  !  call dger( n_forbs, n_forbs, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,i1), 1_int32, &
                  !          & sld_c_e%mat_L_d(1:n_forbs,i1), 1_int32, sgm_c_e%mat_Lmbd(1:n_forbs,1:n_forbs), n_forbs )
                  !else
                  call dger( n_forbs, n_forbs, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,i1), 1_int32, &
                  & sld_c_e%mat_L_u(1:n_forbs,i1), 1_int32, sgm_c_e%mat_Lmbd(1:n_forbs,1:n_forbs), n_forbs )
                  call dger( n_forbs, n_forbs, 1.0_dp, sld_c_e%mat_L_d(1:n_forbs,i1), 1_int32, &
                  & sld_c_e%mat_L_d(1:n_forbs,i1), 1_int32, sgm_c_e%mat_Lmbd(1:n_forbs,1:n_forbs), n_forbs )
                  !endif
               endif
            enddo
         else
            call dgemm('N','T', n_forbs, n_forbs, n_el_d, 1.0_dp, sld_c_e%mat_L_u(1:n_forbs,1:n_el_d), n_forbs, &
            & sld_c_e%mat_L_d(1:n_forbs,1:n_el_d), n_forbs, 0.0_dp, sgm_c_e%mat_Lmbd(1:n_forbs,1:n_forbs), n_forbs )
         endif
         if (n_el_u.gt.n_el_d) then
            sgm_c_e%mat_Lmbd(1:n_forbs,n_forbs+1:n_forbs+n_el_u-n_el_d) = sld_c_e%mat_L_u(1:n_forbs,n_el_d+1:n_el_u)
         endif
      endif
      if (n_po.gt.0_int32) then
         call dgemm('N','T', n_forbs, n_forbs, n_po_d, 1.0_dp, sld_c_p%mat_L_u(1:n_forbs,1:n_po_d), n_forbs, &
         & sld_c_p%mat_L_d(1:n_forbs,1:n_po_d), n_forbs, 0.0_dp, sgm_c_p%mat_Lmbd(1:n_forbs,1:n_forbs), n_forbs )
         if (n_po_u.gt.n_po_d) then
            sgm_c_p%mat_Lmbd(1:n_forbs,n_forbs+1:n_forbs+n_po_u-n_po_d) = sld_c_p%mat_L_u(1:n_forbs,n_po_d+1:n_po_u)
         endif
      endif
   end subroutine sd2sg
end module wvfcnv_mod
