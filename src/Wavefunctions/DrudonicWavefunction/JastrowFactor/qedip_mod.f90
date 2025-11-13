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
module qedip_mod
   use fortran_kinds_v, only: dp, int32
   use qdo_system_v, only: n_qdo
   use jstqepar_var
   implicit none
   real(dp),     allocatable, dimension(:,:) :: mat_Aqe
   real(dp),                  dimension(3)   :: mu_nuc
   public  :: ini_qedip, save_qedip, upd_qedip_par
   private :: read_qedip
contains
   subroutine ini_qedip()
      use openmp_mpi_m
      use quantum_monte_carlo_v,  only: restart
      use molecular_system_v,  only: n_at, atoms, r_at
      integer(int32) :: i1
      allocate( mat_Aqe(1:3,1:3*n_qdo) ) ; mat_Aqe = 0.0_dp
      if(restart) then
         if( mpi_rank.eq.0_int32 ) write(*,'(3X,"Reading QDO-electron dip coupling matrix... ")')
         call read_qedip()
      else
         if( mpi_rank.eq.0_int32 ) write(*,'(3X,"Initialize QDO-electron dip coupling matrix... ")')
      endif
      mu_nuc(:) = 0
      do i1 = 1_int32, n_at
         mu_nuc(:) = - atoms(i1)%atm_z * r_at(:,i1)
      enddo 
   end subroutine ini_qedip
   subroutine save_qedip()
      integer(int32) :: i1, i2
      character(100) :: svf_fle
      svf_fle = 'wvfn.save/jstqedip.sav'
      open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
      write(1,'("# QDO-el dipole matrix coefficients")')
      do i1 = 1_int32, n_qdo
         write(1,"(a,I3,a,I3)") 'electrons-QDO_i= ', i1
         do i2 = 1_int32, 3_int32
            write(1,*) mat_Aqe( 1:3, 3*(i1-1)+i2 )
         enddo 
      enddo 
      close(1)
   end subroutine save_qedip
   subroutine read_qedip()
      use openmp_mpi_m
      integer(int32) :: i1, i2
      character(100) :: svf_fle
      if(mpi_rank.eq.0_int32) then
         svf_fle = 'wvfn.save/jstqedip.sav'
         open(unit=1,file=svf_fle,action='read',form='formatted',status='unknown')
         read(1,*) ! Comment line
         do i1 = 1_int32, n_qdo
            read(1,*) ! Comment line
            do i2 = 1_int32, 3_int32
               read(1,*) mat_Aqe( 1:3, 3*(i1-1)+i2 )
            enddo 
         enddo 
         close(1)
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast( mat_Aqe(1:3,1:3*n_qdo), 9*n_qdo, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr )
      call mpi_barrier(MPI_COMM_WORLD, mpierr)
#endif
   end subroutine read_qedip
   subroutine upd_qedip_par( vec_qedip_fct_var )
      real(dp), dimension(n_par_jstqe), intent(in) :: vec_qedip_fct_var
      integer(int32) :: i1, i2, i3
      i3 = 0_int32
      do i1 = 1_int32, 3_int32
         do i2 = 1_int32, 3_int32*n_qdo
            i3 = i3 + 1_int32
            mat_Aqe(i1,i2) = mat_Aqe(i1,i2) + vec_qedip_fct_var(i3)
         enddo 
      enddo 
   end subroutine upd_qedip_par
end module qedip_mod
