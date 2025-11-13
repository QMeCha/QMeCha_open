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
module jstqe_mod
   use fortran_kinds_v, only: dp, int32, stdout
   use write_lines_m
   use jstqe_cls,  only: jstqefct_t
   use jstqepar_var
   implicit none
   type(jstqefct_t), save, public, allocatable, dimension(:), target :: jstqe_fct
   public  :: ini_jstqe_fct, upd_jstqe_fct_par, save_jstqe_fct_par
contains
   subroutine ini_jstqe_fct( )
      use openmp_mpi_m, only: mpi_rank
      use quantum_monte_carlo_v, only: n_wlk_max
      use qdo_system_v, only: n_qdo
      use qedip_mod,  only: ini_qedip
      integer(int32) :: iw
      n_par_jstqe = 0_int32
      if( mpi_rank.eq.0_int32 ) then
         write(*,'(3X,"               ________________________________                 ")')
         write(*,'(3X,"               Initialize QDO-el Jastrow factor                 ")')
         call write_empty_line(stdout,0,mpi_rank)
      endif
      select case ( jstqe_type )
       case(0) 
         if( mpi_rank.eq.0_int32 ) write(*,'(3X,"WARNING! No QDO-el Jastrow factor is included")')
       case(1) 
         call ini_qedip()
         if( jqe_opt ) n_par_jstqe = n_par_jstqe + 3*n_qdo*3
      end select
      allocate ( jstqe_fct(1:n_wlk_max) )
      do iw = 1_int32, n_wlk_max
         call jstqe_fct(iw)%ini( )
      enddo
   end subroutine ini_jstqe_fct
   subroutine upd_jstqe_fct_par( vec_jstqe_fct_var )
      use qedip_mod, only: upd_qedip_par
      real(dp), dimension(n_par_jstqe), intent(in) :: vec_jstqe_fct_var
      select case (jstqe_type)
       case(1) 
         call upd_qedip_par( vec_jstqe_fct_var )
      end select
   end subroutine upd_jstqe_fct_par
   subroutine save_jstqe_fct_par( )
      use qedip_mod, only: save_qedip
      select case (jstqe_type)
       case(1)
         call save_qedip()
      end select
   end subroutine save_jstqe_fct_par
end module jstqe_mod
