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
module cusps_functions_params_m
   use fortran_kinds_v, only: dp, int32
   use openmp_mpi_m
   implicit none
   type, public :: jcspar_t
      integer(int32)                        :: n_cmp
      integer(int32)                        :: typ, ord
      integer(int32)                        :: n_par, n_opt_par
      real(dp), allocatable, dimension(:,:) :: b
   end type jcspar_t
   logical,        save, public :: jcen_prs, jcpn_prs
   logical,        save, public :: jcee_prs, jcpp_prs, jcep_prs
   type(jcspar_t), save, public, target :: b_en, b_pn
   type(jcspar_t), save, public, target :: b_ee, b_pp, b_ep
   integer(int32), save, public :: n_par_jcen, n_par_jcpn
   integer(int32), save, public :: n_par_jcee, n_par_jcpp, n_par_jcep
   public :: save_cusps_functions_params, read_cusps_functions_params
contains
   subroutine save_cusps_functions_params( f_indx, b_ff )
      integer(int32), intent(in) :: f_indx
      type(jcspar_t), intent(in) :: b_ff
      integer(int32)             :: i1
      write(f_indx,'(5I3)') b_ff%typ, b_ff%ord, b_ff%n_cmp, b_ff%n_par
      if (b_ff%n_par.gt.0) then
         do i1 = 1_int32, b_ff%n_cmp
            if (b_ff%typ.ne.4_int32) then
               write(f_indx,'(100E18.9)') b_ff%b(1,i1)
               if (b_ff%ord.gt.0)then
                  write(f_indx,'(100E18.9)') b_ff%b(2:b_ff%ord+1,i1)
                  if (b_ff%typ.ne.3) write(f_indx,'(100E18.9)') b_ff%b(b_ff%ord+2:2*b_ff%ord+1,i1)
               endif
            else
               write(f_indx,'(100E18.9)') b_ff%b(1,i1), b_ff%b(2,i1)
               if (b_ff%ord.gt.0)then
                  write(f_indx,'(100E18.9)') b_ff%b(3:b_ff%ord+2,i1)
                  if (b_ff%typ.ne.3) write(f_indx,'(100E18.9)') b_ff%b(b_ff%ord+3:2*b_ff%ord+2,i1)
               endif
            endif
         enddo 
      endif
   end subroutine save_cusps_functions_params
   subroutine read_cusps_functions_params( f_indx, b_ff )
      integer(int32), intent(in)    :: f_indx
      type(jcspar_t), intent(inout) :: b_ff
      integer(int32) :: typ, ord, n_cmp, n_par
      integer(int32) :: i1
      if ( mpi_rank.eq.0_int32 ) then
         read(f_indx,*) typ, ord, n_cmp, n_par
         if (n_par.gt.0_int32) then
            if ( b_ff%typ.eq.typ .and. b_ff%n_par.ge.n_par) then
               do i1 = 1_int32, n_cmp
                  if ( typ.ne.4 ) then
                     read(f_indx,*) b_ff%b(1,i1)
                     if (ord.gt.0)then
                        read(f_indx,*) b_ff%b(2:ord+1,i1)
                        if (typ.ne.3) read(f_indx,*) b_ff%b(ord+2:2*ord+1,i1)
                     endif
                  else
                     read(f_indx,*) b_ff%b(1:2,i1)
                     if (ord.gt.0)then
                        read(f_indx,*) b_ff%b(3:ord+2,i1)
                        if (typ.ne.3) read(f_indx,*) b_ff%b(ord+3:2*ord+2,i1)
                     endif
                  endif
               enddo 
            else
               do i1 = 1_int32, n_cmp
                  read(f_indx,*) 
                  if (ord.gt.0)then
                     read(f_indx,*)
                     if (typ.ne.3) read(f_indx,*) 
                  endif
               enddo 
            endif
         endif
      endif
#if defined _MPI || defined _MPIh || defined _MPI08
      call mpi_bcast( b_ff%b(1:b_ff%n_par,1:b_ff%n_cmp), b_ff%n_par*b_ff%n_cmp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr )
#endif
   end subroutine read_cusps_functions_params
end module cusps_functions_params_m
