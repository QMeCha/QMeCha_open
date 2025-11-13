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
module prdqcfs_mod
use fortran_kinds_v, only: dp, int32 
use openmp_mpi_m
use qdo_wavefunction_v, only: ql_opt
use drudonic_orbitals_m, only: n_dorbs
implicit none
type, public :: prdqcfs_t
  real(dp),      allocatable, dimension(:,:)  :: mat_L
end type prdqcfs_t 
type(prdqcfs_t),  public, save, target       :: prdq_c_d
public   :: ini_prdqcfs, upd_prdqcfs, save_prdqcfs
private  :: read_prdqcfs
contains
subroutine ini_prdqcfs( prdq_c, n_qdo, n_par_drd, n_par_drd_s, head_fle )
  use quantum_monte_carlo_v,  only: restart
  type(prdqcfs_t),         intent(inout) :: prdq_c
  integer(int32),         intent(in)    :: n_qdo
  integer(int32),         intent(inout) :: n_par_drd
  integer(int32),         intent(inout) :: n_par_drd_s
  character(4), optional, intent(in)    :: head_fle
  integer(int32)  :: i1 
  allocate( prdq_c%mat_L(1:n_dorbs,1:n_qdo) ) ; prdq_c%mat_L = 0.0_dp
  do i1 = 1_int32, n_qdo
    prdq_c%mat_L(i1,i1) = 1.0_dp
  enddo ! i1 (n_qdo)
  if (restart) call read_prdqcfs( prdq_c, n_qdo, head_fle )
  if (ql_opt) then
    n_par_drd = n_dorbs*n_qdo
    n_par_drd_s = n_dorbs*n_qdo
  else
    n_par_drd = 0_int32
    n_par_drd_s = 0_int32
  endif 
end subroutine ini_prdqcfs
subroutine upd_prdqcfs( prdq_c, n_qdo, n_par_drd, vec_prdq_var )
    type(prdqcfs_t),  target, intent(inout) :: prdq_c 
    integer(int32),       intent(in)    :: n_qdo
    integer(int32),       intent(in)    :: n_par_drd
    real(dp), dimension(n_par_drd), intent(in) :: vec_prdq_var
    integer(int32) :: i1, i2, ip
    ip = 0_int32
    do i1 = 1_int32, n_qdo ; do i2 = 1, n_dorbs 
        prdq_c%mat_L(i2,i1) = prdq_c%mat_L(i2,i1) + vec_prdq_var(ip+i2)
        !prdq_c%mat_L(i2,i1) = max(prdq_c%mat_L(i2,i1) + vec_prdq_var(ip+i2), 0.0_dp )
      enddo
      ip = ip + n_dorbs
    enddo 
end subroutine upd_prdqcfs
subroutine save_prdqcfs( prdq_c, n_qdo, head_fle )
  type(prdqcfs_t),  target, intent(inout) :: prdq_c 
  integer(int32),         intent(in)    :: n_qdo
  character(4), optional, intent(in)    :: head_fle
  integer(int32) :: i1, i2
  character(100) :: svf_fle
  svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
  open(unit=1,file=svf_fle,action='write',form='formatted',status='unknown')
  write(1,'("# QDO product wave function coefficients")')
  i1 = 1_int32
  if (n_qdo.le.6_int32) then
    do i2 = 1_int32, n_dorbs
      write(1,'(6E18.9)') prdq_c%mat_L(i2,1:n_qdo)
    enddo  
  else
    do i1 = 1_int32, int(n_qdo/6_int32)
      do i2 = 1_int32, n_dorbs
        write(1,'(6E18.9)') prdq_c%mat_L(i2,6*(i1-1)+1:6*i1)  
      enddo
    enddo
    if( mod(n_qdo, 6_int32).gt.0_int32 ) then
      do i2 = 1_int32, n_dorbs
        write(1,'(6E18.9)') prdq_c%mat_L(i2,n_qdo-mod(n_qdo, 6_int32)+1:n_qdo)  
      enddo
    endif 
  endif
  close(1)
end subroutine save_prdqcfs
subroutine read_prdqcfs( prdq_c, n_qdo, head_fle )
  type(prdqcfs_t),  target, intent(inout) :: prdq_c 
  integer(int32),         intent(in)    :: n_qdo
  character(4), optional, intent(in)    :: head_fle
  integer(int32) :: i1, i2
  character(100) :: svf_fle
  svf_fle = 'wvfn.save/' // trim(head_fle) // '.sav'
  if( mpi_rank.eq.0_int32 ) then
    open(unit=1,file=svf_fle,action='read',form='formatted',status='old')
    read(1,*) ! '("# QDO product wave function coefficients")'
    i1 = 1_int32
    if (n_qdo.le.6_int32) then
      do i2 = 1_int32, n_dorbs
        read(1,*) prdq_c%mat_L(i2,1:n_qdo)
      enddo  
    else
      do i1 = 1_int32, int(n_qdo/6_int32)
        do i2 = 1_int32, n_dorbs
          read(1,*) prdq_c%mat_L(i2,6*(i1-1)+1:6*i1)  
        enddo
      enddo
      if( mod(n_qdo, 6_int32).gt.0_int32 ) then
        do i2 = 1_int32, n_dorbs
          read(1,*) prdq_c%mat_L(i2,n_qdo-mod(n_qdo, 6_int32)+1:n_qdo)  
        enddo
      endif 
    endif
    close(1)
  endif 
#if defined _MPI || defined _MPIh || defined _MPI08
  call mpi_bcast(prdq_c%mat_L, n_dorbs*n_qdo, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine read_prdqcfs
end module prdqcfs_mod
