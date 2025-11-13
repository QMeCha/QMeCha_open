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
subroutine space_map( d_fn, sm_type, i_at_min, f_spc_map )
   use fortran_kinds_v,    only: dp, int32
   use molecular_system_v, only: n_at, atoms, dist_t
   implicit none
   type(dist_t),   dimension(n_at), intent(in)  :: d_fn
   integer(int32),                  intent(in)  :: sm_type
   integer(int32),                  intent(out) :: i_at_min
   real(dp),                        intent(out) :: f_spc_map
   real(dp)                                     :: d_fn_min, L
   integer(int32) :: ia
   i_at_min  = 0_int32
   d_fn_min  = 100.0_dp
   f_spc_map = 1.0_dp
   if (sm_type.eq.0_int32) return
   if (n_at.gt.0_int32) then
      select case(sm_type)
         case(1)
            do ia = 1_int32, n_at
               if ( atoms(ia)%atm_z.gt.1_int32 ) then
                  L =  (1.0_dp + sqrt(dble(atoms(i_at_min)%atm_z)))
                  if ( d_fn(ia)%m .le. L .and. ( d_fn_min.gt.d_fn(ia)%m ) ) then
                     d_fn_min = d_fn(ia)%m
                     i_at_min = ia
                  endif
               endif
            enddo
         case default
            do ia = 1_int32, n_at
               if ( atoms(ia)%atm_z.gt.1_int32 .and. d_fn_min.gt.d_fn(ia)%m ) then
                  d_fn_min = d_fn(ia)%m
                  i_at_min = ia
               endif
            enddo
      end select
   endif
   if (i_at_min.gt.0_int32) then
      select case(sm_type)
       case(1)
         L =  (1.0_dp + sqrt(dble(atoms(i_at_min)%atm_z)))
         f_spc_map = 1.0_dp - (1.0_dp - 1.0_dp / (dble(atoms(i_at_min)%atm_z))**0.25) * (1.0_dp - d_fn(i_at_min)%m / L)**4
       case(2)
         f_spc_map = dble(atoms(i_at_min)%atm_z)
         f_spc_map = (1.0_dp + f_spc_map * d_fn(i_at_min)%m ) / (1.0_dp + d_fn(i_at_min)%m ) / f_spc_map
      end select
   endif
end subroutine space_map
