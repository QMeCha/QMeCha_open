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
module merge_sort_m
   use fortran_kinds_v, only: int32, dp
   implicit none
   interface mrg_srt
      procedure mrg_srt_int
      procedure mrg_srt_dbl
   end interface mrg_srt
   interface mrg_srt_abs
      procedure mrg_srt_abs_int
   end interface mrg_srt_abs
   interface mrg_srt_col
      procedure mrg_srt_col_dbl
   end interface mrg_srt_col
contains
   recursive subroutine mrg_srt_int( n_ele, lst )
      integer(int32), intent(in)                        :: n_ele
      integer(int32), dimension(n_ele), intent(inout) :: lst
      integer(int32)                                    :: n_mid, tmp_ele
      if (n_ele .eq. 2_int32) then
         if ( lst(1).gt.lst(2) ) then
            tmp_ele = lst(1)
            lst(1) = lst(2)
            lst(2) = tmp_ele
         endif
      elseif (n_ele .gt. 2_int32) then
         n_mid = n_ele / 2_int32
         call mrg_srt_int(n_mid, lst(1:n_mid))
         call mrg_srt_int(n_ele-n_mid, lst(n_mid+1:n_ele))
         if ( lst(n_mid).gt.lst(n_mid+1) ) call mrg_int(n_ele, n_mid, lst )
      endif
   end subroutine mrg_srt_int
   recursive subroutine mrg_srt_dbl( n_ele, lst )
      integer(int32),               intent(in)    :: n_ele
      real(dp), dimension(n_ele), intent(inout) :: lst
      integer(int32)                              :: n_mid
      real(dp)                                    :: tmp_ele
      if (n_ele .eq. 2_int32) then
         if ( lst(1).gt.lst(2) ) then
            tmp_ele = lst(1)
            lst(1) = lst(2)
            lst(2) = tmp_ele
         endif
      elseif (n_ele .gt. 2_int32) then
         n_mid = n_ele / 2_int32
         call mrg_srt_dbl(n_mid, lst(1:n_mid))
         call mrg_srt_dbl(n_ele-n_mid, lst(n_mid+1:n_ele))
         if ( lst(n_mid).gt.lst(n_mid+1) ) call mrg_dbl(n_ele, n_mid, lst )
      endif
   end subroutine mrg_srt_dbl
   recursive subroutine mrg_srt_abs_int( n_ele, lst )
      integer(int32), intent(in)                        :: n_ele
      integer(int32), dimension(n_ele), intent(inout) :: lst
      integer(int32)                                    :: n_mid, tmp_ele
      if (n_ele .eq. 2_int32) then
         if ( abs(lst(1)).gt.abs(lst(2)) ) then
            tmp_ele = lst(1)
            lst(1) = lst(2)
            lst(2) = tmp_ele
         endif
      elseif (n_ele .gt. 2_int32) then
         n_mid = n_ele / 2_int32
         call mrg_srt_abs_int(n_mid, lst(1:n_mid))
         call mrg_srt_abs_int(n_ele-n_mid, lst(n_mid+1:n_ele))
         if ( abs(lst(n_mid)).gt.abs(lst(n_mid+1)) ) call mrg_abs_int(n_ele, n_mid, lst )
      endif
   end subroutine mrg_srt_abs_int
   recursive subroutine mrg_srt_col_dbl( n_ele, n_col, i_col, lst )
      integer(int32), intent(in)                        :: n_ele
      integer(int32), intent(in)                        :: n_col
      integer(int32), intent(in)                        :: i_col
      real(dp), dimension(n_ele,n_col), intent(inout) :: lst
      integer(int32)                                    :: n_mid
      real(dp), dimension(n_col)                      :: tmp_ele
      if (n_ele .eq. 2_int32) then
         if ( lst(1,i_col).gt.lst(2,i_col) ) then
            tmp_ele = lst(1,1:n_col)
            lst(1,1:n_col) = lst(2,1:n_col)
            lst(2,1:n_col) = tmp_ele
         endif
      elseif (n_ele .gt. 2_int32) then
         n_mid = n_ele / 2_int32
         call mrg_srt_col_dbl( n_mid, n_col, i_col, lst(1:n_mid,1:n_col) )
         call mrg_srt_col_dbl( n_ele-n_mid, n_col, i_col, lst(n_mid+1:n_ele,1:n_col) )
         if ( lst(n_mid,i_col).gt.lst(n_mid+1,i_col) ) &
         & call mrg_col_dbl( n_ele, n_col, i_col, n_mid, lst )
      endif
   end subroutine mrg_srt_col_dbl
   subroutine mrg_int( n_ele, n_mid , lst )
      integer(int32), intent(in)                        :: n_ele
      integer(int32), intent(in)                        :: n_mid
      integer(int32), dimension(n_ele), intent(inout) :: lst
      integer(int32), dimension(n_mid)                :: lst_tmp
      integer(int32)                                    :: i1, i2, i3
      i1 = 1_int32
      i2 = n_mid + 1_int32
      i3 = 1_int32
      lst_tmp(1:n_mid) = lst(1:n_mid)
      do while (i1 .le. n_mid .and. i2.le.n_ele )
         if( lst_tmp(i1) .gt. lst(i2) ) then
            lst(i3) = lst(i2)
            i2 = i2 + 1_int32
         else
            lst(i3) = lst_tmp(i1)
            i1 = i1 + 1_int32
         endif
         i3 = i3 + 1_int32
      enddo 
      if (i2.gt.n_ele  ) lst(i3:n_ele) =  lst_tmp(i1:n_mid)
   end subroutine mrg_int
   subroutine mrg_dbl( n_ele, n_mid , lst )
      integer(int32), intent(in)                        :: n_ele
      integer(int32), intent(in)                        :: n_mid
      real(dp),       dimension(n_ele), intent(inout) :: lst
      real(dp),       dimension(n_mid)                :: lst_tmp
      integer(int32)                                    :: i1, i2, i3
      i1 = 1_int32
      i2 = n_mid + 1_int32
      i3 = 1_int32
      lst_tmp(1:n_mid) = lst(1:n_mid)
      do while (i1 .le. n_mid .and. i2.le.n_ele )
         if( lst_tmp(i1) .gt. lst(i2) ) then
            lst(i3) = lst(i2)
            i2 = i2 + 1_int32
         else
            lst(i3) = lst_tmp(i1)
            i1 = i1 + 1_int32
         endif
         i3 = i3 + 1_int32
      enddo 
      if ( i2.gt.n_ele ) lst(i3:n_ele) =  lst_tmp(i1:n_mid)
   end subroutine mrg_dbl
   subroutine mrg_abs_int( n_ele, n_mid , lst )
      integer(int32), intent(in)                        :: n_ele
      integer(int32), intent(in)                        :: n_mid
      integer(int32), dimension(n_ele), intent(inout) :: lst
      integer(int32), dimension(n_mid)                :: lst_tmp
      integer(int32)                                    :: i1, i2, i3
      i1 = 1_int32
      i2 = n_mid + 1_int32
      i3 = 1_int32
      lst_tmp(1:n_mid) = lst(1:n_mid)
      do while (i1 .le. n_mid .and. i2.le.n_ele )
         if( abs(lst_tmp(i1)) .gt. abs(lst(i2)) ) then
            lst(i3) = lst(i2)
            i2 = i2 + 1_int32
         else
            lst(i3) = lst_tmp(i1)
            i1 = i1 + 1_int32
         endif
         i3 = i3 + 1_int32
      enddo 
      if ( i2.gt.n_ele ) lst(i3:n_ele) =  lst_tmp(i1:n_mid)
   end subroutine mrg_abs_int
   subroutine mrg_col_dbl( n_ele, n_col, i_col, n_mid, lst )
      integer(int32), intent(in)                        :: n_ele
      integer(int32), intent(in)                        :: n_col
      integer(int32), intent(in)                        :: i_col
      integer(int32), intent(in)                        :: n_mid
      real(dp), dimension(n_ele, 1:n_col), intent(inout) :: lst
      real(dp), dimension(n_mid, 1:n_col)               :: lst_tmp
      integer(int32)                                    :: i1, i2, i3
      i1 = 1_int32
      i2 = n_mid + 1_int32
      i3 = 1_int32
      lst_tmp(1:n_mid, 1:n_col) = lst(1:n_mid, 1:n_col)
      do while (i1 .le. n_mid .and. i2.le.n_ele )
         if( lst_tmp(i1, i_col) .gt. lst(i2, i_col) ) then
            lst(i3, 1:n_col) = lst(i2, 1:n_col)
            i2 = i2 + 1_int32
         else
            lst(i3, 1:n_col) = lst_tmp(i1, 1:n_col)
            i1 = i1 + 1_int32
         endif
         i3 = i3 + 1_int32
      enddo 
      if ( i2.gt.n_ele ) lst(i3:n_ele, 1:n_col) =  lst_tmp(i1:n_mid, 1:n_col)
   end subroutine mrg_col_dbl
end module merge_sort_m
