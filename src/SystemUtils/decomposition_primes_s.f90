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
subroutine decomposition_primes( x, order_of_primes)
   use fortran_kinds_v, only: int32
   implicit none
   integer(int32), intent(in)  :: x
   integer(int32), intent(out) :: order_of_primes(1:4)
   integer(int32) :: x_tmp 
   order_of_primes(1:4) = 1_int32
   x_tmp = x
   do 
      if (mod(x_tmp,2**order_of_primes(1)).ne.0_int32) then
         order_of_primes(1) = order_of_primes(1) - 1_int32 
         exit
      endif
      order_of_primes(1) = order_of_primes(1) + 1_int32  
   enddo
   x_tmp = x_tmp / 2**order_of_primes(1)
   do 
      if (mod(x_tmp,3**order_of_primes(2)).ne.0_int32) then
         order_of_primes(2) = order_of_primes(2) - 1_int32 
         exit
      endif
      order_of_primes(2) = order_of_primes(2) + 1_int32  
   enddo
   x_tmp = x_tmp / 3**order_of_primes(2)
   do 
      if (mod(x_tmp,5**order_of_primes(3)).ne.0_int32) then
         order_of_primes(3) = order_of_primes(3) - 1_int32 
         exit
      endif
      order_of_primes(3) = order_of_primes(3) + 1_int32  
   enddo
   x_tmp = x_tmp / 5**order_of_primes(3)
   do 
      if (mod(x_tmp,7**order_of_primes(4)).ne.0_int32) then
         order_of_primes(4) = order_of_primes(4) - 1_int32 
         exit
      endif
      order_of_primes(4) = order_of_primes(4) + 1_int32  
   enddo
   x_tmp = x_tmp / 7**order_of_primes(4)
end subroutine decomposition_primes
