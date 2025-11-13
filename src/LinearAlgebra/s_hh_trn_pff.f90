! Adapted from :
! Numeric and symbolic evaluation of the pfaffian of general skew-symmetric matrices
! C. González-Ballestero, L.M. Robledo, G.F. Bertsch 
! Computer Physics Communications 182(10) (2011), 2213-2218
subroutine hh_trn_pff(n, work, ipiv, A, pff_A, A_tmp )
   use fortran_kinds_v, only: dp, int32
   implicit none
   integer(int32),                         intent(in)    :: n
   real(dp),       dimension(2*n),       intent(inout) :: work
   integer(int32), dimension(n),         intent(inout) :: ipiv   ! pivot indices
   real(dp),       dimension(n,n),     intent(inout) :: A_tmp
   real(dp),                               intent(out)   :: pff_A
   real(dp), dimension(n,n), optional, intent(inout) :: A
   real(dp), parameter :: epsln = 1.0d-13 ! smallest number such that 1+ehh_trn_pffpsln=1
   real(dp) :: xmod,u2,ki,al
   integer(int32) :: iden, info
   integer(int32) :: i1, i2
   work = 0.0_dp
   ipiv = 0_int32
   if (present(A)) A_tmp = A
   if ( n.eq.2_int32 ) then
      pff_A = A_tmp(1,2)
   elseif (n.eq.4_int32) then
      pff_A = A_tmp(1,2) * A_tmp(3,4) - A_tmp(1,3) * A_tmp(2,4) + A_tmp(1,4) * A_tmp(2,3)
   else
      pff_A = 1.0_dp
      iden  = 0_int32          ! counts the number of identity transformations
      do i1 = 1_int32, n - 2_int32, 2_int32
         i2 = n - i1
         ! Column n1-i1+1 of A_tmp is copied into column i1 of working vector
         ! This is vector x with dimension i2 of the Householder transformation
         call dcopy(i2,A_tmp(1,n-i1+1),1_int32,work(1:i2),1_int32)
         ! n1orm of x
         xmod = dsqrt( dot_product( work(1:i2),work(1:i2) ) )  ! |x|
         ! if n1orm of x=0 do nothing
         if (xmod.gt.epsln) then
            ki = sign(xmod,A_tmp(i2,i2+1))
            !   u = x -+ sign(A_tmp(n1-i1,n1-i1+1))*xmod e (n1-i1)
            !   e (i2) is the unit cartesian vector with 1 in position i2 and
            !   0.0_dp elsewhere
            work(i2) = work(i2) - ki
            !   n1orm of u
            u2 = dot_product(work(1:i2),work(1:i2))  ! |u|**2
            !  If the norm of u is too small the negative sign in the definition
            !  of u is taken
            if (u2.lt.epsln**2) then
               !  ki changes sign
               ki = -ki
               work(i2) = work(i2) - 2.0_dp*ki
               u2 = dot_product(work(1:i2),work(1:i2))  ! |u|**2
            endif
            !    v = A(n1-1) u --> work(1,n1)
            call dgemv('N',i2,i2,1.0_dp,A_tmp,n,work(1:i2),1_int32,0.0_dp,work(n+1:n+i2),1_int32)
            !      Update of A(n1-i1)
            al = 2.0_dp/u2
            !     A_tmp -> A_tmp + al*u v'
            call dger(i2,i2, al,work(1:i2),1_int32,work(n+1:n+i2),1_int32,A_tmp,n)
            !     A_tmp -> A_tmp - al*v u'
            call dger(i2,i2,-al,work(n+1:n+i2),1_int32,work(1:i2),1_int32,A_tmp,n)
         else ! xmod<=epsln
            ki = 0.0_dp
            iden = iden + 1_int32
         endif  ! xmod>epsln
         if( mod(i1,2_int32).eq.1_int32 ) pff_A = pff_A * ki
      enddo ! i1
      pff_A = pff_A * A_tmp(1,2) * dble(1-2*mod(iden,2)) * dble(1-2*mod(n/2-1,2))
   endif !    n1==2
end subroutine hh_trn_pff
