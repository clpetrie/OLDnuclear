module matrixmod
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: cone = (1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero = (0.0_r8,0.0_r8)

interface matinv
   module procedure rmatinv,cmatinv
end interface

contains
   subroutine rmatinv(a,detl,is,n)
   integer(kind=i4) :: n,info,ipiv(n),i,is
   real(kind=r8) :: a(n,n),detl,b(n,n)
   call dgetrf(n,n,a,n,ipiv,info)
!
! calculate determinant
!
   detl=0.0_r8
   do i=1,n
      detl=detl+log(abs(a(i,i)))
      if (a(i,i).lt.0.0_r8) is=-is
      if (ipiv(i).ne.i) is=-is
   enddo
!
! form inverse
!
   b=0.0_r8
   do i=1,n
      b(i,i)=1.0_r8
   enddo
   call dgetrs('n',n,n,a,n,ipiv,b,n,info)
   a=b
   end subroutine rmatinv

   subroutine eigenrs(eigvec,eigval,n)
!
! compute eigenvalues and eigenvectors of a real symmetric matrix
! input matrix in eigvec (only lower triangle used), output eigenvectors
! in eigvec, output eigenvalues in eigval
!
   integer(kind=i4) :: n,i,j,k
   real(kind=r8) :: eigvec(n,n),eigval(n),aux(2*n),ap((n*(n+1))/2)
   if (n.lt.1) return
!
! pack matrix this uses the lower triangle only
!
   k=0
   do i=1,n
      do j=i,n
         k=k+1
         ap(k)=eigvec(j,i)
      enddo
   enddo
   call dspev(1,ap,eigval,eigvec,n,n,aux,2*n)
   end subroutine eigenrs

   subroutine cmatinv(a,det,n)
   integer(kind=i4) :: n,info,ipiv(n),i
   complex (kind=r8) :: a(n,n),det,b(n,n)
   call zgetrf(n,n,a,n,ipiv,info)
!
! calculate determinant
!
   det=cone
   do i=1,n
      det=det*a(i,i)
      if (ipiv(i).ne.i) det=-det
   enddo
   b=czero
   do i=1,n
      b(i,i)=cone
   enddo
   call zgetrs('n',n,n,a,n,ipiv,b,n,info)
   a=b
   end subroutine cmatinv

   function expmult(a,vec,n)
!
! routine to perform vec=exp(a)*vec for a complex matrix a and vector vec
! a is destroyed
!
   integer(kind=i4) :: n,j,info
   complex(kind=r8) :: a(n,n),vec(n),expmult(n)
   complex(kind=r8) :: val(n),vecr(n,n),vecl(n,n),cwork(3*n),det
   call zgeev(1,a,n,val,vecr,n,.true.,n,cwork,3*n)
   vecl=vecr
   call cmatinv(vecl,det,n)
   val=exp(val)
   expmult=matmul(vecr,val*matmul(vecl,vec))
   end function expmult

end module matrixmod
