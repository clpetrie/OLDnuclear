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
   subroutine rsqsvd(a,u,sigma,vt,n) 
! 
!  svd of a real square matrix a. a=u*sigma*vt.
!
   integer(kind=i4), intent(in) :: n
   real(kind=r8), intent(in) :: a(n,n)
   real(kind=r8), intent(out) :: sigma(n),vt(n,n),u(n,n)
   real(kind=r8) :: work(10*n),atmp(n,n)
   integer(kind=i4) :: info,lwork
   lwork=size(work) !must be >= 5n
   if (n.lt.1) return
   atmp=a !dgesvd destroys a
   call dgesvd('All','All',n,n,atmp,n,sigma,u,n,vt,n,work,lwork,info)
   if (info.ne.0) then
      write (6,'(''error in dgesvd '',i10)') info
      call abort
   endif
   end subroutine rsqsvd
    
   subroutine eigenrs(eigvec,eigval,n)
!
! compute eigenvalues and eigenvectors of a real symmetric matrix
! input matrix in eigvec (only lower triangle used), output eigenvectors
! in eigvec, output eigenvalues in eigval
!
   integer(kind=i4) :: n,info
!
! lwork >= (nb+2)*n so assume block size 128 is more than enough
!
   real(kind=r8) :: eigvec(n,n),eigval(n),work(130*n)
   integer(kind=i4) :: i
   if (n.lt.1) return
   call dsyev('v','l',n,eigvec,n,eigval,work,130*n,info)
   if (info.ne.0) then
      write (6,'(1x,''error in dsyev'',i10)') info
      call abort()
   endif
   do i=1,n
      if (eigvec(1,i).lt.0.0_r8) eigvec(:,i)=-eigvec(:,i)
   enddo
   end subroutine eigenrs

   subroutine cmatinv(a,det,n)
   integer(kind=i4) :: n,ipiv(n),info,i
   complex (kind=r8) :: a(n,n),det,cwork(n,n)
!
! lapack routine for lu factorization
!
   call zgetrf(n,n,a,n,ipiv,info)
   if (info.lt.0) then
      write (6,'(1x,''error in zgetrf'',i10)') info
      call abort()
   endif
   if (info.gt.0) then !determinant is zero -- a is undefined -- set to zero
      det=czero
      a=czero
      return
   endif
!
! calculate determinant
!
   det=cone
   do i=1,n
      det=det*a(i,i)
      if (ipiv(i).ne.i) det=-det
   enddo
!
! lapack routine to calculate inverse from factorization
!
   call zgetri(n,a,n,ipiv,cwork,n*n,info)
   if (info.ne.0) then
      write (6,'(1x,''error in zgetri'',i10)') info
      call abort()
   endif
   end subroutine cmatinv

   subroutine rmatinv(a,detl,is,n)
! calculate inverse, log(abs(det)), and sign.
   integer(kind=i4) :: n,is,ipiv(n),i,info
   real(kind=r8) :: a(n,n),detl,work(n*n)
   real(kind=r8), parameter :: minushuge=-1e30_r8
   call dgetrf(n,n,a,n,ipiv,info)
   if (info.ne.0) then
      write (6,'(1x,''error in dgetrf'',i10)') info
      stop
   endif
   is=1
   detl=0.0_r8
   do i=1,n
      detl=detl+log(abs(a(i,i)))
      is=is*sign(1.0_r8,a(i,i))
      if (ipiv(i).ne.i) is=-is
   enddo
   call dgetri(n,a,n,ipiv,work,n*n,info)
   if (info.ne.0) then
      write (6,'(1x,''error in dgetri'',i10)') info
      stop
   endif
   end subroutine rmatinv

   function expmult(a,vec,n)
!
! routine to perform vec=exp(a)*vec for a complex matrix a and vector vec
!
   integer(kind=i4) :: n,info
   complex(kind=r8) :: a(n,n),vec(n),expmult(n)
   complex(kind=r8) :: det,val(n),vecr(n,n),vecl(n,n),cwork(n,n)
   real(kind=r8) :: work(2*n)
   call zgeev('n','v',n,a,n,val,vecl,n,vecr,n,cwork,n*n,work,info)
   if (info.ne.0) then
      write (6,'(1x,''problem in expmult'',i5)') info
      stop
   endif
   vecl=vecr
   call matinv(vecl,det,n)
   val=exp(val)
   expmult=matmul(vecr,val*matmul(vecl,vec))
   end function expmult

end module matrixmod
