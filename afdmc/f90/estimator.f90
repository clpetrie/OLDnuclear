module estimator
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   character(len=30), private, dimension(:), allocatable, save :: label
   complex(kind=r8), private, dimension(:), allocatable, save :: valblk,valnow,wtblk
   complex(kind=r8), private, dimension(:), allocatable, save :: valblk2,errnow
   integer(kind=i4), private, save :: nest
   integer(kind=i4), private, save, allocatable :: nval(:)
   complex(kind=r8), private, save, allocatable :: cvalblk(:),cvalnow(:),cwtblk(:),cvalblk2(:)
   real(kind=r8), private, save, allocatable :: cerrnow(:)

interface addval
   module procedure addvalrr,addvalcc,addvalrc,addvalcr
end interface addval

contains
   subroutine setestnum(n)
   integer(kind=i4) :: n
   nest=n
   allocate(valblk(n),valnow(n))
   allocate(valblk2(n),errnow(n))
   allocate(wtblk(n))
   allocate(label(n))
   allocate(nval(n))
   allocate(cvalblk(n),cvalnow(n),cwtblk(n),cvalblk2(n),cerrnow(n))
   return
   end subroutine setestnum

   subroutine zerest
   integer(kind=i4) :: i
   do i=1,nest
      valblk(i)=0.0_r8
      valblk2(i)=0.0_r8
      valnow(i)=0.0_r8
      errnow(i)=0.0_r8
      wtblk(i)=0.0_r8
      nval(i)=0
      cvalblk(i)=0.0_r8
      cvalnow(i)=0.0_r8
      cwtblk(i)=0.0_r8
      cvalblk2(i)=0.0_r8
      cerrnow(i)=0.0_r8
   enddo
   return
   end subroutine zerest

   subroutine addest(i,l)
   integer(kind=i4) :: i
   character(len=*) :: l
   integer(kind=i4) :: ln,j
   ln=min(len(l),len(label(i)))
   label(i)(1:ln)=l(1:ln)
   do j=ln+1,len(label(i))
      label(i)(j:j)=" "
   enddo
   return
   end subroutine addest

   subroutine addvalrr(i,val,wt)
   integer(kind=i4) :: i
   real(kind=r8) :: val,wt
   valblk(i)=valblk(i)+val*wt
   valblk2(i)=valblk2(i)+val*val*wt
   wtblk(i)=wtblk(i)+wt
   nval(i)=nval(i)+1
   return
   end subroutine addvalrr

   subroutine addvalcc(i,val,wt)
   integer(kind=i4) :: i
   complex(kind=r8) :: val,wt
   valblk(i)=valblk(i)+val*wt
!  valblk2(i)=valblk2(i)+val*val*wt
   valblk2(i)=valblk2(i)+dconjg(val)*val*wt
   wtblk(i)=wtblk(i)+wt
   nval(i)=nval(i)+1
   return
   end subroutine addvalcc

   subroutine addvalrc(i,val,wt)
   integer(kind=i4) :: i
   real(kind=r8) :: val
   complex(kind=r8) :: wt
   valblk(i)=valblk(i)+val*wt
   valblk2(i)=valblk2(i)+val*val*wt
   wtblk(i)=wtblk(i)+wt
   nval(i)=nval(i)+1
   return
   end subroutine addvalrc

   subroutine addvalcr(i,val,wt)
   integer(kind=i4) :: i
   complex(kind=r8) :: val
   real(kind=r8) :: wt
   valblk(i)=valblk(i)+val*wt
!  valblk2(i)=valblk2(i)+val*val*wt
   valblk2(i)=valblk2(i)+dconjg(val)*val*wt
   wtblk(i)=wtblk(i)+wt
   nval(i)=nval(i)+1
   return
   end subroutine addvalcr

   subroutine update
   use mympi
   integer(kind=i4) :: i,nvalsum(nest)
   complex(kind=r8) :: valsum(nest),wtsum(nest),valsum2(nest)
   complex(kind=r8) :: cvalsum(nest),cwtsum(nest),cvalsum2(nest)
   call addall(valblk,valsum)
   call addall(valblk2,valsum2)
   call addall(wtblk,wtsum)
   call addall(nval,nvalsum)
   call addall(cvalblk,cvalsum)
   call addall(cvalblk2,cvalsum2)
   call addall(cwtblk,cwtsum)
   if (myrank().eq.0) then
      valblk=valsum
      valblk2=valsum2
      wtblk=wtsum
!     wtblk=real(wtsum)
      nval=nvalsum
      cvalblk=cvalsum
      cvalblk2=cvalsum2
      cwtblk=cwtsum
      do i=1,nest
         valnow(i)=valblk(i)/wtblk(i)
         valblk2(i)=valblk2(i)/wtblk(i)
!        errnow(i)=sqrt(abs(valblk2(i)-valnow(i)**2)/wtblk(i))
         errnow(i)=sqrt(abs(valblk2(i)-valnow(i)**2)/nval(i))
         cvalnow(i)=cvalblk(i)/cwtblk(i)
         cvalblk2(i)=cvalblk2(i)/cwtblk(i)
         cerrnow(i)=sqrt(abs(cvalblk2(i)-cvalnow(i)**2)/nval(i))
      enddo
   endif
   valblk=0.0_r8
   valblk2=0.0_r8
   wtblk=0.0_r8
   nval=0
   cvalblk=0.0_r8
   cvalblk2=0.0_r8
   cwtblk=0.0_r8
   return
   end subroutine update
 
   subroutine result(i,val,err)
   integer(kind=i4) :: i
   real(kind=r8) :: val,err
   val=real(valnow(i))
   err=real(errnow(i))
   return
   end subroutine result

   subroutine cresult(i,cval,cerr)
   integer(kind=i4) :: i
   complex(kind=r8) :: cval
   real(kind=r8) :: cerr
   cval=cvalnow(i)
   cerr=cerrnow(i)
   return
   end subroutine cresult

   function resstring(tau)
   character(len=90), dimension(nest) :: resstring
   integer(kind=i4) :: i
   real(kind=r8) :: val,err
   real(kind=r8) :: tau
   do i=1,nest
      call result(i,val,err)
      write (resstring(i),'(a30,t30,1x,1p,e15.8,e19.10,'' +- '',e19.10)') label(i),tau,val,err
   enddo
   return
   end function resstring

end module estimator
