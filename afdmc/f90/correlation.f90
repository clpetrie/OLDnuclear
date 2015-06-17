module correlation
  implicit none
  integer, private, parameter :: i4=selected_int_kind(9)
  integer, private, parameter :: i8=selected_int_kind(15)
  integer, private, parameter :: r8=selected_real_kind(15,9)
  real(kind=r8), private, save, allocatable :: f(:,:),df(:,:),d2f(:,:)
  real(kind=r8), private, save :: dr,dr2
  integer(kind=i4), private, save :: notab
contains
   subroutine corr(r,fc,dfc,d2fc)
   real(kind=r8) :: r,xi,xii,yi,yii,a,b,c,d,ddyi,ddyii
   real(kind=r8) :: fc(:),dfc(:),d2fc(:),r0
   integer(kind=i4) :: m,idx
   do m=1,6
      idx=nint(r/dr)
      xi=idx*dr
      xii=xi+dr
      a=(xii-r)/dr
      b=1.0_r8-a
      c=(a**3-a)*dr2/6.0_r8
      d=(b**3-b)*dr2/6.0_r8
      yi=f(m,idx)
      yii=f(m,idx+1)
      ddyi=d2f(m,idx)
      ddyii=d2f(m,idx+1)
      fc(m)=a*yi+b*yii+c*ddyi+d*ddyii
      dfc(m)=(yii-yi)/dr-(3.0_r8*a**2-1.0_r8)/6.0_r8* &
        dr*ddyi+(3.0_r8*b**2-1.0_r8)/6.0_r8*dr*ddyii
      d2fc(m)=a*ddyi+b*ddyii
   enddo
   end subroutine corr
    
   subroutine readcorr
   use mympi
   use jastrow
   integer(kind=i4) :: i,j,m,n,t
   real(kind=r8), allocatable :: x(:)
   real(kind=r8), allocatable :: orb(:,:),dorb(:,:),d2orb(:,:)
   real(kind=r8) :: tmp,r
   integer(kind=i4) :: itmp
   real(kind=r8) :: alpha,const,const2,xm
   if (myrank().eq.0) then
      open(unit=21,file='corr.deut',status='old')
!     open(unit=21,file='corr.he4',status='old')
      read (21,*) notab
      allocate(orb(6,notab))
      do i=1,notab
         read (21,*) dr,orb(1,i),orb(6,i)
! this is a test to make the deuteron working!
         orb(6,i)=-orb(6,i)/sqrt(8.0_r8)/3.0_r8

!        read (21,*) dr,orb(:,i)
         orb(2,i)=orb(2,i)/orb(1,i)
         orb(3,i)=orb(3,i)/orb(1,i)
         orb(4,i)=orb(4,i)/orb(1,i)
         orb(5,i)=orb(5,i)/orb(1,i)
         orb(6,i)=orb(6,i)/orb(1,i)
!orb(2,i)=0.0_r8
!orb(3,i)=0.0_r8
!orb(4,i)=0.0_r8
!orb(5,i)=0.0_r8
!orb(6,i)=0.5_r8*orb(6,i)
         orb(1,i)=-log(orb(1,i))
      enddo
      close(21)
   endif
   call bcast(notab)
   if (myrank().ne.0) allocate(orb(6,notab))
   call bcast(orb)
   allocate(dorb(6,notab),d2orb(6,notab))
   allocate(f(6,notab),df(6,notab),d2f(6,notab))
   allocate(x(notab))

   dr=0.004_r8  ! this is a test to make the deuteron working!
!  dr=0.01_r8
   dr2=dr**2

   n=notab

   do m=1,6
      do i=1,n
         x(i)=orb(m,i)
      enddo 
      dorb(m,1)=(-3*x(5)+16*x(4)-36*x(3)+48*x(2)-25*x(1))/(12*dr)
      dorb(m,2)=(-3*x(6)+16*x(5)-36*x(4)+48*x(3)-25*x(2))/(12*dr)
      dorb(m,3)=(-3*x(7)+16*x(6)-36*x(5)+48*x(4)-25*x(3))/(12*dr)
  
      d2orb(m,1)=(11*x(5)-56*x(4)+114*x(3)-104*x(2)+35*x(1))/(12*dr*dr)
      d2orb(m,2)=(11*x(6)-56*x(5)+114*x(4)-104*x(3)+35*x(2))/(12*dr*dr)
      d2orb(m,3)=(11*x(7)-56*x(6)+114*x(5)-104*x(4)+35*x(3))/(12*dr*dr)
     
      do i=4,n-3
         dorb(m,i)=(x(i+3)-9*x(i+2)+45*x(i+1)-45*x(i-1)+9*x(i-2)-x(i-3))/(60*dr)
         d2orb(m,i)=(2*x(i+3)-27*x(i+2)+270*x(i+1)-490*x(i)+270*x(i-1)-27*&
            x(i-2)+2*x(i-3))/(180*dr*dr)
      enddo
  
      dorb(m,n-2)=(25*x(n-2)-48*x(n-3)+36*x(n-4)-16*x(n-5)+3*x(n-6))/(12*dr)
      dorb(m,n-1)=(25*x(n-1)-48*x(n-2)+36*x(n-3)-16*x(n-4)+3*x(n-5))/(12*dr)
      dorb(m,n)=(25*x(n)-48*x(n-1)+36*x(n-2)-16*x(n-3)+3*x(n-4))/(12*dr)
      d2orb(m,n-2)=(35*x(n-2)-104*x(n-3)+114*x(n-4)-56*x(n-5)+11*x(n-6))/(dr*dr) 
      d2orb(m,n-1)=(35*x(n-1)-104*x(n-2)+114*x(n-3)-56*x(n-4)+11*x(n-5))/(dr*dr)
      d2orb(m,n)=(35*x(n)-104*x(n-1)+114*x(n-2)-56*x(n-3)+11*x(n-4))/(dr*dr)
   enddo

!  call setjtab(dr,orb(1,:),dorb(1,:),d2orb(1,:),notab)
   orb(1,:)=1.0_r8
   dorb(1,:)=0.0_r8
   d2orb(1,:)=0.0_r8

   f=orb
   df=dorb
   d2f=d2orb

   do i=1,notab
      r=i*dr
      write(41,'(4f25.15)') r,f(1,i),df(1,i),d2f(1,i)
      write(42,'(4f25.15)') r,f(2,i),df(2,i),d2f(2,i)
      write(43,'(4f25.15)') r,f(3,i),df(3,i),d2f(3,i)
      write(44,'(4f25.15)') r,f(4,i),df(4,i),d2f(4,i)
      write(45,'(4f25.15)') r,f(5,i),df(5,i),d2f(5,i)
      write(46,'(4f25.15)') r,f(6,i),df(6,i),d2f(6,i)
   enddo


   return
   end subroutine

end module    
