module density
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntab
   integer(kind=i4), private, save :: nblk
   real(kind=r8), private, save :: dr,scale
   real(kind=r8), private, save, allocatable :: rnow(:),rr(:)
   real(kind=r8), private, save, allocatable :: rerr(:),re(:)
   real(kind=r8), private, save :: wsblk,totalw

contains
   subroutine setuprho(ntabin,npartin,rmax)
   integer(kind=i4) :: ntabin,npartin
   real(kind=r8) :: rmax
   ntab=ntabin
   allocate(rnow(0:ntab),rerr(0:ntab))
   allocate(rr(0:ntab),re(0:ntab))
   npart=npartin
   call zerorho
   dr=rmax/ntab
   scale=1.0_r8/dr
   end subroutine setuprho

   subroutine zerorho
   rnow=0.0_r8
   rr=0.0_r8
   rerr=0.0_r8
   re=0.0_r8
   wsblk=0.0_r8
   nblk=0
   end subroutine zerorho

   subroutine addrho(x,weight)
   real(kind=r8) :: x(3,npart),weight
   integer(kind=i4) :: i,index
   real(kind=r8) :: r,r2
   do i=1,npart
      r2=sum(x(:,i)**2)
      r=sqrt(r2)
      index=scale*r
      index=min(index,ntab)
      rnow(index)=rnow(index)+weight
   enddo
   wsblk=wsblk+weight
   end subroutine addrho

   subroutine updaterho
   use mympi
   real(kind=r8) :: norm
   real(kind=r8) :: rtnow(0:ntab),errtnow(0:ntab),err(0:ntab),wtsblk
   call addall(rnow,rtnow)
   call addall(rerr,errtnow)
   call addall(wsblk,wtsblk)
   if (myrank().eq.0) then
      wsblk=wtsblk
      totalw=wsblk
      norm=1.0_r8/wsblk
      rnow=rtnow
      rerr=errtnow
      rnow=rnow*norm
      rr=rnow
      rerr=rerr*norm
      err=sqrt(abs(rerr-rnow**2)/wsblk)
      re=err
   endif
   wsblk=0.0_r8
   rnow=0.0_r8
   rerr=0.0_r8
   end subroutine updaterho

   subroutine writerho(filename,tau)
   use mympi
   character(len=*) :: filename
   real(kind=r8) :: r0,r1,r
   real(kind=r8) :: vol,pi,tau
   integer(kind=i4) :: i
   if (myrank().ne.0) return
   pi=4.0_r8*atan(1.0_r8)
   open (unit=78,file=filename,position='append')
   write(78,'(''# r (fm), rho (fm**-3) , err'')')
   write(78,'(''# total weight = '',e15.7)') totalw
   write(78,'(''# tau = '',e15.7)') tau
   do i=0,ntab-1
      r0=i*dr
      r1=(i+1)*dr
      r=(i+0.5_r8)*dr
      vol=4.0_r8*pi*(r1**3-r0**3)/3.0_r8
      write (78,'(7e15.7)') r,rr(i)/vol,re(i)/vol
   enddo
   write(78,*) ''
   write(78,*) ''
   close(78)
   end subroutine writerho

   subroutine getradii(x,radii,npart)
   real(kind=r8) :: x(:,:)
   real(kind=r8) :: radii,r2
   integer(kind=i4) :: i,npart
   radii=0.0_r8
   do i=1,npart
      r2=sum(x(:,i)**2)
      radii=radii+r2
   enddo
   radii=sqrt(radii/npart)
!  radii=radii/npart
   end subroutine

end module density
