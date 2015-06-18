module gofr
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntab,npair,nprot,nneut
   real(kind=r8), private, save :: el,dr,scale,wsblk,totalw,eli
   real(kind=r8), private, save, allocatable :: gnow(:,:),g(:,:)
   real(kind=r8), private, save, allocatable :: gerr(:,:),ge(:,:)
   real(kind=r8), private, save, allocatable :: dnow(:,:),d(:,:)
   real(kind=r8), private, save, allocatable :: derr(:,:),de(:,:)
contains
   subroutine setupgofr(npartin,range,nprotin,nneutin)
   integer(kind=i4) :: npartin,nprotin,nneutin
   real(kind=r8) :: range
   ntab=200 ! number of table points
   npart=npartin
   if (range.eq.0.0_r8) then
      range=20.0_r8
      eli=0.0_r8
   else 
      eli=1.0_r8/range
   endif
   el=range
   npair=npart*(npart-1)/2
   allocate(gnow(10,0:ntab),gerr(10,0:ntab))
   allocate(g(10,0:ntab),ge(10,0:ntab))
   allocate(dnow(4,0:ntab),derr(4,0:ntab))
   allocate(d(4,0:ntab),de(4,0:ntab))
   call zerogofr
   dr=0.5_r8*el/ntab
   scale=1.0_r8/dr
   nprot=nprotin
   nneut=nneutin
   end subroutine setupgofr

   subroutine zerogofr
   g=0.0_r8
   gnow=0.0_r8
   gerr=0.0_r8
   ge=0.0_r8
   wsblk=0.0_r8
   d=0.0_r8
   dnow=0.0_r8
   derr=0.0_r8
   de=0.0_r8
   end subroutine zerogofr

   subroutine addgofr(w,weight)
   use stack
   use wavefunction
   real(kind=r8) :: weight
   integer(kind=i4) :: i,j,index,ij,it,ic
   real(kind=r8) :: dx(3),r
   real(kind=r8) :: xcm(3)
   complex(kind=r8), dimension(3,npart) :: sig1,tau1
   complex(kind=r8), dimension(3,3,npart) :: sigtau1
   complex(kind=r8), dimension(3,3,npair) :: sig2,tau2
   complex(kind=r8), dimension(3,3,3,3,npair) :: sigtau2
   complex(kind=r8), dimension(npair) :: tens,tenstau
   complex(kind=r8), dimension(npair) :: np0,np1,pp,nn
   complex(kind=r8) :: ss,tt,st
   type(walker) :: w
   call getsigtauop(tau1,sig1,sigtau1,tau2,sig2,sigtau2,np0,np1,pp,nn,tens,tenstau)
   ij=0
   do i=1,npart-1
      do j=i+1,npart
         ij=ij+1
         dx(:)=w%x(:,i)-w%x(:,j)
         dx=dx-el*nint(dx*eli)
         r=sqrt(sum(dx**2))
         index=scale*r
         index=min(index,ntab-1)
         ss=0.0_r8
         tt=0.0_r8
         st=0.0_r8
         do ic=1,3
            ss=ss+sig2(ic,ic,ij) ! spin dot product
            tt=tt+tau2(ic,ic,ij) ! tau dot product
            do it=1,3
               st=st+sigtau2(ic,ic,it,it,ij)
            enddo
         enddo
         gnow(1,index)=gnow(1,index)+2.0_r8*weight
         gnow(2,index)=gnow(2,index)+2.0_r8*tt*weight
         gnow(3,index)=gnow(3,index)+2.0_r8*ss*weight
         gnow(4,index)=gnow(4,index)+2.0_r8*st*weight
         gnow(5,index)=gnow(5,index)+2.0_r8*tens(ij)*weight
         gnow(6,index)=gnow(6,index)+2.0_r8*tenstau(ij)*weight
         gnow(7,index)=gnow(7,index)+2.0_r8*np0(ij)*weight
         gnow(8,index)=gnow(8,index)+2.0_r8*np1(ij)*weight
         gnow(9,index)=gnow(9,index)+2.0_r8*pp(ij)*weight
         gnow(10,index)=gnow(10,index)+2.0_r8*nn(ij)*weight
         gerr(1,index)=gerr(1,index)+4.0_r8*weight
         gerr(2,index)=gerr(2,index)+4.0_r8*dconjg(tt)*tt*weight
         gerr(3,index)=gerr(3,index)+4.0_r8*dconjg(ss)*ss*weight
         gerr(4,index)=gerr(4,index)+4.0_r8*dconjg(st)*st*weight
         gerr(5,index)=gerr(5,index)+4.0_r8*dconjg(tens(ij))*tens(ij)*weight
         gerr(6,index)=gerr(6,index)+4.0_r8*dconjg(tenstau(ij))*tenstau(ij)*weight
         gerr(7,index)=gerr(7,index)+4.0_r8*dconjg(np0(ij))*np0(ij)*weight
         gerr(8,index)=gerr(8,index)+4.0_r8*dconjg(np1(ij))*np1(ij)*weight
         gerr(9,index)=gerr(9,index)+4.0_r8*dconjg(pp(ij))*pp(ij)*weight
         gerr(10,index)=gerr(10,index)+4.0_r8*dconjg(nn(ij))*nn(ij)*weight
      enddo
   enddo
   wsblk=wsblk+weight
   do i=1,3
      xcm(i)=sum(w%x(i,:))/npart
   enddo
!  xcm=0.0_r8
   do i=1,npart
      dx=w%x(:,i)-xcm
      dx=dx-el*nint(dx*eli)
      r=sqrt(sum(dx**2))
      index=scale*r
      index=min(index,ntab)
      dnow(1,index)=dnow(1,index)+(1.0_r8+tau1(3,i))*(1.0_r8+sig1(3,i))/4.0_r8*weight
      dnow(2,index)=dnow(2,index)+(1.0_r8+tau1(3,i))*(1.0_r8-sig1(3,i))/4.0_r8*weight
      dnow(3,index)=dnow(3,index)+(1.0_r8-tau1(3,i))*(1.0_r8+sig1(3,i))/4.0_r8*weight
      dnow(4,index)=dnow(4,index)+(1.0_r8-tau1(3,i))*(1.0_r8-sig1(3,i))/4.0_r8*weight
      derr(1,index)=derr(1,index)+((1.0_r8+tau1(3,i))*(1.0_r8+sig1(3,i))/4.0_r8)**2*weight
      derr(2,index)=derr(2,index)+((1.0_r8+tau1(3,i))*(1.0_r8-sig1(3,i))/4.0_r8)**2*weight
      derr(3,index)=derr(3,index)+((1.0_r8-tau1(3,i))*(1.0_r8+sig1(3,i))/4.0_r8)**2*weight
      derr(4,index)=derr(4,index)+((1.0_r8-tau1(3,i))*(1.0_r8-sig1(3,i))/4.0_r8)**2*weight
   enddo
   end subroutine addgofr

   subroutine updategofr
   use mympi
   real(kind=r8) :: norm,wtsblk
   real(kind=r8) :: gtnow(10,0:ntab),errtnow(10,0:ntab),err(10,0:ntab)
   real(kind=r8) :: dtnow(4,0:ntab),errdtnow(4,0:ntab),errd(4,0:ntab)
   call addall(gnow,gtnow)
   call addall(gerr,errtnow)
   call addall(wsblk,wtsblk)
   call addall(dnow,dtnow)
   call addall(derr,errdtnow)
   if (myrank().eq.0) then
      wsblk=wtsblk
      totalw=wsblk
      norm=1.0_r8/wsblk
      gnow=gtnow
      gerr=errtnow
      gnow=gnow*norm
      g=gnow
      gerr=gerr*norm
      err=sqrt(abs(gerr-gnow**2)/wsblk)
      ge=err
      dnow=dtnow
      derr=errdtnow
      dnow=dnow*norm
      d=dnow
      derr=derr*norm
      errd=sqrt(abs(derr-dnow**2)/wsblk)
      de=errd
   endif
   wsblk=0.0_r8
   gnow=0.0_r8
   gerr=0.0_r8
   dnow=0.0_r8
   derr=0.0_r8
   end subroutine updategofr

   subroutine writegofr(tau,filext)
   use mympi
   real(kind=r8) :: pi,r0,r1,facl,r
   real(kind=r8) :: vol,rho,tau
   integer(kind=i4) :: i
   character(len=*) :: filext
   call updategofr
   if (myrank().ne.0) return
   pi=4.0_r8*atan(1.0_r8)
   open (unit=77,file='gofr6.'//trim(filext),position='append')
   open (unit=78,file='gofrnp.'//trim(filext),position='append')
   write(77,'(''# r, gcen(r), err, gtau(r), err, gsig(r), err, gsigtau(r), err, gtens(r), err, gtenstau(r), err'')')
   write(77,'(''# total weight = '',e15.7)') totalw
   write(77,'(''# tau = '',e15.7)') tau
   write(78,'(''# r, gnp0(r), err, gnp1(r), err, gpp(r), err, gnn(r), err'')')
   write(78,'(''# total weight = '',e15.7)') totalw
   write(78,'(''# tau = '',e15.7)') tau
   rho=(npart/el**3)
   do i=0,ntab-1
      r0=i*dr
      r1=(i+1)*dr
      r=(i+0.5_r8)*dr
      vol=4.0_r8*pi*(r1**3-r0**3)/3.0_r8
!     facl=1.0_r8/(rho*npart*vol)
      facl=1.0_r8/vol
      g(:,i)=g(:,i)*facl
      ge(:,i)=ge(:,i)*facl
      g(1:6,i)=g(1:6,i)/(npart*(npart-1)/2)
      ge(1:6,i)=ge(1:6,i)/(npart*(npart-1)/2)
      if (nprot.ne.0.and.nneut.ne.0) then
         g(7,i)=g(7,i)/(nprot*nneut)
         ge(7,i)=ge(7,i)/(nprot*nneut)
         g(8,i)=g(8,i)/(nprot*nneut)
         ge(8,i)=ge(8,i)/(nprot*nneut)
      endif
      if (nprot.ne.0) then
         g(9,i)=g(9,i)/(0.5_r8*nprot*(nprot-1))
         ge(9,i)=ge(9,i)/(0.5_r8*nprot*(nprot-1))
      endif
      if (nneut.ne.0) then
         g(10,i)=g(10,i)/(0.5_r8*nneut*(nneut-1))
         ge(10,i)=ge(10,i)/(0.5_r8*nneut*(nneut-1))
      endif
      write (77,'(13e15.7)') r,g(1,i),ge(1,i),g(2,i),ge(2,i),g(3,i),ge(3,i),g(4,i),ge(4,i),g(5,i),ge(5,i),g(6,i),ge(6,i)
      write (78,'(9e15.7)') r,g(7,i),ge(7,i),g(8,i),ge(8,i),g(9,i),ge(9,i),g(10,i),ge(10,i)
   enddo
   write(77,*) ''
   write(77,*) ''
   close(77)
   write(78,*) ''
   write(78,*) ''
   close(78)
   open (unit=77,file='rho.'//trim(filext),position='append')
   write(77,'(''# r, rhopup, err, rhopdn, err, rhonup, err, rhondn, err'')')
   write(77,'(''# total weight = '',e15.7)') totalw
   write(77,'(''# tau = '',e15.7)') tau
   do i=0,ntab
      r0=i*dr
      r1=(i+1)*dr
      r=(i+0.5_r8)*dr
      vol=4.0_r8*pi*(r1**3-r0**3)/3.0_r8
      facl=1.0_r8/vol
      d(:,i)=d(:,i)*facl
      de(:,i)=de(:,i)*facl
      write (77,'(9e15.7)') r,d(1,i),de(1,i),d(2,i),de(2,i),d(3,i),de(3,i),d(4,i),de(4,i)
   enddo
   write(77,*) ''
   write(77,*) ''
   close(77)
   end subroutine writegofr

   subroutine getradii(w,rad,radp,radn,npart)
   use stack
   use wavefunction
   complex(kind=r8), dimension(3,npart) :: sig1,tau1
   complex(kind=r8), dimension(3,3,npart) :: sigtau1
   complex(kind=r8), dimension(3,3,npair) :: sig2,tau2
   complex(kind=r8), dimension(3,3,3,3,npair) :: sigtau2
   complex(kind=r8), dimension(npair) :: tens,tenstau
   complex(kind=r8), dimension(npair) :: np0,np1,pp,nn
   type(walker) :: w
   real(kind=r8) :: rad,radp,radn,r2
   integer(kind=i4) :: i,npart
   call getsigtauop(tau1,sig1,sigtau1,tau2,sig2,sigtau2,np0,np1,pp,nn,tens,tenstau)
   rad=0.0_r8
   radp=0.0_r8
   radn=0.0_r8
   do i=1,npart
      r2=sum(w%x(:,i)**2)
      rad=rad+r2
      radp=radp+r2*(1+tau1(3,i))/2
      radn=radn+r2*(1-tau1(3,i))/2
   enddo
   radp=max(0.0_r8,radp)
   radn=max(0.0_r8,radn)
   rad=rad/npart
   if (nprot.ne.0) radp=radp/nprot
   if (nneut.ne.0) radn=radn/nneut
   end subroutine getradii
end module gofr
