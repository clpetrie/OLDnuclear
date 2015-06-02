module wavefunction
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   integer(kind=i4), private, save :: npart,nbin(4)
   real(kind=r8), private, save :: el
   real(kind=r8), private, save, allocatable :: ak(:,:),ak2(:)
   real(kind=r8), private, save :: a

contains
   subroutine setpsi(npartin,elin,hbar,dolsin)
   use kshellold
   use mympi
   real(kind=r8) :: elin,rho,hbar
   integer(kind=i4) :: npartin,nshell(4),nsh,nk,i,is
   real(kind=r8), allocatable :: vkin(:)
   real(kind=r8), pointer :: aktemp(:,:),ak2temp(:),vk(:)
   real(kind=r8) :: efg,kf,pi
   logical :: filled,dolsin
   if (myrank().eq.0) then
      read (5,*) rho     !number density
      read (5,*) npart   !particle number
      read (5,*) nbin(1) !number of pup
      read (5,*) nbin(2) !number of pdn
      read (5,*) nbin(3) !number of nup
      read (5,*) nbin(4) !number of ndn
      read (5,*) a 
      write (6,'(''density ='',t40,f10.5)') rho
      el=(npart/rho)**(1.0_r8/3.0_r8)
      write (6,'(''L/2 ='',t40,f10.5)') el*0.5_r8
      pi=4.0_r8*atan(1.0_r8)
      kf=(3.0_r8/2.0_r8*pi**2*rho)**(1.0_r8/3.0_r8)
      efg=3.0_r8/5.0_r8*hbar*kf**2
      write (6,'(''k_F ='',t40,f10.5)') kf
      write (6,'(''E_FG ='',t40,f10.5)') efg
      write (6,'(''npart ='',t40,i10)') npart
      write (6,'(''proton up ='',t40,i10)') nbin(1)
      write (6,'(''proton down ='',t40,i10)') nbin(2)
      write (6,'(''neutron up ='',t40,i10)') nbin(3)
      write (6,'(''neutron down ='',t40,i10)') nbin(4)
   endif
   call bcast(npart)
   call bcast(nbin)
   call bcast(a)
   call bcast(el)
   npartin=npart
   elin=el

!
! ugly k-vector setup -- fix me
!
   do is=1,4
      if (nbin(is).eq.0) cycle
      call shell(nbin(is),nshell(is),filled)
      if (.not.filled) then
         write (6,'(1x,''shell'',i3,'' must be filled - stronzo'')') is
      endif
   enddo
   nsh=maxval(nshell)
   allocate (vkin(nsh)) ! not used here but needs to be allocated for setupk
   vkin=0.0_r8
   call setupk(el,nsh,vkin,aktemp,ak2temp,vk,nk)
   deallocate(vkin)
   if (allocated(ak)) deallocate(ak,ak2)
   allocate(ak(3,maxval(nbin)),ak2(maxval(nbin)))
   ak(:,1)=aktemp(:,1)
   ak2(1)=ak2temp(1)
   do i=2,nk
      ak(:,2*i-2)=aktemp(:,i)
      ak(:,2*i-1)=-aktemp(:,i)
      ak2(2*i-2)=ak2temp(i)
      ak2(2*i-1)=ak2temp(i)
   enddo
   if (associated(aktemp)) deallocate(aktemp)
   if (associated(ak2temp)) deallocate(ak2temp)
   if (associated(vk)) deallocate(vk)
   return
   end subroutine setpsi

   subroutine hpsi(w,doall)
   use stack ! to define walker type
   use v6pot
   use jastrow
   use matrixmod
   type (walker) :: w
   complex(kind=r8) :: spx(4,15,npart),sxall(npart,15,npart)
   complex(kind=r8) :: ph(npart,4,npart),dph(npart,4,npart,3)
   complex(kind=r8) :: sxmall(npart,15,npart)
   complex(kind=r8) :: det,ddet(3,npart)
   complex(kind=r8) :: smati(npart,npart)
   complex(kind=r8) :: d1,d2,d3,d4,cvs,cvt,cvst
   real(kind=r8) :: r,dx(3),uj,ujp,ujpp,ak2tot,akr
   real(kind=r8) :: asig(3,npart,3,npart),asigtau(3,npart,3,npart)
   real(kind=r8) :: atau(npart,npart),du(3,npart),u,d2u
   integer(kind=i4) :: i,j,k,kk,is,ic1,ic2,i1,i2,jc1,jc2,it
   integer(kind=i4) :: ic,jc,kc
   logical :: doall
   w%x=w%x-el*nint(w%x/el)
   spx=opmult(w%sp) ! multiply by sigma, tau, and sigma tau
   ph=czero
   dph=czero
   ak2tot=0.0_r8
   kk=0
   do is=1,4
      ak2tot=ak2tot+sum(ak2(1:nbin(is)))
      do k=1,nbin(is)
         kk=kk+1
         do i=1,npart
            akr=w%x(1,i)*ak(1,k)+w%x(2,i)*ak(2,k)+w%x(3,i)*ak(3,k)
            ph(kk,is,i)=cmplx(cos(akr),-sin(akr))
            dph(kk,is,i,:)=-ci*ak(:,k)*ph(kk,is,i)
         enddo
      enddo
   enddo
   do i=1,npart
      smati(:,i)=matmul(ph(:,:,i),w%sp(:,i))
   enddo
   call cmatinv(smati,det,npart)
   u=0.0_r8
   d2u=0.0_r8
   du=0.0_r8
   do i=2,npart
      do j=1,i-1
         dx=w%x(:,i)-w%x(:,j)
         dx=dx-el*nint(dx/el)
         r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         if (r.lt.el*0.5_r8) then
            call jastro(r,uj,ujp,ujpp)
            u=u+uj
            d2u=d2u+ujpp
            du(:,i)=du(:,i)-ujp*dx
            du(:,j)=du(:,j)+ujp*dx
         endif
      enddo
   enddo
   d2u=2.0_r8*d2u
   w%psi=det*exp(-u)
   w%psig=exp(-u)*sqrt(det**2+a)
   call hspot(w%x,w%vc,asig,asigtau,atau,doall)
   if (doall.eq..false.) return

   do i=1,npart
      do ic=1,3
         w%sigma(ic,i)=sum(smati(i,:)*matmul(ph(:,:,i),spx(:,ic,i)))
         w%tau(ic,i)=sum(smati(i,:)*matmul(ph(:,:,i),spx(:,ic+3,i)))
         do jc=1,3
            kc=jc+3*(ic-1)+6
            w%sigmatau(ic,jc,i)=sum(smati(i,:)*matmul(ph(:,:,i),spx(:,kc,i)))
         enddo
      enddo
   enddo

   do i=1,npart
      ddet(1,i)=sum(smati(i,:)*matmul(dph(:,:,i,1),w%sp(:,i)))
      ddet(2,i)=sum(smati(i,:)*matmul(dph(:,:,i,2),w%sp(:,i)))
      ddet(3,i)=sum(smati(i,:)*matmul(dph(:,:,i,3),w%sp(:,i)))
   enddo
   w%dpsi=du+ddet
   w%d2psi=-d2u-ak2tot+sum((du+2.0_r8*ddet)*du)

   do i=1,npart
      sxall(:,:,i)=matmul(ph(:,:,i),spx(:,:,i))
   enddo
   sxmall=reshape(matmul(smati,reshape(sxall,(/npart,15*npart/))),shape(sxall))
!
! sigma term
!
   cvs=czero
   do ic1=1,3
      do i1=1,npart-1
         d1=sxmall(i1,ic1,i1)
         do ic2=1,3
            do i2=i1+1,npart
               d2=sxmall(i2,ic2,i2)
               d3=sxmall(i1,ic2,i2)
               d4=sxmall(i2,ic1,i1)
               cvs=cvs+(d1*d2-d3*d4)*asig(ic1,i1,ic2,i2)
            enddo
         enddo
      enddo
   enddo
!
! sigma-tau term
!
   cvst=czero
   do it=1,3
      do ic1=1,3
         jc1=3*(ic1-1)+6+it
         do i1=1,npart-1
            d1=sxmall(i1,jc1,i1)
            do ic2=1,3
               jc2=3*(ic2-1)+6+it
               do i2=i1+1,npart
                  d2=sxmall(i2,jc2,i2)
                  d3=sxmall(i1,jc2,i2)
                  d4=sxmall(i2,jc1,i1)
                  cvst=cvst+(d1*d2-d3*d4)*asigtau(ic1,i1,ic2,i2)
               enddo
            enddo
         enddo
      enddo
   enddo
!
! tau term
!
   cvt=czero
   do it=4,6
      do i1=1,npart-1
         d1=sxmall(i1,it,i1)
         do i2=i1+1,npart
            d2=sxmall(i2,it,i2)
            d3=sxmall(i1,it,i2)
            d4=sxmall(i2,it,i1)
            cvt=cvt+(d1*d2-d3*d4)*atau(i1,i2)
         enddo
      enddo
   enddo
   w%v=cvs+cvt+cvst+w%vc
   return
   end subroutine hpsi

   function opmult(sp)
   complex(kind=r8) :: sp(4,npart),opmult(4,15,npart)
!
! multiply by sigma
!
   opmult(1,1,:)=sp(2,:)
   opmult(2,1,:)=sp(1,:)
   opmult(3,1,:)=sp(4,:)
   opmult(4,1,:)=sp(3,:)
   opmult(1,2,:)=-ci*sp(2,:)
   opmult(2,2,:)=ci*sp(1,:)
   opmult(3,2,:)=-ci*sp(4,:)
   opmult(4,2,:)=ci*sp(3,:)
   opmult(1,3,:)=sp(1,:)
   opmult(2,3,:)=-sp(2,:)
   opmult(3,3,:)=sp(3,:)
   opmult(4,3,:)=-sp(4,:)
!
! multiply by tau
!
   opmult(1,4,:)=sp(3,:)
   opmult(2,4,:)=sp(4,:)
   opmult(3,4,:)=sp(1,:)
   opmult(4,4,:)=sp(2,:)
   opmult(1,5,:)=-ci*sp(3,:)
   opmult(2,5,:)=-ci*sp(4,:)
   opmult(3,5,:)=ci*sp(1,:)
   opmult(4,5,:)=ci*sp(2,:)
   opmult(1,6,:)=sp(1,:)
   opmult(2,6,:)=sp(2,:)
   opmult(3,6,:)=-sp(3,:)
   opmult(4,6,:)=-sp(4,:)
!
! multiply by sigma tau
!
   opmult(1,7:13:3,:)=opmult(3,1:3:1,:)
   opmult(2,7:13:3,:)=opmult(4,1:3:1,:)
   opmult(3,7:13:3,:)=opmult(1,1:3:1,:)
   opmult(4,7:13:3,:)=opmult(2,1:3:1,:)
   opmult(1,8:14:3,:)=-ci*opmult(3,1:3:1,:)
   opmult(2,8:14:3,:)=-ci*opmult(4,1:3:1,:)
   opmult(3,8:14:3,:)=ci*opmult(1,1:3:1,:)
   opmult(4,8:14:3,:)=ci*opmult(2,1:3:1,:)
   opmult(1,9:15:3,:)=opmult(1,1:3:1,:)
   opmult(2,9:15:3,:)=opmult(2,1:3:1,:)
   opmult(3,9:15:3,:)=-opmult(3,1:3:1,:)
   opmult(4,9:15:3,:)=-opmult(4,1:3:1,:)
   return
   end function opmult

   subroutine chkder(w,dx,error)
   use stack
   type (walker) w
   real(kind=r8) :: dx,error
   integer(kind=i4) :: i,ic
   complex(kind=r8) :: dnum,d2num,psi,dpsi(3,npart),d2psi,d2n
   call hpsi(w,.true.)
   psi=w%psi
   dpsi=w%dpsi
   d2psi=w%d2psi
   d2num=czero
   do i=1,npart
      do ic=1,3
         w%x(ic,i)=w%x(ic,i)+dx
         call hpsi(w,.true.)
         dnum=w%psi
         d2n=w%dpsi(ic,i)*w%psi
         w%x(ic,i)=w%x(ic,i)-2.0_r8*dx
         call hpsi(w,.true.)
         dnum=dnum-w%psi
         d2n=d2n-w%dpsi(ic,i)*w%psi
         w%x(ic,i)=w%x(ic,i)+dx
         dnum=dnum/(2.0_r8*dx*psi)
         if (abs(dpsi(ic,i)-dnum).gt.error) then
            write (6,'(1x,'' ic i analytic numerical '',2i5,1p,4e14.6)') &
               ic,i,dpsi(ic,i),dnum
         endif
         d2num=d2num+d2n/(2.0_r8*dx*psi)
      enddo
   enddo
   if (abs(d2num-d2psi).gt.error) then
      write (6,'(1x,''d2 analytic numerical '',1p,4e14.6)') d2psi,d2num
   endif
   return
   end subroutine chkder

   subroutine checkj(w,dphi,t2,tz,j2,jz)
   use stack
   type (walker) :: w
   real(kind=r8) :: xsave(3,npart),dphi,c,s,c2,s2
   complex(kind=r8) :: aj,ssave(4,npart),psi,psip,psim
   complex(kind=r8) :: spx(4,15,npart)
   integer(kind=i4) :: ic,jc,kc
   integer(kind=i4), parameter :: levi(2,3) = reshape((/2,3, 3,1, 1,2/),(/2,3/))
   complex(kind=r8) :: t2,tz,j2,jz
   c=cos(dphi)
   s=sin(dphi)
   c2=cos(dphi*0.5_r8)
   s2=sin(dphi*0.5_r8)
   xsave=w%x
   ssave=w%sp
   spx=opmult(w%sp)
   call hpsi(w,.false.)
   psi=w%psi
   aj=czero
   do ic=1,3
      w%sp(:,:)=c2*ssave(:,:)-ci*s2*spx(:,ic+3,:)
      call hpsi(w,.false.)
      psip=w%psi
      w%sp(:,:)=c2*ssave(:,:)+ci*s2*spx(:,ic+3,:)
      call hpsi(w,.false.)
      psim=w%psi
      aj=aj-(psip+psim-2.0_r8*psi)/(dphi**2*psi)
!     if (ic.eq.3) write (6,*) ' tau_z = ',-(psip-psim)/(ci*2.0_r8*dphi*psi)
      if (ic.eq.3) tz=-(psip-psim)/(ci*2.0_r8*dphi*psi)
   enddo
!  write (6,*) 'tau^2 value ',aj
   t2=aj
   aj=czero
   do ic=1,3
      jc=levi(1,ic)
      kc=levi(2,ic)
      w%x(ic,:)=xsave(ic,:)
      w%x(jc,:)=xsave(jc,:)*c-xsave(kc,:)*s
      w%x(kc,:)=xsave(kc,:)*c+xsave(jc,:)*s
      w%sp(:,:)=c2*ssave(:,:)-ci*s2*spx(:,ic,:)
      call hpsi(w,.false.)
      psip=w%psi
      w%x(ic,:)=xsave(ic,:)
      w%x(jc,:)=xsave(jc,:)*c+xsave(kc,:)*s
      w%x(kc,:)=xsave(kc,:)*c-xsave(jc,:)*s
      w%sp(:,:)=c2*ssave(:,:)+ci*s2*spx(:,ic,:)
      call hpsi(w,.false.)
      psim=w%psi
      aj=aj-(psip+psim-2.0_r8*psi)/(dphi**2*psi)
!     if (ic.eq.3) write (6,*) ' jz = ',-(psip-psim)/(ci*2.0_r8*dphi*psi)
      if (ic.eq.3) jz=-(psip-psim)/(ci*2.0_r8*dphi*psi)
   enddo
!  write (6,*) 'j2 value ',aj
   j2=aj
   w%x=xsave
   w%sp=ssave
   end subroutine checkj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the following are subrotuines needed by the optimization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine setdetparam(params)
   real(kind=r8) :: params(:)
   end subroutine setdetparam

   subroutine getdetparam(nparam,params)
   integer(kind=i4) :: nparam
   real(kind=r8), pointer :: params(:)
   nparam=0
   allocate(params(nparam))
   end subroutine getdetparam

   subroutine getderpsi(w,dpsi)
   use stack
   type (walker) :: w
   real(kind=r8) :: dpsi(:)
   end subroutine getderpsi

   subroutine setijas(n) ! use a different jastrow
   integer(kind=i4) :: n
   end subroutine setijas

   subroutine getjasder(w,dpsi)
   use stack
   type (walker) :: w
   real(kind=r8) :: dpsi(:)
   end subroutine getjasder


end module wavefunction
