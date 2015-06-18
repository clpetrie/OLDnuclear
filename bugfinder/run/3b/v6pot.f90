!$Id: v6pot.f90,v 1.5 2013/12/10 21:21:53 nuclear Exp $
module v6pot
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)

   integer(kind=i4), private, save :: npart,ntab,ibox
   real(kind=r8), private, allocatable, save, dimension(:,:) :: vtab,vtabls
   real(kind=r8), private, save :: scale,range,el,eli,vcbox
   real(kind=r8), private, parameter :: amu=0.69953054_r8
   real(kind=r8), private, save :: c3b,a2p3b,a2s3b,a3p3b,ar3b

contains
   subroutine v6potinit(npartin,iboxin,elin,lpot,ntabin,vfact)
   use cheft
   use mympi
   real(kind=r8) :: vfact(8),elin
   integer(kind=i4) :: ntabin,lpot,i,npartin,iboxin,j
   integer(kind=i4) :: ix,iy,iz
   real(kind=r8) :: dr,h2m,h2mcsb,vv(18),vp(12),ww(14),dx(3)
   real(kind=r8) :: rr
   real(kind=r8) :: v0r,v0s,v0t,vr,vs,vt,kr,ks,kt
   real(kind=r8), parameter :: verytiny=1.0e-20_r8
   ibox=iboxin
   npart=npartin
   el=elin
   if (el.gt.0.0_r8) then
      eli=1.0_r8/el
      range=el*sqrt(3.0_r8)*0.5_r8*(2.0_r8*ibox+1.0_r8)
   else
      range=abs(el)
      ibox=0
      el=0.0_r8
      eli=0.0_r8
   endif
   ntab=ntabin
   scale=ntab/range
   dr=1.0_r8/scale
   allocate(vtab(6,0:ntab))
   allocate(vtabls(2,0:ntab))
!
! vcbox is the contribution from our own image. This can be interpreted two
! ways. If this is toroidal boundary conditions, then the particle is its
! own image and sigma.sigma = 3 etc. It seems more physical to assume periodic
! boundary conditions so the images are other particles that happen to have
! the same spin/isospin. In that case sigma.sigma = 1 etc. We use that below.
!
   if (abs(lpot).lt.25) then
! density-dependent version not working anymore, careful!!!
      call setpot(lpot,1.0_r8,1.0_r8,1.0_r8,1.0_r8,1.0_r8,1.0_r8,h2m,h2mcsb)
      do i=0,ntab
         rr=i*dr
         call pot(0,rr,vv,vp,ww)
         vtab(:,i)=vfact(1:6:1)*vv(1:6:1)
         vtabls(:,i)=vfact(7:8:1)*vv(7:8:1)
         if (myrank().eq.0) then
            write(79,'(9f25.15)') rr,(vtab(j,i),j=1,6),(vtabls(j,i),j=1,2)
         endif
      enddo
      vcbox=0.0_r8 !self image contribution
      do ix=-ibox,ibox
         dx(1)=el*ix
         do iy=-ibox,ibox
            dx(2)=el*iy
            do iz=-ibox,ibox
               if (abs(ix)+abs(iy)+abs(iz).eq.0) cycle
               dx(3)=el*iz
               rr=sqrt(sum(dx**2))
               call pot(0,rr,vv,vp,ww)
               vcbox=vcbox+vv(1)+vv(2)+vv(3)+vv(4)
            enddo
         enddo
      enddo
      vcbox=0.5_r8*vcbox
      if (myrank().eq.0) then
         write(6,'(''self image energy = '',t40,f25.10)') vcbox
      endif
   elseif (abs(lpot).eq.25) then
! Minnesota potential:
      v0r=200.0_r8  ! MeV
      v0t=178.0_r8  ! MeV
      v0s=91.85_r8  ! MeV
      kr=1.487_r8  ! fm**-2
      kt=0.639_r8  ! fm**-2
      ks=0.465_r8  ! fm**-2
      do i=0,ntab
         rr=i*dr
         vr=v0r*exp(-kr*rr**2)
         vt=-v0t*exp(-kt*rr**2)
         vs=-v0s*exp(-ks*rr**2)
         vtab(1,i)=3.0_r8/8.0_r8*(vr+0.5_r8*vt+0.5_r8*vs)
         vtab(2,i)=1.0_r8/8.0_r8*(-vr-1.5_r8*vt+0.5_r8*vs)
         vtab(3,i)=1.0_r8/8.0_r8*(-vr+0.5_r8*vt-1.5_r8*vs)
         vtab(4,i)=1.0_r8/8.0_r8*(-vr-0.5_r8*vt-0.5_r8*vs)
         vtab(5,i)=0.0_r8
         vtab(6,i)=0.0_r8
         vtabls(1,i)=0.0_r8
         vtabls(2,i)=0.0_r8
      enddo
      vcbox=0.0_r8 !self image contribution
      do ix=-ibox,ibox
         dx(1)=el*ix
         do iy=-ibox,ibox
            dx(2)=el*iy
            do iz=-ibox,ibox
               if (abs(ix)+abs(iy)+abs(iz).eq.0) cycle
               dx(3)=el*iz
               rr=sqrt(sum(dx**2))
               vr=v0r*exp(-kr*rr**2)
               vt=-v0t*exp(-kt*rr**2)
               vs=-v0s*exp(-ks*rr**2)
               vv(1)=3.0_r8/8.0_r8*(vr+0.5_r8*vt+0.5_r8*vs)
               vv(2)=1.0_r8/8.0_r8*(-vr-1.5_r8*vt+0.5_r8*vs)
               vv(3)=1.0_r8/8.0_r8*(-vr+0.5_r8*vt-1.5_r8*vs)
               vv(4)=1.0_r8/8.0_r8*(-vr-0.5_r8*vt-0.5_r8*vs)
               vcbox=vcbox+vv(1)+vv(2)+vv(3)+vv(4)
            enddo
         enddo
      enddo
      vcbox=0.5_r8*vcbox
      if (myrank().eq.0) then
         write(6,'(''self image energy = '',t40,f25.10)') vcbox
      endif
   elseif (abs(lpot).gt.100) then
      vcbox=0.0_r8 !self image contribution
      do ix=-ibox,ibox
         dx(1)=el*ix
         do iy=-ibox,ibox
            dx(2)=el*iy
            do iz=-ibox,ibox
               if (abs(ix)+abs(iy)+abs(iz).eq.0) cycle
               dx(3)=el*iz
               rr=sqrt(sum(dx**2))
               call cheft_pot(abs(lpot),rr,vv)
               vcbox=vcbox+vv(1)+vv(2)+vv(3)+vv(4)
            enddo
         enddo
      enddo
      vcbox=0.5_r8*vcbox
      if (myrank().eq.0) then
         write(6,'(''self image energy = '',t40,f25.10)') vcbox
      endif
      do i=0,ntab
         rr=i*dr
         call cheft_pot(abs(lpot),rr,vv)
         vtab(:,i)=vfact(1:6:1)*vv(1:6:1)
         if (abs(vv(7)).lt.verytiny) vv(7)=0.0_r8
         vtabls(:,i)=vfact(7:8:1)*vv(7:8:1)
         if (myrank().eq.0) then
            write(79,'(9e25.10)') rr,(vtab(j,i),j=1,6),(vtabls(j,i),j=1,2)
         endif
      enddo
   endif
   end subroutine v6potinit

   subroutine hspot(x,vc,v2,v3,v4,v5,v6,doall,ils)
   real(kind=r8), dimension (3,npart) :: x
   real(kind=r8), dimension(3,npart,3,npart) :: v5,v6,sigma,sigmatau
   real(kind=r8), dimension(npart,npart) :: v2,v3,v4,tau
   real(kind=r8) :: vc
   real(kind=r8) :: dx0(3),dx(3),vv(6),vs(3,3),vst(3,3),r,c1,c2,c3,c4,dr
   integer(kind=i4) :: i,j,index,ix,iy,iz,jc,kc
   logical :: doall,ils
   real(kind=r8) :: gls(3,npart,npart),glsum(3,npart)
   real(kind=r8) :: vlsv,vlst,fac,vlsfac,hbar
   hbar=20.73554_r8
   v2=0.0_r8
   v3=0.0_r8
   v4=0.0_r8
   v5=0.0_r8
   v6=0.0_r8
   vc=0.0_r8
   gls=0.0_r8
   glsum=0.0_r8
   do j=2,npart
      do i=1,j-1
         dx0(:)=x(:,i)-x(:,j)
         dx0=dx0-el*nint(dx0*eli)
         vs=0.0_r8
         vst=0.0_r8
         do ix=-ibox,ibox
            do iy=-ibox,ibox
               do iz=-ibox,ibox
                  dx(1)=dx0(1)+el*ix
                  dx(2)=dx0(2)+el*iy
                  dx(3)=dx0(3)+el*iz
                  r=sqrt(dot_product(dx,dx))
                  dx=dx/r
                  if (r.lt.range) then
                     dr=scale*r
                     index=dr
                     index=max(1,min(index,ntab-2))
                     dr=dr-index
                     c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
                     c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
                     c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
                     c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
                     vv(:)=c1*vtab(:,index-1)+c2*vtab(:,index) &
                        +c3*vtab(:,index+1)+c4*vtab(:,index+2)
                     if (.not.ils) vc=vc+vv(1)
                     if (doall) then

                        vs(1,1)=vs(1,1)+vv(5)*(3.0_r8*dx(1)*dx(1)-1.0_r8)
                        vs(1,2)=vs(1,2)+vv(5)*3.0_r8*dx(1)*dx(2)
                        vs(1,3)=vs(1,3)+vv(5)*3.0_r8*dx(1)*dx(3)
                        vs(2,2)=vs(2,2)+vv(5)*(3.0_r8*dx(2)*dx(2)-1.0_r8)
                        vs(2,3)=vs(2,3)+vv(5)*3.0_r8*dx(2)*dx(3)
                        vs(3,3)=vs(3,3)+vv(5)*(3.0_r8*dx(3)*dx(3)-1.0_r8)

!                       vs(1,1)=vs(1,1)+vv(3)+vv(5)*(3.0_r8*dx(1)*dx(1)-1.0_r8)
!                       vs(1,2)=vs(1,2)+vv(5)*3.0_r8*dx(1)*dx(2)
!                       vs(1,3)=vs(1,3)+vv(5)*3.0_r8*dx(1)*dx(3)
!                       vs(2,2)=vs(2,2)+vv(3)+vv(5)*(3.0_r8*dx(2)*dx(2)-1.0_r8)
!                       vs(2,3)=vs(2,3)+vv(5)*3.0_r8*dx(2)*dx(3)
!                       vs(3,3)=vs(3,3)+vv(3)+vv(5)*(3.0_r8*dx(3)*dx(3)-1.0_r8)

                        vst(1,1)=vst(1,1)+vv(6)*(3.0_r8*dx(1)*dx(1)-1.0_r8)
                        vst(1,2)=vst(1,2)+vv(6)*3.0_r8*dx(1)*dx(2)
                        vst(1,3)=vst(1,3)+vv(6)*3.0_r8*dx(1)*dx(3)
                        vst(2,2)=vst(2,2)+vv(6)*(3.0_r8*dx(2)*dx(2)-1.0_r8)
                        vst(2,3)=vst(2,3)+vv(6)*3.0_r8*dx(2)*dx(3)
                        vst(3,3)=vst(3,3)+vv(6)*(3.0_r8*dx(3)*dx(3)-1.0_r8)

!                       vst(1,1)=vst(1,1)+vv(4)+vv(6)*(3.0_r8*dx(1)*dx(1)-1.0_r8)
!                       vst(1,2)=vst(1,2)+vv(6)*3.0_r8*dx(1)*dx(2)
!                       vst(1,3)=vst(1,3)+vv(6)*3.0_r8*dx(1)*dx(3)
!                       vst(2,2)=vst(2,2)+vv(4)+vv(6)*(3.0_r8*dx(2)*dx(2)-1.0_r8)
!                       vst(2,3)=vst(2,3)+vv(6)*3.0_r8*dx(2)*dx(3)
!                       vst(3,3)=vst(3,3)+vv(4)+vv(6)*(3.0_r8*dx(3)*dx(3)-1.0_r8)

                        v2(i,j)=v2(i,j)+vv(2)
                        v3(i,j)=v3(i,j)+vv(3)
                        v4(i,j)=v4(i,j)+vv(4)
!                       if (ils.and.ix.eq.0.and.iy.eq.0.and.iz.eq.0) then
!                          call vls(r,vlsv,vlst)
!                          gls(:,i,j)=vlsv*r*dx(:)
!                          gls(:,j,i)=-gls(:,i,j)
!                          glsum(:,i)=glsum(:,i)+gls(:,i,j)
!                          glsum(:,j)=glsum(:,j)+gls(:,j,i)
!                          vlsfac=(vlsv*r)**2/(16.0_r8*hbar)
!                          vs(1,1)=vs(1,1)+vlsfac*(-1.0_r8+dx(1)*dx(1))
!                          vs(1,2)=vs(1,2)+vlsfac*dx(1)*dx(2)
!                          vs(1,3)=vs(1,3)+vlsfac*dx(1)*dx(3)
!                          vs(2,2)=vs(2,2)+vlsfac*(-1.0_r8+dx(2)*dx(2))
!                          vs(2,3)=vs(2,3)+vlsfac*dx(2)*dx(3)
!                          vs(3,3)=vs(3,3)+vlsfac*(-1.0_r8+dx(3)*dx(3))
!                          vc=vc-2.0_r8*vlsfac
!                       endif
                     endif
                  endif
               enddo
            enddo
         enddo
         if (doall) then
            vs(2,1)=vs(1,2)
            vs(3,1)=vs(1,3)
            vs(3,2)=vs(2,3)
            vst(2,1)=vst(1,2)
            vst(3,1)=vst(1,3)
            vst(3,2)=vst(2,3)
            v2(j,i)=v2(i,j)
            v3(j,i)=v3(i,j)
            v4(j,i)=v4(i,j)
            v5(:,i,:,j)=vs(:,:)
            v6(:,i,:,j)=vst(:,:)
            v6(:,j,:,i)=vst(:,:)
         endif
      enddo
   enddo
!  if (ils) then
!     vlsfac=-1.0_r8/(32.0_r8*hbar)
!     vc=vc+vlsfac*(sum(glsum**2)-sum(gls**2))
!      do i=2,npart
!        do j=1,i-1
!           fac=vlsfac*(sum(gls(:,j,i)*(glsum(:,j)-glsum(:,i))) &
!              +sum(gls(:,:,j)*gls(:,:,i))-2.0_r8*sum(gls(:,j,i)**2))
!            do jc=1,3
!               v5(jc,j,jc,i)=v5(jc,j,jc,i)+fac
!               do kc=1,3
!                   v5(jc,j,kc,i)=v5(jc,j,kc,i)+vlsfac*( &
!                      -gls(jc,j,i)*(glsum(kc,j)-gls(kc,j,i)) &
!                      -gls(kc,i,j)*(glsum(jc,i)-gls(jc,i,j)) &
!                      -sum(gls(jc,i,:)*gls(kc,j,:)))
!               enddo
!            enddo
!        enddo
!     enddo
!  endif
   do i=1,npart
      v5(1,i,:,1:i-1)=v5(:,1:i-1,1,i)
      v5(2,i,:,1:i-1)=v5(:,1:i-1,2,i)
      v5(3,i,:,1:i-1)=v5(:,1:i-1,3,i)
   enddo
 tau=v2
 sigma=0.0_r8
 sigma(1,:,1,:)=v3
 sigma(2,:,2,:)=v3
 sigma(3,:,3,:)=v3
 sigma=sigma+v5
 sigmatau=0.0_r8
 sigmatau(1,:,1,:)=v4
 sigmatau(2,:,2,:)=v4
 sigmatau(3,:,3,:)=v4
 sigmatau=sigmatau+v6
   end subroutine hspot

   subroutine vls(r,vvls,vvlst)
   real(kind=r8) :: r,vvls,vvlst,dr,c1,c2,c3,c4
   integer(kind=i4) :: index
   if (r.lt.range) then
      dr=scale*r
      index=dr
      index=max(1,min(index,ntab-2))
      dr=dr-index
      c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
      c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
      c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
      c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
      vvls=c1*vtabls(1,index-1)+c2*vtabls(1,index) &
         +c3*vtabls(1,index+1)+c4*vtabls(1,index+2)
      vvlst=c1*vtabls(2,index-1)+c2*vtabls(2,index) &
         +c3*vtabls(2,index+1)+c4*vtabls(2,index+2)
   else
      vvls=0.0_r8
      vvlst=0.0_r8
   endif
   end subroutine vls

   subroutine hs3b(x,vc,sigmatau1,sigmatau2,xpi,xxpi,a2p3bout,a3p3bout)
   real(kind=r8), dimension (3,npart) :: x
   real(kind=r8), dimension(3,npart,3,npart) :: sigmatau1,sigmatau2
   real(kind=r8) :: vc
   real(kind=r8) :: dx(3),r
   integer(kind=i4) :: i,j
   real(kind=r8) :: y2,t2,z2
   real(kind=r8) :: xpi(3,npart,3,npart),xxpi(3,npart,3,npart)
   integer(kind=i4) :: k,jc,kc
   real(kind=r8) :: ex,cut
   real(kind=r8) :: a2p3bout,a3p3bout
   real(kind=r8) :: gr3b(npart),g2s3b(3,npart,npart)
   vc=0.0_r8
   sigmatau1=0.0_r8
   sigmatau2=0.0_r8
   gr3b=0.0_r8
   g2s3b=0.0_r8
   xpi=0.0_r8
   do j=2,npart
      do i=1,j-1
         dx=x(:,i)-x(:,j)
         dx=dx-el*nint(dx/el)
         r=sqrt(sum(dx**2))
         dx=dx/r
         ex=exp(-amu*r)
         cut=1.0_r8-exp(-c3b*r**2)
         y2=cut*ex/(amu*r)
         t2=(1.0_r8+3.0_r8/(r*amu)+3.0_r8/(r*amu)**2)*y2*cut
         z2=amu*r*(y2-t2)/3.0_r8
         xpi(:,i,1,j)=3.0_r8*t2*dx(1)*dx
         xpi(:,i,2,j)=3.0_r8*t2*dx(2)*dx
         xpi(:,i,3,j)=3.0_r8*t2*dx(3)*dx
         xpi(1,i,1,j)=xpi(1,i,1,j)+(y2-t2)
         xpi(2,i,2,j)=xpi(2,i,2,j)+(y2-t2)
         xpi(3,i,3,j)=xpi(3,i,3,j)+(y2-t2)
         xpi(:,j,:,i)=xpi(:,i,:,j)
         g2s3b(:,j,i)=dx*z2
         g2s3b(:,i,j)=-g2s3b(:,j,i)
         gr3b(i)=gr3b(i)+t2**2
         gr3b(j)=gr3b(j)+t2**2
         vc=vc-ar3b*t2**4
      enddo
   enddo
   vc=vc+0.5_r8*ar3b*sum(gr3b**2)
   xxpi=reshape( &
      matmul(reshape(xpi,(/3*npart,3*npart/)), &
      reshape(xpi,(/3*npart,3*npart/))),(/3,npart,3,npart/))
   do j=2,npart
      do k=1,j-1
         do jc=1,3
            do kc=1,3
               sigmatau1(kc,k,jc,j)=sigmatau1(kc,k,jc,j) &
                  +4.0_r8*a2p3b*xxpi(kc,k,jc,j)  
               sigmatau2(kc,k,jc,j)=sigmatau2(kc,k,jc,j) &
                  +a2s3b*sum(g2s3b(jc,:,j)*g2s3b(kc,:,k))
            enddo
         enddo
         sigmatau1(:,j,:,k)=sigmatau1(:,k,:,j)
         sigmatau2(:,j,:,k)=sigmatau2(:,k,:,j)
      enddo
   enddo
   a2p3bout=a2p3b
   a3p3bout=a3p3b
   end subroutine
   
   subroutine v3_init(v3b,c,a2p,a2s,a3,ar,v3bflag)
   real(kind=r8) :: c,a2p,a2s,a3,ar
   character(len=3) :: v3b
   integer(kind=i4) :: v3bflag
   if(v3b.eq.'no3') then
      c3b=0.0_r8
      a2p3b=0.0_r8
      a2s3b=0.0_r8
      a3p3b=0.0_r8
      ar3b=0.0_r8
      v3bflag=0
   else if(v3b.eq.'uix') then
      c3b=2.1_r8
      a2p3b=-0.0293_r8
      a2s3b=0.0_r8
      a3p3b=0.0_r8
      ar3b=0.00480_r8
      v3bflag=1
   else if(v3b.eq.'il1') then
      c3b=2.1_r8
      a2p3b=-0.0385_r8
      a2s3b=0.1_r8
      a3p3b=0.0026_r8
      ar3b=0.00705_r8
      v3bflag=2
   else if(v3b.eq.'il2') then
      c3b=2.1_r8
      a2p3b=-0.0370_r8
      a2s3b=-1.0_r8
      a3p3b=0.0026_r8
      ar3b=0.00705_r8
      v3bflag=2
   else if(v3b.eq.'il3') then
      c3b=1.5_r8
      a2p3b=-0.070_r8
      a2s3b=-1.0_r8
      a3p3b=0.0065_r8
      ar3b=0.032_r8
      v3bflag=2
   else if(v3b.eq.'il4')then
      c3b=2.1_r8
      a2p3b=-0.0280_r8
      a2s3b=-1.0_r8
      a3p3b=0.0021_r8
      ar3b=0.0039_r8
      v3bflag=2
   else
      write (6,'(1x,''serious error pot3bpar'')')
      write (6,'(1x,''bad three body potential name  '',a3)') v3b
      write (6,'(1x,''allowed names: '')')
      write (6,'(16x,a3)') 'no3'
      write (6,'(16x,a3)') 'uix'
      write (6,'(16x,a3)') 'il1'
      write (6,'(16x,a3)') 'il2'
      write (6,'(16x,a3)') 'il3'
      write (6,'(16x,a3)') 'il4'
      write (6,'(1x,''stopping the calculation '')')
      stop
   endif
   c=c3b
   a2p=a2p3b
   a2s=a2s3b
   a3=a3p3b
   ar=ar3b
   end subroutine v3_init
   
   subroutine v3_param(c3bin,a2pin,a2sin,a3in,arin)
   real(kind=r8) :: c3bin,a2pin,a2sin,a3in,arin
   c3b=c3bin
   a2p3b=a2pin
   a2s3b=a2sin
   a3p3b=a3in
   ar3b=arin
   end subroutine   
end module v6pot
