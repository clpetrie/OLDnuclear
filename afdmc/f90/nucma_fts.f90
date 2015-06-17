module nucma
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
contains
!*id* nmfts ***********************************************************
! subroutine for finding correlation functions
! in spin-isospin channels for v14 problem
! **********************************************************************
  subroutine nmfts(lcx,lsx,ltx,llx,acn,at,as,ast,atn,att,als,dt,nm,rho,f, &
     &                 fp,fpp,ngrid,h,lpot,h2m)
      use mympi
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
!
      integer*4, parameter :: nout=6
      real*8 :: kf,rho,acn,ast,atn,als,dt,h2m,pi,at,as,att
      real*8 :: rx(ngrid+2),slx(ngrid+2),slpx(ngrid+2),sldpx(ngrid+2),sltpx(ngrid+2)
      real*8 :: rr,vv(18),ww(14),rv(6)
      real*8 :: rsx(ngrid+2),slsx(ngrid+2),slpsx(ngrid+2),chi(ngrid+2),phir(10,ngrid+2),&
    & pm(8,ngrid+2),psi(8,ngrid+2),vx(14,ngrid+2),rlm(8,ngrid+2),rlx(2,ngrid+2),c(3,8),ca(3,8),blm(8)
      dimension f(1:ngrid,8),fp(1:ngrid,8),fpp(1:ngrid,8)
      save
!
      pi=acos(-1.)
      kf=(1.5*nm*pi**2*rho)**(1./3.)
!     
      rt2=sqrt(2.)
      rt5=sqrt(5.)
      ltxp=ltx+1
      ltxm=ltx-1
      lcxp=lcx+1
      lcxm=lcx-1
      llxp=llx+1
      llxm=llx-1
!      h=dt/dble(ltx)
      h2=h*h/12/h2m
!   --------------------
!   set up r,sl,phi,etc.
!   --------------------
      do i=1,ltxp
       rx(i)=h*float(i)
       rsx(i)=rx(i)*rx(i)
       x=kf*rx(i)
       xx=x*x
       if (x.lt.0.125) then
        slx(i)=1.-xx/10.+xx**2/280.
        slpx(i)=(-xx/5.+xx**2/70.-xx**3/2520.)/rx(i)
        sldpx(i)=kf**2*(-1./5.+3.*xx/70.-xx**2/504.)
        sltpx(i)=kf**2*(3.*xx/35.-xx**2/126.+xx**3/3960.)/rx(i)
       else
        y=sin(x)
        z=cos(x)
        slx(i)=3.*(y-x*z)/x**3
        slpx(i)=3.*kf*(-slx(i)+y/x)/x
        sldpx(i)=3.*(-slpx(i)+(slx(i)+z-2.*y/x)/rx(i))/rx(i)
        sltpx(i)=3.*(-sldpx(i)*rx(i)+2.*slpx(i) &
     &  -2.*(slx(i)+z-2.*y/x)/rx(i)-kf*(y+2.*z/x-2.*y/xx))/rsx(i)
       endif
       slsx(i)=slx(i)*slx(i)
       slpsx(i)=slpx(i)**2
       rllp=rx(i)*slx(i)*slpx(i)
       rdls=rllp+rsx(i)*(slpsx(i)+slx(i)*sldpx(i))
       phir(1,i)=sqrt(1.+slsx(i))*rx(i)
       phir(3,i)=sqrt(1.-slsx(i))*rx(i)
       phir(5,i)=phir(3,i)
       phir(7,i)=sqrt((1./5.)*xx-rllp)*rx(i)
       phir(8,i)=sqrt((1./5.)*xx+rllp)*rx(i)
       phir(9,i)=sqrt((12./175.)*xx**2+(2./5.)*xx-rdls)*rx(i)
       pm(1,i)=h2m*(slsx(i)*slpsx(i)/(1.+slsx(i))-slpsx(i) &
     &  -slx(i)*sldpx(i)-2.*slx(i)*slpx(i)/rx(i))/(1.+slsx(i))
       pm(3,i)=h2m*(slsx(i)*slpsx(i)/(1.-slsx(i))+slpsx(i) &
     &  +slx(i)*sldpx(i)+2.*slx(i)*slpx(i)/rx(i))/(1.-slsx(i))
       pm(5,i)=pm(3,i)
       phip=((1./5.)*xx-.5*rdls)/phir(7,i)
       pm(7,i)=h2m*((phip**2*rsx(i)-(1./5.)*xx+rsx(i)*(slpsx(i) &
     &  +slx(i)*sldpx(i)+(1.5*slpx(i)*sldpx(i)+.5*slx(i)*sltpx(i)) &
     &  *rx(i)))/phir(7,i)**2-2.*phip/phir(7,i))
       psi(1,i)=phir(1,i)
       psi(3,i)=phir(3,i)
       psi(5,i)=0.
       psi(7,i)=0.
       rr=rx(i)
       call vpot_op(lpot,0,rr,vv,ww)
       vx(1,i)=acn*vv(1)+at*vv(2)-3.*(as*vv(3)+ast*vv(4))
       vx(3,i)=acn*vv(1)+at*vv(2)+as*vv(3)+ast*vv(4)
       vx(5,i)=atn*vv(5)+att*vv(6)
       vx(7,i)=als*(vv(7)+vv(8))
       vx(9,i)=0.0d0
       vx(10,i)=0.0d0
       vx(11,i)=0.0d0
       vx(12,i)=0.0d0
       vx(13,i)=0.0d0
       vx(14,i)=0.0d0
       if (nm.eq.1) then
        phir(2,i)=phir(3,i)
        phir(4,i)=phir(1,i)
        phir(6,i)=phir(1,i)
        phir(10,i)=sqrt((12./175.)*xx**2+(2./5.)*xx+rdls)*rx(i)
        pm(2,i)=pm(3,i)
        pm(4,i)=pm(1,i)
        pm(6,i)=pm(1,i)
        phip=((1./5.)*xx+.5*rdls)/phir(8,i)
        pm(8,i)=h2m*((phip**2*rsx(i)-(1./5.)*xx-rsx(i)*(slpsx(i) &
     &   +slx(i)*sldpx(i)+(1.5*slpx(i)*sldpx(i)+.5*slx(i)*sltpx(i)) &
     &   *rx(i)))/phir(8,i)**2-2.*phip/phir(8,i))
        psi(2,i)=phir(2,i)
        psi(4,i)=phir(4,i)
        psi(6,i)=0.
        psi(8,i)=0.
        vx(2,i)=acn*vv(1)-3.*(at*vv(2)+as*vv(3)-3.*ast*vv(4))
        vx(4,i)=acn*vv(1)+as*vv(3)-3.*(at*vv(2)+ast*vv(4))
        vx(6,i)=atn*vv(5)-3.*att*vv(6)
        vx(8,i)=als*(vv(7)-3.*vv(8))
        rlm(2,i)=vx(2,i)
        rlx(2,i)=vx(4,i)
       end if
       rlm(1,i)=vx(1,i)
       rlx(1,i)=vx(3,i)
      enddo
      rr=0.
      call vpot_op(lpot,1,rr,vv,ww)


      if (nm.eq.1) then
       rv(1)=acn*vv(1)+at*vv(2)-3*(as*vv(3)+ast*vv(4))
       rv(4)=acn*vv(1)+as*vv(3)-3*(at*vv(2)+ast*vv(4))
       rv(6)=atn*vv(5)-3.*att*vv(6)
      else if (nm.eq.2) then
       rv(1)=acn*vv(1)+at*vv(2)-3.*(as*vv(3)+ast*vv(4))
      end if
!   ---------------------------
!   single-channel psi equation
!   ---------------------------
      do 190 k=1,2,nm
        kp=9-k
        kq=k+8
        blm(k)=0.
        small=1.e-7
        do 140 loop=1,20
          do 110 i=1,lcx
            rlm(k,i)=blm(k)
  110     continue
          psim=0.
          fcm=3.*psi(k,1)/phir(k,1)-3.*psi(k,2)/phir(k,2) &
     &         +psi(k,3)/phir(k,3)
          if (k.eq.1) gpsim=h2*rt2*rv(k)*fcm
          if (k.eq.2) gpsim=h2*h2m*2*kf*fcm/rt5
          psi0=psi(k,1)
          gpsi0=h2*(vx(k,1)-pm(k,1)-blm(k)+vx(kq,1)*(phir(kp,1)/phir(k,1))**2)*psi0
          do 120 j=2,lcxp
            gp=h2*(vx(k,j)-pm(k,j)-blm(k)+vx(kq,j)*(phir(kp,j)/phir(k,j))**2)
            psi(k,j)=(2.*psi0+10.*gpsi0-psim+gpsim)/(1.-gp)
            psim=psi0
            gpsim=gpsi0
            psi0=psi(k,j)
            gpsi0=gp*psi0
  120     continue
          dldif=phir(k,lcx)*(psi(k,lcxp)-psi(k,lcxm)) &
     &         -psi(k,lcx)*(phir(k,lcxp)-phir(k,lcxm))
          if (loop.eq.1) then
            dldifo=dldif
            blmo=blm(k)
            blm(k)=(-1.)**k
          else
            if (abs(dldifo-dldif).le.small) go to 170
            blmn=(dldifo*blm(k)-dldif*blmo)/(dldifo-dldif)
            dldifo=dldif
            blmo=blm(k)
            blm(k)=blmn
          end if
  140   continue
  170   fac=phir(k,lcx)/psi(k,lcx)
        do 180 j=1,lcx
          psi(k,j)=fac*psi(k,j)
  180   continue
        psi(k,lcxp)=phir(k,lcxp)
  190 continue


! -----------------------------
! coupled-channel psi equations
! -----------------------------
      do 510 l=3,4,nm
        m=l+2
        n=l+4
        k=l-2
        kqq=l+6
        kq=l+8
        kbb=l+10
        blm(l)=0.
        blm(m)=0.
        blm(n)=0.
        do 201 i=1,lcx
          rlm(l,i)=0.
  201   continue
        do 202 i=lcxp,ltxp
          rlm(l,i)=rlx(l-2,i)
  202   continue
        do 203 i=1,ltx
          rlm(m,i)=0.
  203   continue
        pqpcs=(phir(n,ltxp)/phir(l,ltxp))**2
        rlm(m,ltxp)=vx(m,ltxp)-(1./12.)*pqpcs*vx(kbb,ltxp)
        do 204 i=1,llx
          rlm(n,i)=0.
  204   continue
        do 205 i=llxp,ltxp
          rlm(n,i)=vx(n,i)-.5*vx(kbb,i)
  205   continue
        do 500 iloop=1,10
! ---------------
! central channel
! ---------------
          do 240 loop=1,20
            do 210 i=1,lcx
              rlm(l,i)=blm(l)
  210       continue
            psim=0.
            chim=0.
            fcm=3.*psi(l,1)/phir(l,1)-3.*psi(l,2)/phir(l,2)+psi(l,3)/phir(l,3)
            if (l.eq.3) gpsim=h2*h2m*2.*kf*fcm/rt5
            if (l.eq.4) gpsim=h2*rt2*rv(l)*fcm
            gchim=gpsim
            hm=0.
            psi0=psi(l,1)
            chi0=phir(l,1)
            chi(1)=chi0
            yc=vx(l,1)-rlm(l,1)
            yt=vx(m,1)-rlm(m,1)
            yb=vx(n,1)-rlm(n,1)
            yq=vx(kq,1)
            ybb=vx(kbb,1)
            pqpc=phir(n,1)/phir(l,1)
            pqpcs=pqpc*pqpc
            gpsi0=h2*(yc+pqpcs*(yq+(2./3.)*ybb)-pm(l,1))*psi0
            gchi0=h2*(yc+pqpcs*(yq+(2./3.)*ybb)-pm(l,1))*chi0
            h0=h2*((8.*yt-(2./3.)*pqpcs*ybb)*psi(m,1)+(2./3.)*pqpc*(yb-.5*ybb)*psi(n,1))
            do 220 j=2,lcxp
              yc=vx(l,j)-rlm(l,j)
              yt=vx(m,j)-rlm(m,j)
              yb=vx(n,j)-rlm(n,j)
              yq=vx(kq,j)
              ybb=vx(kbb,j)
              pqpc=phir(n,j)/phir(l,j)
              pqpcs=pqpc*pqpc
              gp=h2*(yc+pqpcs*(yq+(2./3.)*ybb)-pm(l,j))
              hp=h2*((8.*yt-(2./3.)*pqpcs*ybb)*psi(m,j)+(2./3.)*pqpc*(yb-.5*ybb)*psi(n,j))
              psi(l,j)=(2.*psi0+10.*gpsi0-psim+gpsim+hp+10.*h0+hm)/(1.-gp)
              psim=psi0
              gpsim=gpsi0
              hm=h0
              psi0=psi(l,j)
              gpsi0=gp*psi0
              h0=hp
              chi(j)=(2.*chi0+10.*gchi0-chim+gchim)/(1.-gp)
              chim=chi0
              gchim=gchi0
              chi0=chi(j)
              gchi0=gp*chi0
  220       continue
            fac=(phir(l,lcx)-psi(l,lcx))/chi(lcx)
            do 230 j=1,lcxp
              psi(l,j)=psi(l,j)+fac*chi(j)
  230       continue
            dldif=phir(l,lcx)*(psi(l,lcxp)-psi(l,lcxm))-psi(l,lcx)*(phir(l,lcxp)-phir(l,lcxm))
            if (loop.eq.1) then
              dldifo=dldif
              blmo=blm(l)
              blm(l)=blmo+(-1.)**(l-1)
            else
              if (abs(dldifo-dldif).le.small) go to 290
              blmn=(dldifo*blm(l)-dldif*blmo)/(dldifo-dldif)
              dldifo=dldif
              blmo=blm(l)
              blm(l)=blmn
            end if
  240     continue
  290     psi(l,lcxp)=phir(l,lcxp)
! --------------
! tensor channel
! --------------
          do 360 loop=1,20
            do 310 i=lcxp,ltxp
              xft=psi(m,i)/phir(m,i)
              xfb=psi(n,i)/phir(n,i)
              rlm(l,i)=rlx(l-2,i)+8.*(vx(m,i)-rlm(m,i))*xft &
     &         +(2./3.)*((vx(n,i)-rlm(n,i)-.5*vx(kbb,i))*xfb &
     &         -vx(kbb,i)*xft)*(phir(n,i)/phir(l,i))**2
  310       continue
            do 320 i=1,ltx
              rlm(m,i)=blm(m)
  320       continue
            do 330 i=llxp,ltxp
              xft=psi(m,i)/phir(m,i)
              rlm(n,i)=vx(n,i)-.5*vx(kbb,i)*(1.-4.*xft)/(1.-xft)
  330       continue
            psim=0.
            chim=0.
            gpsim=0.
            gchim=0.
            fcm=3*psi(l,1)/phir(l,1)-3*psi(l,2)/phir(l,2)+psi(l,3)/phir(l,3)
            if (l.eq.3) hm=0.
            if (l.eq.4) hm=h2*rv(m)*fcm*rt2
            psi0=psi(m,1)
            chi0=phir(m,1)
            chi(1)=chi0
            yc=vx(l,1)-rlm(l,1)
            yt=vx(m,1)-rlm(m,1)
            yb=vx(n,1)-rlm(n,1)
            yq=vx(kq,1)
            ybb=vx(kbb,1)
            pqpc=phir(n,1)/phir(l,1)
            pqpcs=pqpc*pqpc
            gpsi0=h2*(h2m*6./rsx(1)+yc-2.*yt-3.*yb+6.*yq+9.*ybb &
     &       +pqpcs*(yq+(5./6.)*ybb)-pm(m,1))*psi0
            gchi0=h2*(h2m*6./rsx(1)+yc-2.*yt-3.*yb+6.*yq+9.*ybb &
     &       +pqpcs*(yq+(5./6.)*ybb)-pm(m,1))*chi0
            h0=h2*((yt-(1./12.)*pqpcs*ybb)*psi(l,1) &
     &          -(1./12.)*pqpc*(yb-2.*ybb)*psi(n,1))
            do 340 j=2,ltxp
              yc=vx(l,j)-rlm(l,j)
              yt=vx(m,j)-rlm(m,j)
              yb=vx(n,j)-rlm(n,j)
              yq=vx(kq,j)
              ybb=vx(kbb,j)
              pqpc=phir(n,j)/phir(l,j)
              pqpcs=pqpc*pqpc
              gp=h2*(h2m*6./rsx(j)+yc-2.*yt-3.*yb+6.*yq+9.*ybb &
     &              +pqpcs*(yq+(5./6.)*ybb)-pm(m,j))
              hp=h2*((yt-(1./12.)*pqpcs*ybb)*psi(l,j) &
     &            -(1./12.)*pqpc*(yb-2.*ybb)*psi(n,j))
              psi(m,j)=(2.*psi0+10.*gpsi0-psim+gpsim+hp+10.*h0+hm)/(1.-gp)
              psim=psi0
              gpsim=gpsi0
              hm=h0
              psi0=psi(m,j)
              gpsi0=gp*psi0
              h0=hp
              chi(j)=(2*chi0+10*gchi0-chim+gchim)/(1-gp)
              chim=chi0
              gchim=gchi0
              chi0=chi(j)
              gchi0=gp*chi0
  340       continue
            fac=-psi(m,ltx)/chi(ltx)
            do 350 j=1,ltxp
              psi(m,j)=psi(m,j)+fac*chi(j)
  350       continue
            dldif=phir(m,ltx)*(psi(m,ltxp)-psi(m,ltxm)) &
     &           -psi(m,ltx)*(phir(m,ltxp)-phir(m,ltxm))
            if (loop.eq.1) then
              dldifo=dldif
              blmo=blm(m)
              blm(m)=blmo+(-1.)**(l-1)
            else
              if (abs(dldifo-dldif).le.small) go to 390
              blmn=(dldifo*blm(m)-dldif*blmo)/(dldifo-dldif)
              dldifo=dldif
              blmo=blm(m)
              blm(m)=blmn
            end if
  360     continue
  390     psi(m,ltxp)=0.
! ------------------
! spin-orbit channel
! ------------------
          do 460 loop=1,20
            do 410 i=lcxp,ltxp
              xft=psi(m,i)/phir(m,i)
              xfb=psi(n,i)/phir(n,i)
              rlm(l,i)=rlx(l-2,i)+8*(vx(m,i)-rlm(m,i))*xft &
     &         +(2./3.)*((vx(n,i)-rlm(n,i)-.5*vx(kbb,i))*xfb &
     &         -vx(kbb,i)*xft)*(phir(n,i)/phir(l,i))**2
  410       continue
            do 420 i=1,llx
              rlm(n,i)=blm(n)
  420       continue
            psim=0
            chim=0
            fbm=3*psi(n,1)/phir(n,1)-3*psi(n,2)/phir(n,2)+psi(n,3)/phir(n,3)
            if (l.eq.3) gpsim=h2*h2m*2*kf*fbm*rt2/rt5
            if (l.eq.4) gpsim=0
            gchim=gpsim
            hm=0
            psi0=psi(n,1)
            chi0=phir(n,1)
            chi(1)=chi0
            yc=vx(l,1)-rlm(l,1)
            yt=vx(m,1)-rlm(m,1)
            yb=vx(n,1)-rlm(n,1)
            yq=vx(kq,1)
            ybb=vx(kbb,1)
            pqpc=phir(n,1)/phir(l,1)
            pqpcs=pqpc*pqpc
            pqqpqs=(phir(kqq,1)/phir(n,1))**2
            gpsi0=h2*(yc-yt-.5*yb+pqqpqs*(yq+ybb)-pm(n,1))*psi0
            gchi0=h2*(yc-yt-.5*yb+pqqpqs*(yq+ybb)-pm(n,1))*chi0
            h0=h2*pqpc*((yb-.5*ybb)*psi(l,1)-(yb-2*ybb)*psi(m,1))
            do 440 j=2,llxp
              yc=vx(l,j)-rlm(l,j)
              yt=vx(m,j)-rlm(m,j)
              yb=vx(n,j)-rlm(n,j)
              yq=vx(kq,j)
              ybb=vx(kbb,j)
              pqpc=phir(n,j)/phir(l,j)
              pqpcs=pqpc*pqpc
              pqqpqs=(phir(kqq,j)/phir(n,j))**2
              gp=h2*(yc-yt-.5*yb+pqqpqs*(yq+ybb)-pm(n,j))
              hp=h2*pqpc*((yb-.5*ybb)*psi(l,j)-(yb-2*ybb)*psi(m,j))
              psi(n,j)=(2*psi0+10*gpsi0-psim+gpsim+hp+10*h0+hm)/(1-gp)
              psim=psi0
              gpsim=gpsi0
              hm=h0
              psi0=psi(n,j)
              gpsi0=gp*psi0
              h0=hp
              chi(j)=(2*chi0+10*gchi0-chim+gchim)/(1-gp)
              chim=chi0
              gchim=gchi0
              chi0=chi(j)
              gchi0=gp*chi0
  440       continue
            fac=-psi(n,llx)/chi(llx)
            do 450 j=1,llxp
              psi(n,j)=psi(n,j)+fac*chi(j)
  450       continue
            dldif=phir(n,llx)*(psi(n,llxp)-psi(n,llxm)) &
     &           -psi(n,llx)*(phir(n,llxp)-phir(n,llxm))
            if (loop.eq.1) then
              dldifo=dldif
              blmo=blm(n)
              blm(n)=blmo+(-1)**k
            else
              if (abs(dldifo-dldif).le.small) go to 490
              blmn=(dldifo*blm(n)-dldif*blmo)/(dldifo-dldif)
              dldifo=dldif
              blmo=blm(n)
              blm(n)=blmn
            end if
  460     continue
  490     psi(n,llxp)=0.
  500   continue
  510 continue


! print ================================================================
!     if (myrank().eq.0) write(nout,950) ltx
  950 format(/4x,'wave equations solved in t,s projected channels' &
     & /4x,'fine grid has ',i5,' points in dt')
!     if (myrank().eq.0) write(nout,951) blm
  951 format (/4x,'lambda'/4x,'1s',6x,'1p',4x,'3p(c)',3x,'3s(c)',3x &
     & ,'3p(t)',3x,'3s(t)',3x,'3p(b)',3x,'3s(b)'/8f8.3)
! ======================================================================
! -------------------
! projections for r<d
! -------------------

      lfit=lcx+1
      ltfit=ltx
      llfit=llx+1
      do i=1,ltfit
       j=i
       jp=j+1
       jm=j-1
       do k=1,8,nm
        blm(k)=rlm(k,j)
        cp=psi(k,jp)/phir(k,jp)
        c0=psi(k,j)/phir(k,j)
        cm=0.!psi(k,jm)/phir(k,jm)
        c(1,k)=c0
        c(2,k)=.5*(cp-cm)/h
        c(3,k)=(cp-2*c0+cm)/h**2+2*c(2,k)/rx(j)
       enddo



! I canali dovrebbero essere #=(S,T) :
!              1=(0,1) 2=(0,0) 3=(1,1) 4=(1,0)
! tensore:     5=(1,1) 6=(1,0)
! spin-orbita: 7=(?,?) 8=(?,?)
         
! Reprojection on operator basis
       do k=1,3
        ca(k,1)=.25*(3*c(k,3)+c(k,1))
        ca(k,3)=.25*(c(k,3)-c(k,1))
        ca(k,5)=c(k,5)
        ca(k,7)=c(k,7)
        if (nm.eq.1) then 
         ca(k,2)=.25*ca(k,1)-.0625*(3*c(k,4)+c(k,2))
         ca(k,4)=.25*ca(k,3)-.0625*(c(k,4)-c(k,2))
         ca(k,6)=.25*(ca(k,5)-c(k,6))
         ca(k,8)=.25*(ca(k,7)-c(k,8))
         ca(k,1)=.75*ca(k,1)+.0625*(3*c(k,4)+c(k,2))
         ca(k,3)=.75*ca(k,3)+.0625*(c(k,4)-c(k,2))
         ca(k,5)=.25*(3*ca(k,5)+c(k,6))
         ca(k,7)=.25*(3*ca(k,7)+c(k,8))
        endif
       enddo

       do k=1,8,nm
        f(i,k)=ca(1,k)
        fp(i,k)=ca(2,k)
        if (k.le.4) then
         fpp(i,k)=ca(3,k)
        else if (k.ge.7) then
         fpp(i,k)=ca(3,k)+2.*ca(2,k)/rx(i)
        else
         fpp(i,k)=ca(3,k)-6.*ca(1,k)/rsx(i)
        end if
       enddo
       
       if(i.gt.lfit)then
        do k=1,4,nm
         f(i,k)=0.
         fp(i,k)=0.
         fpp(i,k)=0.
         f(i,1)=1.
        enddo
       end if

       if (i.gt.llfit) then
        do k=7,8,nm
         f(i,k)=0.
         fp(i,k)=0.
         fpp(i,k)=0.
        enddo
       end if

!     if (myrank().eq.0) then
!     if(nm.eq.2) then
!      write(1021,'(5(E13.6,1x))') rx(i),f(i,1),f(i,3),f(i,5),f(i,7)
!     elseif(nm.eq.1)then
!      write(1021,'(7(E13.6,1x))') rx(i),f(i,1),f(i,2),f(i,3),f(i,4),f(i,5),f(i,6)
!     endif
!     endif

     enddo

!      fpp(1,:)=2.*fpp(2,:)-fpp(3,:)
!      fp(1,:)=2.*fp(2,:)-fp(3,:)
!      f(1,:)=2.*f(2,:)-f(3,:)

!
! end tabulation of the correlations
!

      return
      end subroutine

      subroutine vpot_op(lpot,lr,r,vv,ww)
      use cheft
      real(kind=r8) :: r,rr,vv(18),ww(14),vp(12)
      integer(kind=i4) :: lr
      vv=0.0d0
      ww=0.0d0
      rr=r
      if (lpot.lt.26) then
         call pot(0,rr,vv,vp,ww)
      elseif (lpot.eq.26) then
         v0r=200.0_r8  ! MeV
         v0t=178.0_r8  ! MeV
         v0s=91.85_r8  ! MeV
         kr=1.487_r8  ! fm**-2
         kt=0.639_r8  ! fm**-2
         ks=0.465_r8  ! fm**-2
         vr=v0r*exp(-kr*r**2)
         vt=-v0t*exp(-kt*r**2)
         vs=-v0s*exp(-ks*r**2)
         vv(1)=3.0_r8/8.0_r8*(vr+0.5_r8*vt+0.5_r8*vs)
         vv(2)=1.0_r8/8.0_r8*(-vr-1.5_r8*vt+0.5_r8*vs)
         vv(3)=1.0_r8/8.0_r8*(-vr+0.5_r8*vt-1.5_r8*vs)
         vv(4)=1.0_r8/8.0_r8*(-vr-0.5_r8*vt-0.5_r8*vs)
      elseif (lpot.gt.100) then
         call cheft_pot(lpot,rr,vv)
      endif
      vv(:)=vv(:)*r**lr
      ww(:)=ww(:)*r**lr
      return
      end  subroutine
end module
