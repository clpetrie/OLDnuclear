module orbital
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   integer(kind=i4), private,save :: nk,nex
   real(kind=r8), private, save, pointer :: ak(:,:),alphak(:),ak2(:)
   real(kind=r8), private, save, allocatable :: akex(:,:),ak2ex(:)
   complex(kind=r8), private, save, allocatable :: spex(:,:),sp(:,:)

contains
   subroutine setorbital(nsh,vsh,el,nex,ikex,spexin,spin)
   use kshellold
   use correlation
   integer(kind=i4) :: i
   integer(kind=i4) :: nsh,nex,ikex(3,nex)
   complex(kind=r8) :: spin(4,4),spexin(4,nex)
   real(kind=r8) :: el,pi,vsh(:)
   pi=4.0_r8*atan(1.0_r8)
   nk=0
   if (nsh.ne.0) call setupk(el,nsh,vsh,ak,ak2,alphak,nk)
   allocate(akex(3,nex),ak2ex(nex),spex(4,nex))
   akex=2.0_r8*pi*ikex/el
!  do i=1,nex
!     ak2ex(i)=sum(akex(:,i)**2)
!  enddo
   spex=spexin
   allocate(sp(4,4))
   sp=spin
!  call readcorr
   end subroutine setorbital

!  subroutine pairfn(x,pf,dpf,d2pf,el)
   subroutine pairfn(x,pf,el)
   use correlation
   real(kind=r8) :: x(3)
   complex(kind=r8) :: pf(4,4),dpf(4,4,3),d2pf(4,4),ex,phip(4,4),phim(4,4)
   real(kind=r8) :: arg,r
   integer(kind=i4) :: i,j
   real(kind=r8) :: f(6),df(6),d2f(6),el
   pf=czero
!  dpf=czero
!  d2pf=czero
   do i=1,nk
      arg=ak(1,i)*x(1)+ak(2,i)*x(2)+ak(3,i)*x(3)
      ex=cmplx(cos(arg),sin(arg),r8)
      phip=alphak(i)*ex*sp(:,:)
      phim=-alphak(i)*conjg(ex)*transpose(sp(:,:))
      pf=pf+phip+phim
!     d2pf=d2pf-ak2(i)*(phip+phim)
!     dpf(:,:,1)=dpf(:,:,1)+ci*ak(1,i)*(phip-phim)
!     dpf(:,:,2)=dpf(:,:,2)+ci*ak(2,i)*(phip-phim)
!     dpf(:,:,3)=dpf(:,:,3)+ci*ak(3,i)*(phip-phim)
   enddo
!  r=sqrt(sum(x**2))
!  if (r.lt.el*0.5_r8) then
!     call corr(r,f,df,d2f)
!     call correlatepair(pf,dpf,d2pf,f,df,d2f,x)
!  else
!     f=0.0_r8
!     df=0.0_r8
!     d2f=0.0_r8
!  endif
   end subroutine pairfn

!  subroutine singlefn(ifun,x,ph,dph,d2ph)
   subroutine singlefn(ifun,x,ph)
   integer(kind=i4) :: ifun
   real(kind=r8) :: x(3)
   complex(kind=r8) :: ph(4),dph(4,3),d2ph(4)
   complex(kind=r8) :: p,dp(3),d2p
   real(kind=r8) :: arg
   arg=akex(1,ifun)*x(1)+akex(2,ifun)*x(2)+akex(3,ifun)*x(3)
   p=cmplx(cos(arg),-sin(arg))
!  dp=ci*akex(:,ifun)*p
!  d2p=-ak2ex(ifun)*p
   ph=spex(:,ifun)*p
!  dph(:,1)=spex(:,ifun)*dp(1)
!  dph(:,2)=spex(:,ifun)*dp(2)
!  dph(:,3)=spex(:,ifun)*dp(3)
!  d2ph=spex(:,ifun)*d2p
   return
   end subroutine singlefn

   subroutine correlatepair(pf,dpf,d2pf,f,fp,fpp,dx)
   complex(kind=r8) :: pf(4,4),dpf(4,4,3),d2pf(4,4)
   complex(kind=r8) :: pfout(4,4),dpfout(4,4,3),d2pfout(4,4),xp
   complex(kind=r8) :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15
   complex(kind=r8) :: c16,c17,c18,c19,c20,c30,c31,c32,c33,c34,c35,c36,c37
   complex(kind=r8) :: c21,c22,c23,c24,c25,c26,c27,c28,c29
   complex(kind=r8) :: c38,c39,c40,c41,c42,dc1(3),dc2(3),dc3(3),dc4(3)
   complex(kind=r8) :: dc5(3),dc6(3),dc7(3),dc8(3),dc9(3),dc10(3),dc11(3)
   complex(kind=r8) :: dc12(3),dc13(3),dc14(3),dc15(3),dc16(3),dc17(3)
   complex(kind=r8) :: dc18(3),dc19(3),dc20(3),dc21(3),dc22(3)
   complex(kind=r8) :: dc23(3),dc24(3),dc25(3),dc26(3),dc27(3),dc28(3)
   complex(kind=r8) :: dc29(3),dc30(3),dc31(3),dc32(3),dc33(3),dc34(3)
   complex(kind=r8) :: dc35(3),dc36(3),dc37(3),dc38(3),dc39(3),dc40(3)
   complex(kind=r8) :: dc41(3),dc42(3),d2c1,d2c2,d2c3,d2c4,d2c5,d2c6
   complex(kind=r8) :: d2c7,d2c8,d2c9,d2c10,d2c11,d2c12,d2c13,d2c14
   complex(kind=r8) :: d2c15,d2c16,d2c17,d2c18,d2c19,d2c20,d2c21,d2c22
   complex(kind=r8) :: d2c23,d2c24,d2c25,d2c26,d2c27,d2c28,d2c29,d2c30
   complex(kind=r8) :: d2c31,d2c32,d2c33,d2c34,d2c35,d2c36,d2c37,d2c38
   complex(kind=r8) :: d2c39,d2c40,d2c41,d2c42,ft(-2:2),ftt(-2:2),dft(3,-2:2)
   complex(kind=r8) :: dftt(3,-2:2),d2ft(-2:2),d2ftt(-2:2),y(-2:2),dy(3,-2:2)
   real(kind=r8) :: f(6),fp(6),fpp(6),dx(3),xhat(3),r,ri,df(3,4),d2f(4)
   integer(kind=i4) :: ipu,ipd,inu,ind,i
   r=sqrt(sum(dx**2))
   ri=1.0_r8/r
   xhat=dx*ri
   xp=cmplx(xhat(1),-xhat(2))
   do i=1,4
      df(:,i)=xhat(:)*fp(i)
      d2f(i)=fpp(i)+2.0*ri*fp(i)
   enddo
   y(0)=3.0_r8*xhat(3)**2-1.0_r8
   y(1)=3.0_r8*xhat(3)*xp
   y(2)=3.0_r8*xp**2
   do i=0,2
      dy(:,i)=-xhat(:)*2.0_r8*y(i)*ri
   enddo
   dy(1,0)=dy(1,0)-2.0_r8*xhat(1)*ri
   dy(2,0)=dy(2,0)-2.0_r8*xhat(2)*ri
   dy(3,0)=dy(3,0)+4.0_r8*xhat(3)*ri
   dy(1,1)=dy(1,1)+3.0_r8*xhat(3)*ri
   dy(2,1)=dy(2,1)-3.0_r8*ci*xhat(3)*ri
   dy(3,1)=dy(3,1)+3.0_r8*xp*ri
   dy(1,2)=dy(1,2)+6.0_r8*xp*ri
   dy(2,2)=dy(2,2)-6.0_r8*ci*xp*ri
   do i=1,2
      y(-i)=conjg(y(i))
      dy(:,-i)=conjg(dy(:,i))
   enddo
   ft=f(5)*y
   ftt=f(6)*y
   d2ft=(fpp(5)+2.0_r8*ri*fp(5)-6.0_r8*ri**2*f(5))*y
   d2ftt=(fpp(6)+2.0_r8*ri*fp(6)-6.0_r8*ri**2*f(6))*y
   do i=-2,2
      dft(:,i)=xhat(:)*fp(5)*y(i)+f(5)*dy(:,i)
      dftt(:,i)=xhat(:)*fp(6)*y(i)+f(6)*dy(:,i)
   enddo

   c1=f(1)+f(2)+f(3)+f(4)+ft(0)+ftt(0)
   dc1=df(:,1)+df(:,2)+df(:,3)+df(:,4)+dft(:,0)+dftt(:,0)
   d2c1=d2f(1)+d2f(2)+d2f(3)+d2f(4)+d2ft(0)+d2ftt(0)
   c2=ft(-1)+ftt(-1)
   dc2=dft(:,-1)+dftt(:,-1)
   d2c2=d2ft(-1)+d2ftt(-1)
   c3=ft(-2)+ftt(-2)
   dc3=dft(:,-2)+dftt(:,-2)
   d2c3=d2ft(-2)+d2ftt(-2)

   c4=f(1)+f(2)-f(3)-f(4)-ft(0)-ftt(0)
   dc4=df(:,1)+df(:,2)-df(:,3)-df(:,4)-dft(:,0)-dftt(:,0)
   d2c4=d2f(1)+d2f(2)-d2f(3)-d2f(4)-d2ft(0)-d2ftt(0)
   c5=2.0_r8*(f(3)+f(4))-ft(0)-ftt(0)
   dc5=2.0_r8*(df(:,3)+df(:,4))-dft(:,0)-dftt(:,0)
   d2c5=2.0_r8*(d2f(3)+d2f(4))-d2ft(0)-d2ftt(0)
   c6=ft(1)+ftt(1)
   dc6=dft(:,1)+dftt(:,1)
   d2c6=d2ft(1)+d2ftt(1)
   c7=-ft(-1)-ftt(-1)
   dc7=-dft(:,-1)-dftt(:,-1)
   d2c7=-d2ft(-1)-d2ftt(-1)

   c8=f(1)-f(2)+f(3)-f(4)+ft(0)-ftt(0)
   dc8=df(:,1)-df(:,2)+df(:,3)-df(:,4)+dft(:,0)-dftt(:,0)
   d2c8=d2f(1)-d2f(2)+d2f(3)-d2f(4)+d2ft(0)-d2ftt(0)
   c9=ft(-1)-ftt(-1)
   dc9=dft(:,-1)-dftt(:,-1)
   d2c9=d2ft(-1)-d2ftt(-1)
   c10=ft(-2)-ftt(-2)
   dc10=dft(:,-2)-dftt(:,-2)
   d2c10=d2ft(-2)-d2ftt(-2)
   c11=2.0_r8*(f(2)+f(4)+ftt(0))
   dc11=2.0_r8*(df(:,2)+df(:,4)+dftt(:,0))
   d2c11=2.0_r8*(d2f(2)+d2f(4)+d2ftt(0))
   c12=2.0_r8*ftt(-1)
   dc12=2.0_r8*dftt(:,-1)
   d2c12=2.0_r8*d2ftt(-1)
   c13=2.0_r8*ftt(-2)
   dc13=2.0_r8*dftt(:,-2)
   d2c13=2.0_r8*d2ftt(-2)

   c14=f(1)-f(2)-f(3)+f(4)-ft(0)+ftt(0)
   dc14=df(:,1)-df(:,2)-df(:,3)+df(:,4)-dft(:,0)+dftt(:,0)
   d2c14=d2f(1)-d2f(2)-d2f(3)+d2f(4)-d2ft(0)+d2ftt(0)
   c15=2.0_r8*(f(3)-f(4))-ft(0)+ftt(0)
   dc15=2.0_r8*(df(:,3)-df(:,4))-dft(:,0)+dftt(:,0)
   d2c15=2.0_r8*(d2f(3)-d2f(4))-d2ft(0)+d2ftt(0)
   c16=ft(1)-ftt(1)
   dc16=dft(:,1)-dftt(:,1)
   d2c16=d2ft(1)-d2ftt(1)
   c17=-ft(-1)+ftt(-1)
   dc17=-dft(:,-1)+dftt(:,-1)
   d2c17=-d2ft(-1)+d2ftt(-1)
   c18=2.0_r8*(f(2)-f(4)-ftt(0))
   dc18=2.0_r8*(df(:,2)-df(:,4)-dftt(:,0))
   d2c18=2.0_r8*(d2f(2)-d2f(4)-d2ftt(0))
   c19=4.0_r8*f(4)-2.0_r8*ftt(0)
   dc19=4.0_r8*df(:,4)-2.0_r8*dftt(:,0)
   d2c19=4.0_r8*d2f(4)-2.0_r8*d2ftt(0)
   c20=2.0_r8*ftt(1)
   dc20=2.0_r8*dftt(:,1)
   d2c20=2.0_r8*d2ftt(1)
   c21=-2.0_r8*ftt(-1)
   dc21=-2.0_r8*dftt(:,-1)
   d2c21=-2.0_r8*d2ftt(-1)

   c22=dconjg(c1)
   dc22=dconjg(dc1)
   d2c22=dconjg(d2c1)
   c23=-dconjg(c2)
   dc23=-dconjg(dc2)
   d2c23=-dconjg(d2c2)
   c24=dconjg(c3)
   dc24=dconjg(dc3)
   d2c24=dconjg(d2c3)
   c25=dconjg(c4)
   dc25=dconjg(dc4)
   d2c25=dconjg(d2c4)
   c26=dconjg(c5)
   dc26=dconjg(dc5)
   d2c26=dconjg(d2c5)
   c27=-dconjg(c6)
   dc27=-dconjg(dc6)
   d2c27=-dconjg(d2c6)
   c28=-dconjg(c7)
   dc28=-dconjg(dc7)
   d2c28=-dconjg(d2c7)
   c29=dconjg(c8)
   dc29=dconjg(dc8)
   d2c29=dconjg(d2c8)
   c30=-dconjg(c9)
   dc30=-dconjg(dc9)
   d2c30=-dconjg(d2c9)
   c31=dconjg(c10)
   dc31=dconjg(dc10)
   d2c31=dconjg(d2c10)
   c32=dconjg(c11)
   dc32=dconjg(dc11)
   d2c32=dconjg(d2c11)
   c33=-dconjg(c12)
   dc33=-dconjg(dc12)
   d2c33=-dconjg(d2c12)
   c34=dconjg(c13)
   dc34=dconjg(dc13)
   d2c34=dconjg(d2c13)
   c35=dconjg(c14)
   dc35=dconjg(dc14)
   d2c35=dconjg(d2c14)
   c36=dconjg(c15)
   dc36=dconjg(dc15)
   d2c36=dconjg(d2c15)
   c37=-dconjg(c16)
   dc37=-dconjg(dc16)
   d2c37=-dconjg(d2c16)
   c38=-dconjg(c17)
   dc38=-dconjg(dc17)
   d2c38=-dconjg(d2c17)
   c39=dconjg(c18)
   dc39=dconjg(dc18)
   d2c39=dconjg(d2c18)
   c40=dconjg(c19)
   dc40=dconjg(dc19)
   d2c40=dconjg(d2c19)
   c41=-dconjg(c20)
   dc41=-dconjg(dc20)
   d2c41=-dconjg(d2c20)
   c42=-dconjg(c21)
   dc42=-dconjg(dc21)
   d2c42=-dconjg(d2c21)
   
   do i=1,2
      select case (i)
      case (1) ! as written below
         ipu=1
         ipd=2
         inu=3
         ind=4
      case (2) ! interchange proton and neutron
         ipu=3
         ipd=4
         inu=1
         ind=2
      end select
      pfout(ipu,ipu)=c1*pf(ipu,ipu)+c2*(pf(ipu,ipd)+pf(ipd,ipu))+c3*pf(ipd,ipd)
      dpfout(ipu,ipu,:)=&
          c1*dpf(ipu,ipu,:)+dc1*pf(ipu,ipu)&
         +c2*(dpf(ipu,ipd,:)+dpf(ipd,ipu,:))+dc2*(pf(ipu,ipd)+pf(ipd,ipu))&
         +c3*dpf(ipd,ipd,:)+dc3*pf(ipd,ipd)
      d2pfout(ipu,ipu)=&
          c1*d2pf(ipu,ipu)+d2c1*pf(ipu,ipu)&
         +c2*(d2pf(ipu,ipd)+d2pf(ipd,ipu))+d2c2*(pf(ipu,ipd)+pf(ipd,ipu))&
         +c3*d2pf(ipd,ipd)+d2c3*pf(ipd,ipd)&
         +2.0_r8*sum(&
          dc1*dpf(ipu,ipu,:)&
         +dc2*(dpf(ipu,ipd,:)+dpf(ipd,ipu,:))&
         +dc3*dpf(ipd,ipd,:))
      pfout(ipu,ipd)=c4*pf(ipu,ipd)+c5*pf(ipd,ipu)+c6*pf(ipu,ipu)+c7*pf(ipd,ipd)
      dpfout(ipu,ipd,:)=&
          c4*dpf(ipu,ipd,:)+dc4*pf(ipu,ipd)&
         +c5*dpf(ipd,ipu,:)+dc5*pf(ipd,ipu)&
         +c6*dpf(ipu,ipu,:)+dc6*pf(ipu,ipu)&
         +c7*dpf(ipd,ipd,:)+dc7*pf(ipd,ipd)
      d2pfout(ipu,ipd)=&
          c4*d2pf(ipu,ipd)+d2c4*pf(ipu,ipd)&
         +c5*d2pf(ipd,ipu)+d2c5*pf(ipd,ipu)&
         +c6*d2pf(ipu,ipu)+d2c6*pf(ipu,ipu)&
         +c7*d2pf(ipd,ipd)+d2c7*pf(ipd,ipd)&
         +2.0_r8*sum(&
          dc4*dpf(ipu,ipd,:)&
         +dc5*dpf(ipd,ipu,:)&
         +dc6*dpf(ipu,ipu,:)&
         +dc7*dpf(ipd,ipd,:))
      pfout(ipu,inu)=c8*pf(ipu,inu)+c9*(pf(ipu,ind)+pf(ipd,inu))+c10*pf(ipd,ind)&
         +c11*pf(inu,ipu)+c12*(pf(inu,ipd)+pf(ind,ipu))+c13*pf(ind,ipd)
      dpfout(ipu,inu,:)=&
          c8*dpf(ipu,inu,:)+dc8*pf(ipu,inu)&
         +c9*(dpf(ipu,ind,:)+dpf(ipd,inu,:))+dc9*(pf(ipu,ind)+pf(ipd,inu))&
         +c10*dpf(ipd,ind,:)+dc10*pf(ipd,ind)&
         +c11*dpf(inu,ipu,:)+dc11*pf(inu,ipu)&
         +c12*(dpf(inu,ipd,:)+dpf(ind,ipu,:))+dc12*(pf(inu,ipd)+pf(ind,ipu))&
         +c13*dpf(ind,ipd,:)+dc13*pf(ind,ipd)
      d2pfout(ipu,inu)=&
          c8*d2pf(ipu,inu)+d2c8*pf(ipu,inu)&
         +c9*(d2pf(ipu,ind)+d2pf(ipd,inu))+d2c9*(pf(ipu,ind)+pf(ipd,inu))&
         +c10*d2pf(ipd,ind)+d2c10*pf(ipd,ind)&
         +c11*d2pf(inu,ipu)+d2c11*pf(inu,ipu)&
         +c12*(d2pf(inu,ipd)+d2pf(ind,ipu))+d2c12*(pf(inu,ipd)+pf(ind,ipu))&
         +c13*d2pf(ind,ipd)+d2c13*pf(ind,ipd)&
         +2.0_r8*sum(&
          dc8*dpf(ipu,inu,:)&
         +dc9*(dpf(ipu,ind,:)+dpf(ipd,inu,:))&
         +dc10*dpf(ipd,ind,:)&
         +dc11*dpf(inu,ipu,:)&
         +dc12*(dpf(inu,ipd,:)+dpf(ind,ipu,:))&
         +dc13*dpf(ind,ipd,:))
      pfout(ipu,ind)=c14*pf(ipu,ind)+c15*pf(ipd,inu)+c16*pf(ipu,inu)+c17*pf(ipd,ind)&
         +c18*pf(inu,ipd)+c19*pf(ind,ipu)+c20*pf(inu,ipu)+c21*pf(ind,ipd)
      dpfout(ipu,ind,:)=&
          c14*dpf(ipu,ind,:)+dc14*pf(ipu,ind)&
         +c15*dpf(ipd,inu,:)+dc15*pf(ipd,inu)&
         +c16*dpf(ipu,inu,:)+dc16*pf(ipu,inu)&
         +c17*dpf(ipd,ind,:)+dc17*pf(ipd,ind)&
         +c18*dpf(inu,ipd,:)+dc18*pf(inu,ipd)&
         +c19*dpf(ind,ipu,:)+dc19*pf(ind,ipu)&
         +c20*dpf(inu,ipu,:)+dc20*pf(inu,ipu)&
         +c21*dpf(ind,ipd,:)+dc21*pf(ind,ipd)
      d2pfout(ipu,ind)=&
          c14*d2pf(ipu,ind)+d2c14*pf(ipu,ind)&
         +c15*d2pf(ipd,inu)+d2c15*pf(ipd,inu)&
         +c16*d2pf(ipu,inu)+d2c16*pf(ipu,inu)&
         +c17*d2pf(ipd,ind)+d2c17*pf(ipd,ind)&
         +c18*d2pf(inu,ipd)+d2c18*pf(inu,ipd)&
         +c19*d2pf(ind,ipu)+d2c19*pf(ind,ipu)&
         +c20*d2pf(inu,ipu)+d2c20*pf(inu,ipu)&
         +c21*d2pf(ind,ipd)+d2c21*pf(ind,ipd)&
         +2.0_r8*sum(&
          dc14*dpf(ipu,ind,:)&
         +dc15*dpf(ipd,inu,:)&
         +dc16*dpf(ipu,inu,:)&
         +dc17*dpf(ipd,ind,:)&
         +dc18*dpf(inu,ipd,:)&
         +dc19*dpf(ind,ipu,:)&
         +dc20*dpf(inu,ipu,:)&
         +dc21*dpf(ind,ipd,:))
      pfout(ipd,ipd)=c22*pf(ipd,ipd)+c23*(pf(ipd,ipu)+pf(ipu,ipd))+c24*pf(ipu,ipu)
      dpfout(ipd,ipd,:)=&
          c22*dpf(ipd,ipd,:)+dc22*pf(ipd,ipd)&
         +c23*(dpf(ipd,ipu,:)+dpf(ipu,ipd,:))+dc23*(pf(ipd,ipu)+pf(ipu,ipd))&
         +c24*dpf(ipu,ipu,:)+dc24*pf(ipu,ipu)
      d2pfout(ipd,ipd)=&
          c22*d2pf(ipd,ipd)+d2c22*pf(ipd,ipd)&
         +c23*(d2pf(ipd,ipu)+d2pf(ipu,ipd))+d2c23*(pf(ipd,ipu)+pf(ipu,ipd))&
         +c24*d2pf(ipu,ipu)+d2c24*pf(ipu,ipu)&
         +2.0_r8*sum(&
          dc22*dpf(ipd,ipd,:)&
         +dc23*(dpf(ipd,ipu,:)+dpf(ipu,ipd,:))&
         +dc24*dpf(ipu,ipu,:))
      pfout(ipd,ipu)=c25*pf(ipd,ipu)+c26*pf(ipu,ipd)+c27*pf(ipd,ipd)+c28*pf(ipu,ipu)
      dpfout(ipd,ipu,:)=&
          c25*dpf(ipd,ipu,:)+dc25*pf(ipd,ipu)&
         +c26*dpf(ipu,ipd,:)+dc26*pf(ipu,ipd)&
         +c27*dpf(ipd,ipd,:)+dc27*pf(ipd,ipd)&
         +c28*dpf(ipu,ipu,:)+dc28*pf(ipu,ipu)
      d2pfout(ipd,ipu)=&
          c25*d2pf(ipd,ipu)+d2c25*pf(ipd,ipu)&
         +c26*d2pf(ipu,ipd)+d2c26*pf(ipu,ipd)&
         +c27*d2pf(ipd,ipd)+d2c27*pf(ipd,ipd)&
         +c28*d2pf(ipu,ipu)+d2c28*pf(ipu,ipu)&
         +2.0_r8*sum(&
          dc25*dpf(ipd,ipu,:)&
         +dc26*dpf(ipu,ipd,:)&
         +dc27*dpf(ipd,ipd,:)&
         +dc28*dpf(ipu,ipu,:))
      pfout(ipd,ind)=c29*pf(ipd,ind)+c30*(pf(ipd,inu)+pf(ipu,ind))+c31*pf(ipu,inu)&
         +c32*pf(ind,ipd)+c33*(pf(ind,ipu)+pf(inu,ipd))+c34*pf(inu,ipu)
      dpfout(ipd,ind,:)=&
          c29*dpf(ipd,ind,:)+dc29*pf(ipd,ind)&
         +c30*(dpf(ipd,inu,:)+dpf(ipu,ind,:))+dc30*(pf(ipd,inu)+pf(ipu,ind))&
         +c31*dpf(ipu,inu,:)+dc31*pf(ipu,inu)&
         +c32*dpf(ind,ipd,:)+dc32*pf(ind,ipd)&
         +c33*(dpf(ind,ipu,:)+dpf(inu,ipd,:))+dc33*(pf(ind,ipu)+pf(inu,ipd))&
         +c34*dpf(inu,ipu,:)+dc34*pf(inu,ipu)
      d2pfout(ipd,ind)=&
          c29*d2pf(ipd,ind)+d2c29*pf(ipd,ind)&
         +c30*(d2pf(ipd,inu)+d2pf(ipu,ind))+d2c30*(pf(ipd,inu)+pf(ipu,ind))&
         +c31*d2pf(ipu,inu)+d2c31*pf(ipu,inu)&
         +c32*d2pf(ind,ipd)+d2c32*pf(ind,ipd)&
         +c33*(d2pf(ind,ipu)+d2pf(inu,ipd))+d2c33*(pf(ind,ipu)+pf(inu,ipd))&
         +c34*d2pf(inu,ipu)+d2c34*pf(inu,ipu)&
         +2.0_r8*sum(&
          dc29*dpf(ipd,ind,:)&
         +dc30*(dpf(ipd,inu,:)+dpf(ipu,ind,:))&
         +dc31*dpf(ipu,inu,:)&
         +dc32*dpf(ind,ipd,:)&
         +dc33*(dpf(ind,ipu,:)+dpf(inu,ipd,:))&
         +dc34*dpf(inu,ipu,:))
      pfout(ipd,inu)=c35*pf(ipd,inu)+c36*pf(ipu,ind)+c37*pf(ipd,ind)+c38*pf(ipu,inu)&
         +c39*pf(ind,ipu)+c40*pf(inu,ipd)+c41*pf(ind,ipd)+c42*pf(inu,ipu)
      dpfout(ipd,inu,:)=&
          c35*dpf(ipd,inu,:)+dc35*pf(ipd,inu)&
         +c36*dpf(ipu,ind,:)+dc36*pf(ipu,ind)&
         +c37*dpf(ipd,ind,:)+dc37*pf(ipd,ind)&
         +c38*dpf(ipu,inu,:)+dc38*pf(ipu,inu)&
         +c39*dpf(ind,ipu,:)+dc39*pf(ind,ipu)&
         +c40*dpf(inu,ipd,:)+dc40*pf(inu,ipd)&
         +c41*dpf(ind,ipd,:)+dc41*pf(ind,ipd)&
         +c42*dpf(inu,ipu,:)+dc42*pf(inu,ipu)
      d2pfout(ipd,inu)=&
          c35*d2pf(ipd,inu)+d2c35*pf(ipd,inu)&
         +c36*d2pf(ipu,ind)+d2c36*pf(ipu,ind)&
         +c37*d2pf(ipd,ind)+d2c37*pf(ipd,ind)&
         +c38*d2pf(ipu,inu)+d2c38*pf(ipu,inu)&
         +c39*d2pf(ind,ipu)+d2c39*pf(ind,ipu)&
         +c40*d2pf(inu,ipd)+d2c40*pf(inu,ipd)&
         +c41*d2pf(ind,ipd)+d2c41*pf(ind,ipd)&
         +c42*d2pf(inu,ipu)+d2c42*pf(inu,ipu)&
         +2.0_r8*sum(&
          dc35*dpf(ipd,inu,:)&
         +dc36*dpf(ipu,ind,:)&
         +dc37*dpf(ipd,ind,:)&
         +dc38*dpf(ipu,inu,:)&
         +dc39*dpf(ind,ipu,:)&
         +dc40*dpf(inu,ipd,:)&
         +dc41*dpf(ind,ipd,:)&
         +dc42*dpf(inu,ipu,:))
   enddo
   pf=pfout
   dpf=dpfout
   d2pf=d2pfout
   end subroutine correlatepair

end module orbital
