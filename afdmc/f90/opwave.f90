module wavefunction
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   real(kind=r8), private, parameter :: tiny=1.0e-8_r8
   integer(kind=i4), private, save :: npart,nbin(4)
   real(kind=r8), private, save :: el
   real(kind=r8), private, save, allocatable :: ak(:,:,:)
   logical, private, save :: dols
   integer(kind=i4), private, save :: npdet
   real(kind=r8), private, save :: rcut,acut,a
   complex(kind=r8), private, save, allocatable :: tau1(:,:),sig1(:,:),sigtau1(:,:,:)
   complex(kind=r8), private, save, allocatable :: tau2(:,:,:),sig2(:,:,:),sigtau2(:,:,:,:,:)
   complex(kind=r8), private, save, allocatable :: np0(:),np1(:),pp(:),nn(:)
   complex(kind=r8), private, save, allocatable :: tens(:),tenstau(:)
   logical, private, save :: noprot
   logical, private, save :: optv3,optcsb
contains
   subroutine setpsi(npartin,elin,hbar,vlsin,noprotin,nprot,nneut)
   use mympi
   use backflow
   use correlator
   use jastrow
   real(kind=r8) :: elin,rho,hbar
   integer(kind=i4) :: npartin,i,j,npair,nprot,nneut
   integer(kind=i4), allocatable :: iktemp(:,:,:)
   real(kind=r8) :: efg,kf,pi,fac,ak2
   logical :: vlsin,noprotin
   noprot=.false.
   allocate(iktemp(4,3,200)) ! this should be large enough!
   if (myrank().eq.0) then
      read (5,*) rho     !number density
      read (5,*) npart   !particle number
      read (5,*) nbin(1) !number of pup
      do i=1,nbin(1)
         read (5,*) iktemp(1,:,i)
      enddo
      read (5,*) nbin(2) !number of pdn
      do i=1,nbin(2)
         read (5,*) iktemp(2,:,i)
      enddo
      read (5,*) nbin(3) !number of nup
      do i=1,nbin(3)
         read (5,*) iktemp(3,:,i)
      enddo
      read (5,*) nbin(4) !number of ndn
      do i=1,nbin(4)
         read (5,*) iktemp(4,:,i)
      enddo
      read (5,*) rcut
      read (5,*) acut
      read (5,*) a
      if (rho.le.0.0_r8) then
         el=-2.0_r8*rho
         rho=npart/el**3
      endif
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
      write(6,'(''k values ='',(t35,3i5))') (iktemp(1,:,i),i=1,nbin(1))
      write (6,'(''proton down ='',t40,i10)') nbin(2)
      write(6,'(''k values ='',(t35,3i5))') (iktemp(2,:,i),i=1,nbin(2))
      write (6,'(''neutron up ='',t40,i10)') nbin(3)
      write(6,'(''k values ='',(t35,3i5))') (iktemp(3,:,i),i=1,nbin(3))
      write (6,'(''neutron down ='',t40,i10)') nbin(4)
      write(6,'(''k values ='',(t35,3i5))') (iktemp(4,:,i),i=1,nbin(4))
      write(6,'(''rcut in pair sampling ='',t40,f10.5)') rcut
      write(6,'(''acut in pair sampling ='',t40,f10.5)') acut
      write(6,'(''a constant in psig ='',t40,e12.5)') a
      if ((nbin(1)+nbin(2)).eq.0) then
!        noprot=.true. ! still need some work to calculate observables, FIX IT!
         write(6,'(''Pure neutron matter'')')
      endif
   endif
   call bcast(npart)
   call bcast(nbin)
   call bcast(iktemp)
   call bcast(rcut)
   call bcast(acut)
   call bcast(a)
   call bcast(noprot)
   call bcast(el)
   npartin=npart
   elin=el
   noprotin=noprot
   allocate(ak(4,3,maxval(nbin)))
   pi=4.0_r8*atan(1.0_r8)
   fac=2.0_r8*pi/el
   if (nbin(1).ne.0) ak(1,:,1:nbin(1))=fac*iktemp(1,:,1:nbin(1))
   if (nbin(2).ne.0) ak(2,:,1:nbin(2))=fac*iktemp(2,:,1:nbin(2))
   if (nbin(3).ne.0) ak(3,:,1:nbin(3))=fac*iktemp(3,:,1:nbin(3))
   if (nbin(4).ne.0) ak(4,:,1:nbin(4))=fac*iktemp(4,:,1:nbin(4))
   ak2=0.0_r8
   do i=1,4
      do j=1,nbin(i)
         ak2=ak2-sum(ak(i,:,j)**2)
      enddo
   enddo
   ak2=-hbar*ak2
   if (myrank().eq.0) then
      write(6,'(''calculating LS ='',t40,l10)') dols
      write(6,'(''energy of uncorrelated state ='',t40,f10.5)') ak2/npart
   endif
!  call setbf(el)
   dols=vlsin
   call setijas(1)
   call initcormod(npart,elin)
   npair=npart*(npart-1)/2
   allocate(tau1(3,npart),sig1(3,npart),sigtau1(3,3,npart))
   allocate(tau2(3,3,npair),sig2(3,3,npair),sigtau2(3,3,3,3,npair))
   allocate(np0(npair),np1(npair),pp(npair),nn(npair))
   allocate(tens(npair),tenstau(npair))
   nprot=nbin(1)+nbin(2)
   nneut=nbin(3)+nbin(4)
   end subroutine setpsi

   subroutine hpsi(w,doall)
   use stack
   use v6pot
   use random
   integer(kind=i4) :: i,j,is,it,kc
   real(kind=r8) :: ddx(3),r,vvls,vvlst,vvlspr,vvlstpr
   complex(kind=r8) :: psi,dpsi(3,npart),d2psi
   complex(kind=r8) :: dpsix(3,npart)
   complex(kind=r8) :: cvls,cvlst
   complex(kind=r8) :: cvlspr,cvlstpr
   type (walker) :: w
   complex(kind=r8) :: spsave(4,npart),spx(4,15,npart),p(3,3),pt(3,3)
   logical :: doall
   real(kind=r8) :: prob,rn(1)
   real(kind=r8) :: xsave(3,npart)
   call vnpsi2(w,.false.)
   if (doall) then
      xsave=w%x
      psi=w%psi
      dpsi=czero
      d2psi=czero
      cvls=czero
      cvlst=czero
      cvlspr=czero
      cvlstpr=czero
      call tpsi(w,dpsi,d2psi)
      if (dols) then
         call setrn(w%irn)
         spsave=w%sp
         spx=opmult(w%sp)
         do i=2,npart
            do j=1,i-1
               p=0.0_r8
               pt=0.0_r8
! note that neighbour boxes are not included in V_LS !!!
               ddx=w%x(:,i)-w%x(:,j)
               ddx=ddx-el*nint(ddx/el)
               r=sqrt(sum(ddx**2))
               call vls(r,vvls,vvlst,1)
               call vls(r,vvlspr,vvlstpr,2)
               if (noprot) then
                  vvls=vvls+vvlst
                  vvlspr=vvlspr+vvlstpr
                  vvlst=0.0_r8
                  vvlstpr=0.0_r8
               endif
               if (r.le.rcut) then
                  prob=1.0_r8
               else
                  prob=exp(-acut*(r-rcut))
               endif
               rn=randn(1)
               if (rn(1).lt.prob) then
                  do is=1,3
                     if (vvls.ne.0.0_r8.or.vvlspr.ne.0.0_r8) then
                        w%sp(:,i)=spx(:,is,i)
                        call gradpsione(w,dpsix(:,i),i)
                        call gradpsione(w,dpsix(:,j),j)
                        p(:,is)=dpsix(:,i)-dpsix(:,j)
                        w%sp(:,i)=spsave(:,i)
                        w%sp(:,j)=spx(:,is,j)
                        call gradpsione(w,dpsix(:,i),i)
                        call gradpsione(w,dpsix(:,j),j)
                        p(:,is)=p(:,is)+dpsix(:,i)-dpsix(:,j)
                        w%sp(:,j)=spsave(:,j)
                     endif
                     do it=1,3
                        if (vvlst.ne.0.0_r8.or.vvlstpr.ne.0.0_r8) then
                           kc=it+3*(is-1)+6
                           w%sp(:,i)=spx(:,kc,i)
                           w%sp(:,j)=spx(:,it+3,j)
                           call gradpsione(w,dpsix(:,i),i)
                           call gradpsione(w,dpsix(:,j),j)
                           pt(:,is)=pt(:,is)+dpsix(:,i)-dpsix(:,j)
                           w%sp(:,i)=spx(:,it+3,i)
                           w%sp(:,j)=spx(:,kc,j)
                           call gradpsione(w,dpsix(:,i),i)
                           call gradpsione(w,dpsix(:,j),j)
                           pt(:,is)=pt(:,is)+dpsix(:,i)-dpsix(:,j)
                           w%sp(:,i)=spsave(:,i)
                           w%sp(:,j)=spsave(:,j)
                        endif
                     enddo
                  enddo
                  cvls=cvls+vvls*(ddx(1)*(p(2,3)-p(3,2))+ddx(2)*(p(3,1)-p(1,3))+ddx(3)*(p(1,2)-p(2,1)))/prob
                  cvlst=cvlst+vvlst*(ddx(1)*(pt(2,3)-pt(3,2))+ddx(2)*(pt(3,1)-pt(1,3))+ddx(3)*(pt(1,2)-pt(2,1)))/prob
                  cvlspr=cvlspr+vvlspr*(ddx(1)*(p(2,3)-p(3,2))+ddx(2)*(p(3,1)-p(1,3))+ddx(3)*(p(1,2)-p(2,1)))/prob
                  cvlstpr=cvlstpr+vvlstpr*(ddx(1)*(pt(2,3)-pt(3,2))+ddx(2)*(pt(3,1)-pt(1,3))+ddx(3)*(pt(1,2)-pt(2,1)))/prob
               endif
            enddo
         enddo
         cvls=0.25_r8*ci*cvls/psi
         cvlst=0.25_r8*ci*cvlst/psi
         cvlspr=0.25_r8*ci*cvlspr/psi
         cvlstpr=0.25_r8*ci*cvlstpr/psi
         call savern(w%irn)
      endif
      w%x=xsave
      call vnpsi2(w,.true.)
      w%dpsi=dpsi
      w%d2psi=d2psi
      w%v8all(7)=cvls
      w%v8all(8)=cvlst
      w%v8allpr(7)=cvlspr
      w%v8allpr(8)=cvlstpr
   endif
   return
   end subroutine hpsi

   subroutine tpsi(w,dpsi,d2psi)
   use stack
   type (walker) :: w
   integer(kind=i4) :: i,ic
   complex(kind=r8) :: psi,dpsi(3,npart),d2psi
   complex(kind=r8) :: d2num,dnum
!  real(kind=r8), parameter :: dx=0.0005_r8
   real(kind=r8), parameter :: dx=0.01_r8
   psi=w%psi
   dpsi=czero
   d2psi=czero
   do i=1,npart
      do ic=1,3
         d2num=-2.0_r8*psi
         w%x(ic,i)=w%x(ic,i)+dx
         call vnpsi2(w,.false.)
         d2num=d2num+w%psi
         dnum=w%psi
         w%x(ic,i)=w%x(ic,i)-2.0_r8*dx
         call vnpsi2(w,.false.)
         d2num=d2num+w%psi
         dnum=dnum-w%psi
         w%x(ic,i)=w%x(ic,i)+dx
         d2psi=d2psi+d2num/dx**2/psi
         dpsi(ic,i)=dnum/2.0_r8/dx/psi
      enddo
   enddo
   w%psi=psi ! restore psi in the original positions
   end subroutine tpsi

   subroutine gradpsione(w,dpsi,i)
   use stack
   type (walker) :: w
   integer(kind=i4) :: i,ic
   complex(kind=r8) :: psi,dpsi(:)
   complex(kind=r8) :: dnum
!  real(kind=r8), parameter :: dx=0.0005_r8
   real(kind=r8), parameter :: dx=0.01_r8
   psi=w%psi
   do ic=1,3
      w%x(ic,i)=w%x(ic,i)+dx
      call vnpsi2(w,.false.)
      dnum=w%psi
      w%x(ic,i)=w%x(ic,i)-2.0_r8*dx
      call vnpsi2(w,.false.)
      dnum=dnum-w%psi
      w%x(ic,i)=w%x(ic,i)+dx
! careful, here I calculate dpsi, NOT dpsi/psi
      dpsi(ic)=dnum/2.0_r8/dx
   enddo
   w%psi=psi
   end subroutine gradpsione

   subroutine vnpsi2(w,dopot)
   use stack  ! to define walker type
   use jastrow
   use v6pot
   use matrixmod
   use correlator
   use v3bpot
   use random
   integer(kind=i4) :: i,k,kk,is,j
   real(kind=r8) :: u,akr,fls(npart,npart),utau
   real(kind=r8) :: fasig(3,npart,3,npart),fasigtau(3,npart,3,npart),fatau(npart,npart)
   real(kind=r8) :: fataupp(npart,npart),fataunn(npart,npart)
   real(kind=r8) :: fasigtautni(3,npart,3,npart),fatautni(npart,npart)
   real(kind=r8) :: fcsigtautni(3,npart,3,npart)
   real(kind=r8), dimension(npart,npart) :: v2,v3,v4
   real(kind=r8), dimension(3,npart,3,npart) :: v5,v6
   real(kind=r8), dimension(3,npart,3,npart) :: xpi
   real(kind=r8), dimension(3,npart,npart) :: g2s3b
   real(kind=r8), dimension(npart,npart) :: delta
   complex(kind=r8) :: ph(npart,4,npart)
   complex(kind=r8) :: als(4,3)
   complex(kind=r8) :: sxmallz(npart,4,npart),sxz(4,npart,npart)
   complex(kind=r8) :: det
   complex(kind=r8) :: tdet
   complex(kind=r8) :: smati(npart,npart)
   complex(kind=r8) :: cvs(6)
   complex(kind=r8) :: detrat
   real(kind=r8) :: dummy1
   type (walker) :: w
   logical :: dopot
   real(kind=r8) :: tnic,vc,u3c
   complex(kind=r8) :: tni2pia,tni2pitm,tni2pic,tni2picxd,tni2picdd,tnivd1,tnivd2,tnive,tni2piaxd,tni2piadd
   complex(kind=r8) :: tni2piapr,tni2piaxdpr,tni2piaddpr
   real(kind=r8) :: rij(3),cf3
   w%x=w%x-el*nint(w%x/el)
!  call addbf(w%x,xbig)
   call hspot(w%x,w%vc,v2,v3,v4,v5,v6,dopot,1) !always need vc
   w%vcoulc=0.0_r8
   w%vcoul=0.0_r8
   w%vext=0.0_r8
   if (noprot) then
      w%vc=w%vc+0.5_r8*sum(v2)
      v2=0.0_r8
      v3=v3+v4
      v4=0.0_r8
      v5=v5+v6
      v6=0.0_r8
   endif
   if (dopot) then
      call setxspx(w,3.0_r8*rcut,acut,noprot)
   endif
! TNI stuff
   call hstnipot(w%x,tnic,xpi,g2s3b,delta,dummy1,dummy1,dummy1,dummy1,dummy1,dummy1,dummy1,dummy1,dummy1,1.0_r8)
   w%tnic=tnic
!
   call hscor(npart,w%x,u,fasig,fasigtau,fatau,fls,fataupp,fataunn)
   u3c=0.0_r8
   call tniacor(npart,w%x,u3c,fasigtautni,fatautni,fcsigtautni,cf3)
   fasigtau=fasigtau+fasigtautni
   fatau=fatau+fatautni
   utau=0.0_r8
   if (noprot) then
      utau=0.5_r8*sum(fatau)
      fatau=0.0_r8
      fasig=fasig+fasigtau
      fasigtau=0.0_r8
   endif
   call calfop(fatau,fasig,fasigtau,fataupp,fataunn,fcsigtautni,cf3,w%sp,0.0_r8,dopot)
   ph=czero
   kk=0
   do is=1,4
      do k=1,nbin(is)
         kk=kk+1
         do i=1,npart
            akr=w%x(1,i)*ak(is,1,k)+w%x(2,i)*ak(is,2,k)+w%x(3,i)*ak(is,3,k)
            ph(kk,is,i)=cmplx(cos(akr),-sin(akr))
! add LS backflow correlations
            als=czero
            do j=1,npart
               if (j.ne.i) then
                  rij=w%x(:,i)-w%x(:,j)
                  rij=rij-el*nint(rij/el)
                  als(is,1)=als(is,1)+0.5_r8*fls(i,j)*(rij(2)*ak(is,3,k)-rij(3)*ak(is,2,k))
                  als(is,2)=als(is,2)+0.5_r8*fls(i,j)*(rij(3)*ak(is,1,k)-rij(1)*ak(is,3,k))
                  als(is,3)=als(is,3)+0.5_r8*fls(i,j)*(rij(1)*ak(is,2,k)-rij(2)*ak(is,1,k))
               endif
            enddo
            als=als*ph(kk,is,i)
! sigma_y has a minus sign because we use psi*
            ph(kk,1,i)=ph(kk,1,i)+als(2,1)+ci*als(2,2)+als(1,3)
            ph(kk,2,i)=ph(kk,2,i)+als(1,1)-ci*als(1,2)-als(2,3)
            ph(kk,3,i)=ph(kk,3,i)+als(4,1)+ci*als(4,2)+als(3,3)
            ph(kk,4,i)=ph(kk,4,i)+als(3,1)-ci*als(3,2)-als(4,3)
         enddo
      enddo
   enddo
   tdet=czero
   w%v8all=czero
   w%v8allpr=czero
   w%vcoul=czero
   do i=1,npart
      smati(:,i)=matmul(ph(:,:,i),w%sp(:,i))
   enddo
   call matinv(smati,det,npart)
   sxmallz=reshape(matmul(smati,reshape(ph(:,:,:),(/npart,4*npart/))),shape(sxmallz))
   sxz=reshape(transpose(reshape(sxmallz,(/npart,4*npart/))),shape(sxz))
   call cordet(detrat,sxz)
   tdet=det*detrat
      if (dopot) then
         call v6tot(w%x,w%sp,v2,v3,v4,v5,v6,cvs,tnic,tni2pia,tni2pitm,&
           tni2pic,tni2picxd,tni2picdd,tnivd1,tnivd2,tnive,tni2piaxd,&
           tni2piadd,tni2piapr,tni2piaxdpr,tni2piaddpr)
         cvs=det*cvs/tdet
         call getop1(sig1,tau1,sigtau1)
         sig1=det*sig1/tdet
         tau1=det*tau1/tdet
         sigtau1=det*sigtau1/tdet
         call getop2(sig2,tau2,sigtau2,np0,np1,pp,nn)
         sig2=det*sig2/tdet
         tau2=det*tau2/tdet
         sigtau2=det*sigtau2/tdet
         np0=det*np0/tdet
         np1=det*np1/tdet
         pp=det*pp/tdet
         nn=det*nn/tdet
         call getoptens(w%x,tens,tenstau)
         tens=det*tens/tdet
         tenstau=det*tenstau/tdet
         w%tni2pia=det*tni2pia/tdet
         w%tni2piaxd=det*tni2piaxd/tdet
         w%tni2piadd=det*tni2piadd/tdet
         w%tni2pic=det*tni2pic/tdet
         w%tni2picxd=det*tni2picxd/tdet
         w%tni2picdd=det*tni2picdd/tdet
         w%tni2pitm=det*tni2pitm/tdet
         w%tnivd1=det*tnivd1/tdet
         w%tnivd2=det*tnivd2/tdet
         w%tnive=det*tnive/tdet
!print*,w%x
!print*,'!!!!!!!'
!print*,w%sp
!print*,'!!!!!!!'
!print*,'vc=',w%tnic
!print*,'2pia=',w%tni2pia
!print*,'2pixd=',w%tni2piaxd
!print*,'2pidd=',w%tni2piadd
!print*,'tm=',w%tni2pitm
!print*,'vd1=',w%tnivd1
!print*,'vd2=',w%tnivd2
!print*,'vde=',w%tnive
!print*,'!!!!!!!'
!print*,'2pitot=',w%tni2pia+w%tni2piaxd+w%tni2piadd
!print*,'vdtot=',w%tnivd1+w%tnivd2
!stop
         w%tni2piapr=det*tni2piapr/tdet
         w%tni2piaxdpr=det*tni2piaxdpr/tdet
         w%tni2piaddpr=det*tni2piaddpr/tdet
         w%v8all(1)=w%vc
         w%v8all(2:6)=cvs(2:6)
         call hspot(w%x,vc,v2,v3,v4,v5,v6,dopot,2)
         call calpot(cvs,v2,v3,v4,v5,v6)
         cvs=det*cvs/tdet
         w%v8allpr(1)=vc
         w%v8allpr(2:6)=cvs(2:6)
      endif
   w%psi=(tdet+utau*det)*exp(-u)*(1.0_r8+u3c)
   w%psig=abs(real(w%psi))+a*abs(aimag(w%psi))
   end subroutine vnpsi2

   subroutine getsigtauop(tau1in,sig1in,sigtau1in,tau2in,sig2in,sigtau2in,np0in,np1in,ppin,nnin,tensin,tenstauin)
   complex(kind=r8) :: tau1in(:,:),sig1in(:,:),sigtau1in(:,:,:)
   complex(kind=r8) :: tau2in(:,:,:),sig2in(:,:,:),sigtau2in(:,:,:,:,:)
   complex(kind=r8) :: np0in(:),np1in(:),ppin(:),nnin(:)
   complex(kind=r8) :: tensin(:),tenstauin(:)
   tau1in=tau1
   sig1in=sig1
   sigtau1in=sigtau1
   tau2in=tau2
   sig2in=sig2
   sigtau2in=sigtau2
   np0in=np0
   np1in=np1
   ppin=pp
   nnin=nn
   tensin=tens
   tenstauin=tenstau
   end subroutine getsigtauop

   function opmult(sp)
   complex(kind=r8) :: sp(:,:),opmult(4,15,npart)
!
! The order is 1-3 sx,sy,sz, 4-6 tx,ty,tx, 7-9 sx*(tx,ty,tz)
! 10-12 sy*(tx,ty,tz), 13-15 sz*(tx,ty,tz)
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
   integer(kind=i4) :: ic
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
      if (ic.eq.3) tz=-(psip-psim)/(ci*2.0_r8*dphi*psi)
   enddo
   t2=aj
   aj=czero
   jz=czero
   j2=czero
   w%x=xsave
   w%sp=ssave
   end subroutine checkj

   subroutine chkpot(w)
   use stack
   use v6pot
   type (walker) :: w
   complex(kind=r8) :: spx(4,15,npart),psi,spsave(4,npart)
   complex(kind=r8) :: vsig,vtau,vsigtau
   real(kind=r8) :: vc
   real(kind=r8) :: asig(3,npart,3,npart),asigtau(3,npart,3,npart),atau(npart,npart)
   real(kind=r8), dimension(npart,npart) :: v2,v3,v4
   real(kind=r8), dimension(3,npart,3,npart) :: v5,v6
   integer(kind=i4) :: i,j,ic,jc,kc,it
   call vnpsi2(w,.true.)
   spsave=w%sp
   psi=w%psi
   spx=opmult(w%sp)
   call hspot(w%x,vc,v2,v3,v4,v5,v6,.true.,1)
   asig=0.0_r8
   asigtau=0.0_r8
   atau=0.0_r8
   do i=1,npart
      do j=1,npart
         atau(i,j)=atau(i,j)+v2(i,j)
         do ic=1,3
            asig(ic,i,ic,j)=asig(ic,i,ic,j)+v3(i,j)
            asigtau(ic,i,ic,j)=asigtau(ic,i,ic,j)+v4(i,j)
            do jc=1,3
               asig(ic,i,jc,j)=asig(ic,i,jc,j)+v5(ic,i,jc,j)
               asigtau(ic,i,jc,j)=asigtau(ic,i,jc,j)+v6(ic,i,jc,j)
            enddo
         enddo
      enddo
   enddo
   vsig=czero
   vtau=czero
   vsigtau=czero
   do i=2,npart
      do j=1,i-1
         do ic=1,3
            w%sp(:,i)=spx(:,ic,i)
            call vnpsi2(w,.true.)
            vsig=vsig+sum(asig(:,j,ic,i)*sig1(:,j))*w%psi
            w%sp(:,i)=spx(:,ic+3,i)
            call vnpsi2(w,.true.)
            vtau=vtau+atau(j,i)*tau1(ic,j)*w%psi
            do it=1,3
               kc=it+3*(ic-1)+6
               w%sp(:,i)=spx(:,kc,i)
               call vnpsi2(w,.true.)
               do jc=1,3
                  vsigtau=vsigtau+asigtau(jc,j,ic,i)*sigtau1(jc,it,j)*w%psi
               enddo
               w%sp(:,i)=spsave(:,i)
            enddo
         enddo
      enddo
   enddo
   vsig=vsig/psi
   vtau=vtau/psi
   vsigtau=vsigtau/psi
   write(6,'(''chkpot, v= '',6e15.7)') vc+vtau+vsig+vsigtau
   call vnpsi2(w,.true.)
   write(6,'(''hpsi,   v= '',6e15.7)') sum(w%v8all(1:6))
   end subroutine chkpot

   subroutine chkls(w)
   use stack
   use v6pot
   type (walker) :: w
   complex(kind=r8) :: spsave(4,npart),psi,dpsi(3,npart)
   complex(kind=r8) :: spx(4,15,npart),cvls,cvlstau,p(3,3)
   real(kind=r8) :: dx(3),r,vvls,vvlstau
   integer(kind=i4) :: i,j,is
!  call vnpsi2(w,.false.,.false.)
! careful, here I need the derivatives already stored in stack!
   spsave=w%sp
   psi=w%psi
   dpsi=w%dpsi
   spx=opmult(w%sp)
   cvls=czero
   cvlstau=czero
   do i=2,npart
      do j=1,i-1
         dx=w%x(:,i)-w%x(:,j)
         r=sqrt(sum(dx**2))
         call vls(r,vvls,vvlstau,1)
         do is=1,3
            w%sp(:,i)=spx(:,is,i)
            call vnpsi2(w,.false.)
            p(:,is)=(dpsi(:,i)-dpsi(:,j))*w%psi
            w%sp(:,i)=spsave(:,i)
            w%sp(:,j)=spx(:,is,j)
            call vnpsi2(w,.false.)
            p(:,is)=p(:,is)+(dpsi(:,i)-dpsi(:,j))*w%psi
            w%sp(:,j)=spsave(:,j)
         enddo
         cvls=cvls+vvls*(dx(1)*(p(2,3)-p(3,2))+dx(2)*(p(3,1)-p(1,3))+dx(3)*(p(1,2)-p(2,1)))
      enddo
   enddo
   cvls=0.25_r8*ci*cvls/psi
!  write(6,'(''chkls, vls= '',6e15.7)') cvls
! cvlstau to be done!!!!!!!!!!!!!
   end subroutine chkls

   subroutine chkop(w)
   use stack
   type (walker) :: w
   complex(kind=r8) :: spsave(4,npart),spx(4,15,npart),psi
   complex(kind=r8) :: sig,tau,sigtau(3)
   integer(kind=i4) :: i,ic,jc,kc
   spsave=w%sp
   spx=opmult(w%sp)
   do i=1,npart
      do ic=1,3
         call vnpsi2(w,.true.)
         psi=w%psi
         sig=sig1(ic,i)
         tau=tau1(ic,i)
         sigtau=sigtau1(ic,:,i)
         w%sp(:,i)=spx(:,ic,i)
         call vnpsi2(w,.true.)
         w%sp(:,i)=spsave(:,i)
         write(6,'(''sigma ana ic '',t25,i5,2f15.10)') ic,sig
         write(6,'(''sigma num ic '',t25,i5,2f15.10)') ic,w%psi/psi
         w%sp(:,i)=spx(:,ic+3,i)
         call vnpsi2(w,.true.)
         w%sp(:,i)=spsave(:,i)
         write(6,'(''tau ana ic '',t25,i5,2f15.10)') ic,tau
         write(6,'(''tau num ic '',t25,i5,2f15.10)') ic,w%psi/psi
         do jc=1,3
            kc=jc+3*(ic-1)+6
            w%sp(:,i)=spx(:,kc,i)
            call vnpsi2(w,.true.)
            w%sp(:,i)=spsave(:,i)
            write(6,'(''sigtau ana ic jc '',t20,2i5,2f15.10)') ic,jc,sigtau(jc)
            write(6,'(''sigtau num ic jc '',t20,2i5,2f15.10)') ic,jc,w%psi/psi
         enddo
      enddo
   enddo
   end subroutine chkop

   subroutine addbf(x,xbig)
   use backflow
   real(kind=r8) :: x(:,:),xbig(:,:)
   real(kind=r8) :: r,dx(3),bf,dbf,d2bf
   integer(kind=i4) :: i,j
   xbig=x
   do i=1,npart-1
      do j=i+1,npart
         dx=x(:,i)-x(:,j)
         dx=dx-el*nint(dx/el)
         r=sqrt(sum(dx**2))
         call getbf(r,bf,dbf,d2bf)
         xbig(:,i)=xbig(:,i)+bf*dx
         xbig(:,j)=xbig(:,j)-bf*dx
      enddo
   enddo
   end subroutine addbf

   subroutine setoptv3(optv3in)
   logical :: optv3in
   optv3=optv3in
   end subroutine setoptv3
 
   subroutine setoptcsb(optcsbin)
   logical :: optcsbin
   optcsb=optcsbin
   end subroutine setoptcsb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the following are subrotuines needed by the optimization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine setdetparam(params)
   use backflow
   real(kind=r8) :: params(:),pbf(4)
   pbf=params
   call initbf(pbf)
   end subroutine setdetparam

   subroutine getdetparam(nparam,params)
   use backflow
   integer(kind=i4) :: nparam
   real(kind=r8), pointer :: params(:)
   real(kind=r8) :: pbf(4)
   nparam=4
   allocate(params(nparam))
   call getbfpar(pbf)
   params=pbf
   npdet=nparam
   end subroutine getdetparam

   subroutine getderpsi(w,dpsi)
   use stack
   use backflow
   type (walker) :: w
   real(kind=r8) :: dpsi(:)
   real(kind=r8) :: pbf(4)
   complex(kind=r8) :: psi0
   integer(kind=i4) :: i
   real(kind=r8), parameter :: dp=0.001_r8
   dpsi=0.0_r8
   call getbfpar(pbf)
   call hpsi(w,.false.)
   psi0=w%psi
   do i=1,npdet
      pbf(i)=pbf(i)+dp
      call initbf(pbf)
      call hpsi(w,.false.)
      dpsi(i)=(w%psi-psi0)/dp
      pbf(i)=pbf(i)-dp
      call initbf(pbf)
   enddo
   call hpsi(w,.false.)
   dpsi=dpsi/w%psi
   end subroutine getderpsi

   subroutine getjasder(w,dpsi)
   use stack
   use jastrow
   type (walker) :: w
   real(kind=r8) :: dpsi(:)
   real(kind=r8), parameter :: dp=0.0025_r8
   complex(kind=r8) :: psi0
   integer(kind=i4) :: i,ipar
   real(kind=r8) :: rscpp,rscnn
   dpsi=0.0_r8
   call setijas(1)
   call hpsi(w,.false.)
   psi0=w%psi
   ipar=1
   do i=1,size(dpsi)-17 ! the last parameters are q1c,q2c,q1p,q2p,rscal,a3c,a3a,a3atm,a3avd1,a3avd2,a3avdc3,a3avec3,a3ve,rscpp,rscnn,bpp,bnn
      call setijas(i+1)
      call hpsi(w,.false.)
      dpsi(i)=(w%psi-psi0)/dp
      call setijas(1)
      ipar=ipar+1
   enddo
   call setijas(1)
   do i=1,4 ! q1c,q2c,q1p,q2p
      call f3p(i)
      call hpsi(w,.false.)
      call f3m(i)
      dpsi(ipar)=(w%psi-psi0)/dp
      ipar=ipar+1
   enddo
   if (optv3) then
      do i=5,13 ! rsctni(1),ptni(8)
         call f3p(i)
         call hpsi(w,.false.)
         call f3m(i)
         dpsi(ipar)=(w%psi-psi0)/dp
         ipar=ipar+1
      enddo
   else
      ipar=ipar+9
   endif
   if (optcsb) then
! now do rscpp and rscnn
      call getcsb(rscpp,rscnn)
      rscpp=rscpp+dp
      call setcsb(rscpp,rscnn)
      call hpsi(w,.false.)
      rscpp=rscpp-dp
      call setcsb(rscpp,rscnn)
      dpsi(ipar)=(w%psi-psi0)/dp
      ipar=ipar+1  
      rscnn=rscnn+dp
      call setcsb(rscpp,rscnn)
      call hpsi(w,.false.)
      rscnn=rscnn-dp
      call setcsb(rscpp,rscnn)
      dpsi(ipar)=(w%psi-psi0)/dp
      ipar=ipar+1  
! and finally bpp and bnn
      do i=16,17 ! ugly, think how to fix this!
         call setijas(i+1)
         call hpsi(w,.false.)
         dpsi(ipar)=(w%psi-psi0)/dp
         call setijas(1)
         ipar=ipar+1
      enddo
   endif
   call hpsi(w,.false.)
   dpsi=dpsi/w%psi
   end subroutine getjasder
end module wavefunction
