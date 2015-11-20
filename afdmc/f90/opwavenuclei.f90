module wavefunction
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   real(kind=r8), private, parameter :: tiny=1.0e-8_r8
   integer(kind=i4), private, save :: npart,ndet,npair
   complex(kind=r8), private, save, allocatable :: cdet(:),cdet0(:)
   integer(kind=i4), private, save :: nconf
   integer(kind=i4), private, save, allocatable :: norb(:)
   real(kind=r8), private, save, allocatable :: vorb(:)
   complex(kind=r8), private, save, allocatable :: totdet(:)
   logical, private, save :: dols
   integer(kind=i4), private, save :: npdet
   real(kind=r8), private, save :: rcut,acut,a
   complex(kind=r8), private, save, allocatable :: tau1(:,:),sig1(:,:),sigtau1(:,:,:)
   complex(kind=r8), private, save, allocatable :: tau2(:,:,:),sig2(:,:,:),sigtau2(:,:,:,:,:)
   complex(kind=r8), private, save, allocatable :: np0(:),np1(:),pp(:),nn(:)
   complex(kind=r8), private, save, allocatable :: tens(:),tenstau(:)
   logical, private, save :: noprot
   logical, private, save :: optv3,optcsb
   integer(kind=i4), private, save :: vexid
   real(kind=r8), private, save :: wsv0,wsrad,wss,hbari,hbarom
contains
   subroutine setpsi(npartin,elin,hbar,vlsin,noprotin,nprot,nneut)
   use mympi
   use phimod
   use correlator
   use jastrow
   integer(kind=i4) :: npartin,nprot,nneut
   real(kind=r8) :: elin,hbar
   character(len=30) :: orbnow
   character(len=80) :: line
   character(len=80), allocatable :: lines(:,:)
   character(len=30), allocatable :: orbfil(:)
   real(kind=r8), allocatable :: orb(:,:)
   integer(kind=i4) :: jmax,nrad,notab,j2,m2,ll,itau
   integer(kind=i4) :: i,j,k
   integer(kind=i4), allocatable :: lorb(:),iorb(:,:),jpval(:,:),mpval(:,:),lval(:,:),isoval(:,:)
   real(kind=r8) :: romax,dr,dummy
   real(kind=r8), allocatable :: fscal(:)
   logical :: rdiv,vlsin,noprotin
   npart=npartin
   noprot=.false.
   noprotin=noprot
   nprot=0
   nneut=0
   if (myrank().eq.0) then
      read (5,*) notab   !table points for orbitals
      read (5,*) romax   !rmax for orbitals
      read (5,*) rdiv    !divido orbitali per r?
      read (5,*) npart   !particle number
      read (5,*) vexid   !external potential
      if (vexid.eq.1) then
         read (5,*) wsv0    ! wood-saxon strength (this should be negative)
         read (5,*) wsrad   ! wood-saxon radius
         read (5,*) wss     ! wood-saxon smoothness (length)
      else if (vexid.eq.2) then
         read (5,*) hbarom
         hbari=1.0_r8/(2.0_r8*hbar)
      endif
      read (5,*) rcut
      read (5,*) acut
      read (5,*) a
      read (5,*) nrad    !total number of radial functions for run
      allocate(orb(notab,nrad),orbfil(nrad),lorb(nrad),fscal(nrad))
!
! each radial function is given by the file name followed by
! the angular momentum L value
!
      do i=1,nrad
         read (5,'(a80)') line
         line=adjustl(line)
         orbfil(i)=line(1:index(line,' ')-1)
         read(line(index(line,' '):80),*) lorb(i),fscal(i)
         open(unit=22,file=orbfil(i),status='old')
         do j=1,notab
            read (22,*) dr,orb(j,i)
            if (rdiv) orb(j,i)=orb(j,i)/max(dr,0.00001_r8)
         enddo
         read (22,*,end=10,err=10) dummy
         write (6,'(''too many lines in radial file'')')
         call abort
     10  close(22)
      enddo
      read (5,*) ndet   !number of determinants
      read (5,*) nconf  !number of components
      allocate(norb(nconf),vorb(nconf))
      do i=1,nconf
         read (5,*) norb(i),vorb(i) ! degeneracy and amplitude
      enddo
      if (sum(norb).ne.ndet) then
         write (6,'(''Wrong total number of determinants'')')
         write (6,*) ndet,norb
         call abort
      endif         
      allocate(jpval(npart,ndet),mpval(npart,ndet),lval(npart,ndet))
      allocate(isoval(npart,ndet),iorb(npart,ndet),cdet(ndet),cdet0(ndet))
      allocate(lines(npart,ndet))
      do i=1,ndet
         read (5,*) cdet0(i) ! coefficient
         do j=1,npart
            read (5,'(a80)') line
            lines(j,i)=line
            read (line,*) j2,m2,ll,itau ! 2*J,2*M,L,2*tau then file
            line=adjustl(line)
            do k=1,4
               line=line(index(line,' '):80)
               line=adjustl(line)
            enddo
            orbnow=line(1:index(line,' ')-1)
            do k=1,nrad
               iorb(j,i)=k
               if (orbnow.eq.orbfil(k)) exit
               if (k.eq.nrad) then
                  write(6,'(''Orbital file not found '',a30)') orbnow
                  call abort
               endif
            enddo
            if (mod(j2,2).ne.1.or.mod(abs(m2),2).ne.1.or.abs(j2-2*ll).ne.1 &
               .or.abs(m2).gt.j2.or.abs(itau).ne.1 &
               .or.ll.ne.lorb(iorb(j,i))) then
               write (6,'(''Orbital quantum numbers wrong '',6i5)') &
                  j2,m2,ll,iorb(j,i),lorb(iorb(j,i)),itau
               call abort
            endif
            jpval(j,i)=(j2+1)/2     ! jp = J+1/2
            mpval(j,i)=(j2+m2+2)/2  ! mp = J+M+1
            lval(j,i)=(2*ll-j2+1)/2-1 ! -1 for L=J-1/2 , 0 for L=J+1/2
            isoval(j,i)=2-itau  ! 3 for neutron, 1 for proton
            if (i.eq.1) then
               if (itau.eq.1) then
                  nprot=nprot+1
               else 
                  nneut=nneut+1
               endif
            endif
         enddo
      enddo
      jmax=maxval(jpval)
      do j=1,nrad 
         write (6,'(''scale orbital '',i3,'' with factor '',f10.5)') j,fscal(j)
      enddo
      do j=1,ndet
         write (6,'(''determinant #'',i2)') j
         write (6,'(''coefficient'',t30,2f10.5)') cdet0(j)
         do i=1,npart
            write (6,'(a80)') lines(i,j)
         enddo
      enddo
      k=1
      do i=1,nconf
         do j=1,norb(i)
            cdet(k)=vorb(i)*cdet0(k)
            write (6,'(''rescaled coefficient'',t30,2f10.5)') cdet(k)
            k=k+1
         enddo
      enddo
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      select case (vexid)
         case (0)
            write (6,'(''No external potential'')')
         case (1)
            write (6,'(''external WS V0 ='',t40,f10.5)') wsv0
            write (6,'(''external WS width ='',t40,f10.5)') wsrad
            write (6,'(''external WS smoothness ='',t40,f10.5)') wss
         case (2)
            write (6,'(''external harmonic hbar*omega ='',t40,f10.5)') wsv0
      end select
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''rcut in pair sampling ='',t40,f10.5)') rcut
      write (6,'(''acut in pair sampling ='',t40,f10.5)') acut
      write (6,'(''a constant in psig ='',t40,f10.5)') a
   endif
   call bcast(npart)
   call bcast(vexid)
   call bcast(wsv0)
   call bcast(wsrad)
   call bcast(wss)
   call bcast(hbarom)
   call bcast(hbari)
   call bcast(rcut)
   call bcast(acut)
   call bcast(a)
   call bcast(notab)
   call bcast(romax)
   call bcast(nrad)
   call bcast(jmax)
   call bcast(ndet)
   call bcast(nconf)
   if (myrank().ne.0) then
      allocate(norb(nconf),vorb(nconf))
      allocate(jpval(npart,ndet),mpval(npart,ndet),lval(npart,ndet))
      allocate(isoval(npart,ndet),iorb(npart,ndet),cdet(ndet),cdet0(ndet))
      allocate(orb(notab,nrad),lorb(nrad),fscal(nrad))
   endif
   call bcast(cdet0)
   call bcast(cdet)
   call bcast(jpval)
   call bcast(mpval)
   call bcast(lval)
   call bcast(isoval)
   call bcast(iorb)
   call bcast(lorb)
   call bcast(fscal)
   call bcast(orb)
   call bcast(norb)
   call bcast(vorb)
   call setphi(notab,jmax,nrad,npart,ndet,jpval,mpval,lval,isoval, &
      iorb,lorb,orb,romax,fscal)
   npartin=npart
   elin=-romax
   allocate(totdet(ndet))
   dols=vlsin
   call setijas(1)
   call initcormod(npart,elin)
   npair=npart*(npart-1)/2
   allocate(tau1(3,npart),sig1(3,npart),sigtau1(3,3,npart))
   allocate(tau2(3,3,npair),sig2(3,3,npair),sigtau2(3,3,3,3,npair))
   allocate(np0(npair),np1(npair),pp(npair),nn(npair))
   allocate(tens(npair),tenstau(npair))
   call bcast(nprot)
   call bcast(nneut)
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
         spsave=w%sp
         spx=opmult(w%sp)
         do i=2,npart
            do j=1,i-1
               p=0.0_r8
               pt=0.0_r8
               ddx=w%x(:,i)-w%x(:,j)
               r=sqrt(sum(ddx**2))
               call vls(r,vvls,vvlst,1)
               call vls(r,vvlspr,vvlstpr,2)
               if (r.le.rcut) then
                  prob=1.0_r8
               else
                  prob=exp(-acut*(r-rcut))
               endif
               call setrn(w%irn)
               rn=randn(1)
               call savern(w%irn)
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
   real(kind=r8), parameter :: dx=0.0005_r8
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
   real(kind=r8), parameter :: dx=0.0005_r8
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
   use phimod
   use correlator
   use v3bpot
   use random
   integer(kind=i4) :: i,j,k,idet
   real(kind=r8) :: eni
   real(kind=r8) :: u
   real(kind=r8) :: fasig(3,npart,3,npart),fasigtau(3,npart,3,npart),fatau(npart,npart)
   real(kind=r8) :: fataupp(npart,npart),fataunn(npart,npart)
   real(kind=r8) :: fasigtautni(3,npart,3,npart),fatautni(npart,npart)
   real(kind=r8) :: xpi(3,npart,3,npart),xd(3,npart,3,npart)
   real(kind=r8) :: fls(npart,npart)
   real(kind=r8), dimension(npart,npart) :: v2,v3,v4
   real(kind=r8), dimension(3,npart,3,npart) :: v5,v6
   real(kind=r8), dimension(3,npart,npart) :: g2s3b
   real(kind=r8), dimension(npart,npart) :: delta
   complex(kind=r8) :: ph(npart,4,npart,ndet),dph(npart,3,4,npart,ndet)
   complex(kind=r8) :: als(4,3,ndet)
   complex(kind=r8) :: sxmallz(npart,4,npart),sxz(4,npart,npart)
   complex(kind=r8) :: det
   complex(kind=r8) :: tdet
   complex(kind=r8) :: smati(npart,npart),smat(npart,npart)
   complex(kind=r8) :: cvs(6)
   complex(kind=r8) :: detrat,ccvs(6),ccvspr(6)
   complex(kind=r8) :: vcoul
   real(kind=r8) :: vcoulc
   complex(kind=r8) :: cvcoul
   real(kind=r8) :: cvcoulc
   complex(kind=r8) :: csig1(3,npart),csigtau1(3,3,npart),ctau1(3,npart)
   complex(kind=r8) :: csig2(3,3,npair),csigtau2(3,3,3,3,npair),ctau2(3,3,npair)
   complex(kind=r8) :: t(npair),ttau(npair)
   complex(kind=r8) :: np0a(npair),np1a(npair),ppa(npair),nna(npair)
   real(kind=r8) :: vcoulmat(npart,npart)
   real(kind=r8) :: phig
   type (walker) :: w
   logical :: dopot
   real(kind=r8) :: xcm(3),v,r
   real(kind=r8) :: tnic,dummy1,vc,u3c
   complex(kind=r8) :: tni2pia,tni2pitm,tni2pic,tni2picxd,tni2picdd,tnivd1,tnivd2,tnive,tni2piaxd,tni2piadd
   complex(kind=r8) :: ttni2pia,ttni2pitm,ttni2pic,ttni2picxd,ttni2picdd,ttnivd1,ttnivd2,ttnive,ttni2piaxd,ttni2piadd
   complex(kind=r8) :: tni2piapr,ttni2piapr,tni2piaxdpr,ttni2piaxdpr,tni2piaddpr,ttni2piaddpr
   real(kind=r8) :: rij(3),cf3xx,cf3xd,cf3dd
   real(kind=r8) :: rscal(3) ! match the dimension as in jastrowtabop
   rscal=1.0_r8
   eni=1.0_r8/npart
   do i=1,3
      xcm(i)=sum(w%x(i,:))*eni
      w%x(i,:)=w%x(i,:)-xcm(i)
   enddo
   call hspot(w%x,w%vc,v2,v3,v4,v5,v6,dopot,1) !always need vc
   vcoulmat=0.0_r8
   w%vcoulc=0.0_r8
   call getvcsb(w%x,vcoulmat)
   w%vcoulc=0.5_r8*0.25_r8*sum(vcoulmat)
   w%vext=0.0_r8
   if (vexid.ne.0) then
      do i=1,npart
         r=sqrt(sum(w%x(:,i)**2))
         call vpotex(r,v)
         w%vext=w%vext+v
      enddo
   endif
   if (dopot) then
      call setxspx(w,3.0_r8*rcut,acut,noprot)
   endif
! TNI stuff
   call hstnipot(w%x,tnic,xpi,g2s3b,delta,dummy1,dummy1,dummy1,dummy1,dummy1,dummy1,dummy1,dummy1,dummy1,rscal)
   w%tnic=tnic
!
   call hscor(npart,w%x,u,fasig,fasigtau,fatau,fls,fataupp,fataunn)
   u3c=0.0_r8
   call tniacor(npart,w%x,u3c,fasigtautni,fatautni,xpi,xd,cf3xx,cf3xd,cf3dd)
   fasigtau=fasigtau+fasigtautni
   fatau=fatau+fatautni
   call calfop(fatau,fasig,fasigtau,fataupp,fataunn,xpi,xd,cf3xx,cf3xd,cf3dd,w%sp,0.0_r8,dopot)
   do i=1,npart
      call getphi(w%x(:,i),ph(:,:,i,:),dph(:,:,:,i,:))
   enddo
   do i=1,npart
      do k=1,npart
! add linear LS single particle correlations   
         als=czero
         do j=1,npart
            if (j.ne.i.and.fls(i,j).ne.0.0_r8) then
               rij=w%x(:,i)-w%x(:,j)
               als(:,1,:)=als(:,1,:)+fls(i,j)*(rij(2)*dph(k,3,:,i,:)-rij(3)*dph(k,2,:,i,:))
               als(:,2,:)=als(:,2,:)+fls(i,j)*(rij(3)*dph(k,1,:,i,:)-rij(1)*dph(k,3,:,i,:))
               als(:,3,:)=als(:,3,:)+fls(i,j)*(rij(1)*dph(k,2,:,i,:)-rij(2)*dph(k,1,:,i,:))
            endif
         enddo
         als=-ci*als
! sigma_y has a minus sign because we use psi*
         ph(k,1,i,:)=ph(k,1,i,:)+als(2,1,:)+ci*als(2,2,:)+als(1,3,:)
         ph(k,2,i,:)=ph(k,2,i,:)+als(1,1,:)-ci*als(1,2,:)-als(2,3,:)
         ph(k,3,i,:)=ph(k,3,i,:)+als(4,1,:)+ci*als(4,2,:)+als(3,3,:)
         ph(k,4,i,:)=ph(k,4,i,:)+als(3,1,:)-ci*als(3,2,:)-als(4,3,:)
      enddo
   enddo
!
   tdet=czero
   w%v8all=czero
   w%v8allpr=czero
   w%vcoul=czero
   ccvs=czero
   ccvspr=czero
   vcoul=czero
   vcoulc=czero
   tau1=czero
   sig1=czero
   sigtau1=czero
   tau2=czero
   sig2=czero
   sigtau2=czero
   np0=czero
   np1=czero
   pp=czero
   nn=czero
   tens=czero
   tenstau=czero
   ttni2pia=czero
   ttni2pic=czero
   ttni2picxd=czero
   ttni2picdd=czero
   ttni2pitm=czero
   ttnivd1=czero
   ttnivd2=czero
   ttnive=czero
   ttni2piaxd=czero
   ttni2piadd=czero
   ttni2piapr=czero
   ttni2piaxdpr=czero
   ttni2piaddpr=czero
   do idet=1,ndet
      do i=1,npart
         smati(:,i)=matmul(ph(:,:,i,idet),w%sp(:,i))
      enddo
      if (idet.eq.1) smat=smati
      call matinv(smati,det,npart) !smati is now the inverse of the earlier smati
      det=cdet(idet)*det
      sxmallz=reshape(matmul(smati,reshape(ph(:,:,:,idet) &
           ,(/npart,4*npart/))),shape(sxmallz))
      sxz=reshape(transpose(reshape(sxmallz,(/npart,4*npart/))),shape(sxz))
      call cordet(detrat,sxz)
      !call cordet(detrat,sxz,w%sp)
      tdet=tdet+det*detrat
      totdet(idet)=cdet0(idet)*det*detrat/cdet(idet)
      if (dopot) then
         call v6tot(w%x,w%sp,v2,v3,v4,v5,v6,cvs,tnic,tni2pia,tni2pitm,&
           tni2pic,tni2picxd,tni2picdd,tnivd1,tnivd2,tnive,tni2piaxd,&
           tni2piadd,tni2piapr,tni2piaxdpr,tni2piaddpr)
         ccvs=ccvs+det*cvs
         call getop1(csig1,ctau1,csigtau1)
         sig1=sig1+det*csig1
         tau1=tau1+det*ctau1
         sigtau1=sigtau1+det*csigtau1
         call getop2(csig2,ctau2,csigtau2,np0a,np1a,ppa,nna)
         sig2=sig2+det*csig2
         tau2=tau2+det*ctau2
         sigtau2=sigtau2+det*csigtau2
         np0=np0+det*np0a
         np1=np1+det*np1a
         pp=pp+det*ppa
         nn=nn+det*nna
         call getoptens(w%x,t,ttau)
         tens=tens+det*t
         tenstau=tenstau+det*ttau
         ttni2pia=ttni2pia+det*tni2pia
         ttni2piaxd=ttni2piaxd+det*tni2piaxd
         ttni2piadd=ttni2piadd+det*tni2piadd
         ttni2pic=ttni2pic+det*tni2pic
         ttni2picxd=ttni2picxd+det*tni2picxd
         ttni2picdd=ttni2picdd+det*tni2picdd
         ttni2pitm=ttni2pitm+det*tni2pitm
         ttnivd1=ttnivd1+det*tnivd1
         ttnivd2=ttnivd2+det*tnivd2
         ttnive=ttnive+det*tnive
         ttni2piapr=ttni2piapr+det*tni2piapr
         ttni2piaxdpr=ttni2piaxdpr+det*tni2piaxdpr
         ttni2piaddpr=ttni2piaddpr+det*tni2piaddpr
         call hspot(w%x,vc,v2,v3,v4,v5,v6,dopot,2)
         call calpot(cvs,v2,v3,v4,v5,v6)
         ccvspr=ccvspr+det*cvs
         call calpotcoul(cvcoul,cvcoulc,vcoulmat)
         vcoul=vcoul+det*cvcoul
      endif
   enddo
   if (dopot) then
      ccvs=ccvs/tdet
      w%v8all(1)=w%vc
      w%v8all(2:6)=ccvs(2:6)
      tau1=tau1/tdet
      sig1=sig1/tdet
      sigtau1=sigtau1/tdet
      tau2=tau2/tdet
      sig2=sig2/tdet
      sigtau2=sigtau2/tdet
      np0=np0/tdet
      np1=np1/tdet
      pp=pp/tdet
      nn=nn/tdet
      tens=tens/tdet
      tenstau=tenstau/tdet
      w%tni2pia=ttni2pia/tdet
      w%tni2piaxd=ttni2piaxd/tdet
      w%tni2piadd=ttni2piadd/tdet
      w%tni2pic=ttni2pic/tdet
      w%tni2picxd=ttni2picxd/tdet
      w%tni2picdd=ttni2picdd/tdet
      w%tni2pitm=ttni2pitm/tdet
      w%tnivd1=ttnivd1/tdet
      w%tnivd2=ttnivd2/tdet
      w%tnive=ttnive/tdet
      w%tni2piapr=ttni2piapr/tdet
      w%tni2piaxdpr=ttni2piaxdpr/tdet
      w%tni2piaddpr=ttni2piaddpr/tdet
      w%v8allpr(1)=vc
      w%vcoul=vcoul/tdet+w%vcoulc
      ccvspr=ccvspr/tdet
      w%v8allpr(2:6)=ccvspr(2:6)
   endif
   w%psi=tdet*exp(-u)*(1.0_r8+u3c)
   totdet=totdet*exp(-u)*(1.0_r8+u3c)
   phig=1.0_r8
   do i=1,npart
      phig=phig*sum(dconjg(smat(:,i))*smat(:,i))
   enddo
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

   subroutine vpotex(r,v)
   real(kind=r8) :: r,v
   if (vexid.eq.0) then
      return
   else if (vexid.eq.1) then
         v=wsv0/(1.0_r8+exp((r-wsrad)/wss))
   else if (vexid.eq.2) then
         v=0.5_r8*hbari*(hbarom*r)**2
   endif
   end subroutine vpotex

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
      if (ic.eq.3) tz=-(psip-psim)/(ci*2.0_r8*dphi*psi)
   enddo
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
      if (ic.eq.3) jz=-(psip-psim)/(ci*2.0_r8*dphi*psi)
   enddo
   j2=aj
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
!   open(10293,file='chkpotold.dat',position='append') !DELETE
!   open(102932,file='hpsiold.dat',position='append') !DELETE
!   write(10293,'(6e15.7)') vc+vtau+vsig+vsigtau !DELETE
!   write(102932,'(6e15.7)') sum(w%v8all(1:6)) !DELETE
!   close(10293) !DELETE
!   close(102932) !DELETE
   end subroutine chkpot

   subroutine chkls(w)
   use stack
   use v6pot
   type (walker) :: w
   complex(kind=r8) :: spsave(4,npart),psi,dpsi(3,npart)
   complex(kind=r8) :: spx(4,15,npart),cvls,cvlstau,p(3,3)
   real(kind=r8) :: dx(3),r,vvls,vvlstau
   integer(kind=i4) :: i,j,is
!  call vnpsi2(w,.false.)
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
   real(kind=r8) :: params(:)
   integer(kind=i4) :: i,j,k
   k=1
   do i=1,nconf
      do j=1,norb(i)
         cdet(k)=params(i)*cdet0(k)
         k=k+1
      enddo
   enddo
   end subroutine setdetparam

   subroutine getdetparam(nparam,params)
   integer(kind=i4) :: nparam
   real(kind=r8), pointer :: params(:)
   nparam=nconf
   allocate(params(nparam))
   params=vorb
   end subroutine getdetparam

   subroutine getderpsi(w,dpsi)
   use stack
   type (walker) :: w
   real(kind=r8) :: dpsi(:)
   integer(kind=i4) :: i,j,k
   dpsi=0.0_r8
   k=1
   do i=1,nconf
      do j=1,norb(i)
         dpsi(i)=dpsi(i)+totdet(k)
         k=k+1
      enddo
   enddo
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
   do i=1,size(dpsi)-22 ! the last parameters are q1c,q2c,q1p,q2p,rscal,a3c,a3a,a3atm,a3avd1,a3avd2,a3avdc3,a3avec3,a3ve,a3comm,rscpp,rscnn,bpp,bnn
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
      do i=5,18 ! rsctni(3),ptni(11)
         call f3p(i)
         call hpsi(w,.false.)
         call f3m(i)
         dpsi(ipar)=(w%psi-psi0)/dp
         ipar=ipar+1
      enddo
   else
      ipar=ipar+14
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
