module wavefunction
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   integer(kind=i4), private, save :: npart,ndet
   complex(kind=r8), private, save, allocatable :: cdet(:),cdet0(:)
   real(kind=r8), private, save :: a
   integer(kind=i4), private, save :: nconf
   integer(kind=i4), private, save, allocatable :: norb(:)
   real(kind=r8), private, save, allocatable :: vorb(:)
   integer (kind=i4), private, parameter :: levi(2,3)= &
      reshape((/2,3, 3,1, 1,2/),(/2,3/))
contains
   subroutine setpsi(npartin,elin,hbar)
   use mympi
   use phimod
   integer(kind=i4) :: npartin
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
   logical :: rdiv
   npart=npartin
   if (myrank().eq.0) then
      read (5,*) notab   !table points for orbitals
      read (5,*) romax   !rmax for orbitals
      read (5,*) rdiv    !divido orbitali per r?
      read (5,*) npart   !particle number
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
   endif
   call bcast(npart)
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
   end subroutine setpsi



   subroutine hpsi(w,doall)
   use stack ! to define walker type
   use v6pot
   use jastrow
   use matrixmod
   use phimod
   type (walker) :: w
   complex(kind=r8) :: spx(4,15,npart),sxall(npart,15,npart)
   complex(kind=r8) :: ph(npart,4,npart,ndet),dph(npart,3,4,npart,ndet)
   complex(kind=r8) :: d2ph(npart,4,npart,ndet),d2det
   complex(kind=r8) :: sxmall(npart,15,npart),sall(npart,3,npart)
   complex(kind=r8) :: small(npart,3,npart)
   complex(kind=r8) :: det,ddet(3,npart),d2cor
   complex(kind=r8) :: tdet,tddet(3,npart),td2det
   complex(kind=r8) :: smati(npart,npart),dsum(3)
   complex(kind=r8) :: d1,d2,d3,d4,cvs,cvt,cvst
   complex(kind=r8) :: ccvs,ccvt,ccvst
   complex(kind=r8) :: sigma(3,npart),tau(3,npart),sigmatau(3,3,npart)
   complex(kind=r8) :: dpall(npart,3,0:15,npart),dphall(npart,3,0:15,npart)
   complex(kind=r8) :: cvls,cvlst,c1,c2,c3,cdot,c1t,c2t,c3t,cdott
   real(kind=r8) :: r,dx(3),uj,ujp,ujpp,eni
   real(kind=r8) :: asig(3,npart,3,npart),asigtau(3,npart,3,npart)
   real(kind=r8) :: atau(npart,npart),du(3,npart),u,d2u
   real(kind=r8) :: vvls,vvlst
   integer(kind=i4) :: i,j,ic1,ic2,i1,i2,jc1,jc2,it
   integer(kind=i4) :: ic,jc,kc,idet,i3
   logical :: doall
   complex(kind=r8) :: smat(npart,npart),phig
   integer(kind=i4), parameter :: levi(2,3)=reshape((/2,3, 3,1, 1,2/),(/2,3/))
!
! hpsi calculates < Psi_T | R S > or < Psi_T | O | R S >.
! < Psi_T | R S > = Psi_T^*(R,S) so orbitals from getphi with L_z = 1
! should have exp(-i phi) angular dependence.
!
! zero center of mass
!
   eni=1.0_r8/npart
   w%x(1,:)=w%x(1,:)-sum(w%x(1,:))*eni
   w%x(2,:)=w%x(2,:)-sum(w%x(2,:))*eni
   w%x(3,:)=w%x(3,:)-sum(w%x(3,:))*eni
   spx=opmult(w%sp) ! multiply by sigma, tau, and sigma tau
   call hspot(w%x,w%vc,asig,asigtau,atau,doall) ! calculate potential
   u=0.0_r8
   d2u=0.0_r8
   du=0.0_r8
   do i=2,npart
      do j=1,i-1
         dx=w%x(:,i)-w%x(:,j)
         r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         call jastro(r,uj,ujp,ujpp)
         u=u+uj
         d2u=d2u+ujpp
         du(:,i)=du(:,i)-ujp*dx
         du(:,j)=du(:,j)+ujp*dx
      enddo
   enddo
   d2u=2.0_r8*d2u
!
! calculate orbitals -- this is ugly with the i in the middle. The compiler
! will have to copy these to make this work. Think about fixing this.
!
   do i=1,npart
      call getphi(w%x(:,i),ph(:,:,i,:),dph(:,:,:,i,:),d2ph(:,:,i,:))
   enddo
!
! calculate determinants
!
   tdet=czero
   tddet=czero
   td2det=czero
   w%sigma=czero
   w%tau=czero
   w%sigmatau=czero
   w%v=czero
   w%vls=czero
   w%vlst=czero
   ccvs=czero
   ccvt=czero
   ccvst=czero
   do idet=1,ndet
      do i=1,npart
         smati(:,i)=matmul(ph(:,:,i,idet),w%sp(:,i))
      enddo
      smat=smati
      call matinv(smati,det,npart)
      det=cdet(idet)*det
      tdet=tdet+det
      d2det=czero
      ddet=czero
      if (doall) then
         do i=1,npart
            ddet(1,i)=sum(smati(i,:)*matmul(dph(:,1,:,i,idet),w%sp(:,i)))
            ddet(2,i)=sum(smati(i,:)*matmul(dph(:,2,:,i,idet),w%sp(:,i)))
            ddet(3,i)=sum(smati(i,:)*matmul(dph(:,3,:,i,idet),w%sp(:,i)))
            d2det=d2det+sum(smati(i,:)*matmul(d2ph(:,:,i,idet),w%sp(:,i)))
         enddo
         dsum(1)=sum(ddet(1,:))*eni
         dsum(2)=sum(ddet(2,:))*eni
         dsum(3)=sum(ddet(3,:))*eni
         ddet(1,:)=ddet(1,:)-dsum(1)
         ddet(2,:)=ddet(2,:)-dsum(2)
         ddet(3,:)=ddet(3,:)-dsum(3)
!
! two orbital second derivative terms for center of mass correction
!
         do i=1,npart
            sall(:,:,i)=reshape(matmul( &
               reshape(dph(:,:,:,i,idet),(/3*npart,4/)),w%sp(:,i)),(/npart,3/))
         enddo
         small=reshape(matmul(smati,reshape(sall,(/npart,3*npart/))),shape(small))
         d2cor=czero
         do ic1=1,3
            do i1=1,npart-1
               d1=small(i1,ic1,i1)
               do i2=i1+1,npart
                  d2=small(i2,ic1,i2)
                  d3=small(i1,ic1,i2)
                  d4=small(i2,ic1,i1)
                  d2cor=d2cor+d1*d2-d3*d4
               enddo
            enddo
         enddo
         d2det=(cone-eni)*d2det-2.0_r8*eni*d2cor
         tddet=tddet+det*ddet
         td2det=td2det+det*d2det
!
! expectation values
!
         do i=1,npart
            do ic=1,3
               sigma(ic,i)=sum(smati(i,:)*matmul(ph(:,:,i,idet),spx(:,ic,i)))
               tau(ic,i)=sum(smati(i,:)*matmul(ph(:,:,i,idet),spx(:,ic+3,i)))
               do jc=1,3
                  kc=jc+3*(ic-1)+6
                  sigmatau(ic,jc,i)=sum(smati(i,:)*matmul(ph(:,:,i,idet),spx(:,kc,i)))
               enddo
            enddo
         enddo
         w%sigma=w%sigma+det*sigma
         w%tau=w%tau+det*tau
         w%sigmatau=w%sigmatau+det*sigmatau
!
! potential terms
!
         do i=1,npart
            sxall(:,:,i)=matmul(ph(:,:,i,idet),spx(:,:,i))
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
! spin-orbit terms
         do i=1,npart
            dpall(:,:,0,i)=reshape(matmul(    &
               reshape(dph(:,:,:,i,idet),(/3*npart,4/)),w%sp(:,i)),(/npart,3/))
            dpall(:,:,1:15,i)=reshape(matmul( &
               reshape(dph(:,:,:,i,idet),(/3*npart,4/)),spx(:,:,i)),(/npart,3,15/))
         enddo
         dphall=reshape(matmul(smati,reshape(dpall,(/npart,48*npart/))),shape(dphall))
         cvls=czero
         cvlst=czero
         do i=2,npart
            do j=1,i-1
               dx=w%x(:,i)-w%x(:,j)
               r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
               call vls(r,vvls,vvlst)
               c1=czero
               c2=czero
               c3=czero
               cdot=czero
               c1t=czero
               c2t=czero
               c3t=czero
               cdott=czero
               do i1=1,3
                  i2=levi(1,i1)
                  i3=levi(2,i1)
                  c1=c1+dx(i1)*( &
                    dphall(i,i2,0,i)*sxmall(j,i3,j) &
                   -dphall(j,i2,0,i)*sxmall(i,i3,j) &
                   -dphall(i,i3,0,i)*sxmall(j,i2,j) &
                   +dphall(j,i3,0,i)*sxmall(i,i2,j))
                  c2=c2+dx(i1)*( &
                   -dphall(j,i2,0,j)*sxmall(i,i3,i) &
                   +dphall(i,i2,0,j)*sxmall(j,i3,i) &
                   +dphall(j,i3,0,j)*sxmall(i,i2,i) &
                   -dphall(i,i3,0,j)*sxmall(j,i2,i))
                  cdot=cdot+(sxmall(i,i1,i)+sxmall(j,i1,j))*( &
                     dx(i2)*(du(i3,i)-du(i3,j))-dx(i3)*(du(i2,i)-du(i2,j)))
                  c3=c3+dx(i1)*(dphall(i,i2,i3,i)-dphall(i,i3,i2,i) &
                    -dphall(j,i2,i3,j)+dphall(j,i3,i2,j))
                  do it=1,3
                     c1t=c1t+dx(i1)*( &
                        dphall(i,i2,3+it,i)*sxmall(j,3*i3+3+it,j) &
                       -dphall(j,i2,3+it,i)*sxmall(i,3*i3+3+it,j) &
                       -dphall(i,i3,3+it,i)*sxmall(j,3*i2+3+it,j) &
                       +dphall(j,i3,3+it,i)*sxmall(i,3*i2+3+it,j))
                     c2t=c2t+dx(i1)*( &
                       -dphall(j,i2,3+it,j)*sxmall(i,3*i3+3+it,i) &
                       +dphall(i,i2,3+it,j)*sxmall(j,3*i3+3+it,i) &
                       +dphall(j,i3,3+it,j)*sxmall(i,3*i2+3+it,i) &
                       -dphall(i,i3,3+it,j)*sxmall(j,3*i2+3+it,i))
                     cdott=cdott+( &
                         sxmall(i,3*i1+3+it,i)*sxmall(j,3+it,j) &
                        -sxmall(j,3*i1+3+it,i)*sxmall(i,3+it,j) &
                        +sxmall(j,3*i1+3+it,j)*sxmall(i,3+it,i) &
                        -sxmall(i,3*i1+3+it,j)*sxmall(j,3+it,i))*( &
                        dx(i2)*(du(i3,i)-du(i3,j))-dx(i3)*(du(i2,i)-du(i2,j)))
                     c3t=c3t+dx(i1)*( &
                        dphall(i,i2,3*i3+3+it,i)*sxmall(j,3+it,j) &
                       -dphall(j,i2,3*i3+3+it,i)*sxmall(i,3+it,j) &
                       -dphall(i,i3,3*i2+3+it,i)*sxmall(j,3+it,j) &
                       +dphall(j,i3,3*i2+3+it,i)*sxmall(i,3+it,j) &
                       -dphall(j,i2,3*i3+3+it,j)*sxmall(i,3+it,i) &
                       +dphall(i,i2,3*i3+3+it,j)*sxmall(j,3+it,i) &
                       +dphall(j,i3,3*i2+3+it,j)*sxmall(i,3+it,i) &
                       -dphall(i,i3,3*i2+3+it,j)*sxmall(j,3+it,i))
                  enddo
               enddo
               cvls=cvls+vvls*(c1+c2+c3+cdot)
               cvlst=cvlst+vvlst*(c1t+c2t+c3t+cdott)
            enddo
         enddo
         cvls=cvls*ci*0.25_r8
         cvlst=cvlst*ci*0.25_r8
         w%v=w%v+det*(cvs+cvt+cvst+cvls+cvlst)
         w%vls=w%vls+det*cvls
         w%vlst=w%vlst+det*cvlst
         ccvs=ccvs+det*cvs
         ccvt=ccvt+det*cvt
         ccvst=ccvst+det*cvst
      endif
   enddo
   if (doall) then
      w%sigma=w%sigma/tdet
      w%tau=w%tau/tdet
      w%sigmatau=w%sigmatau/tdet
      w%v=w%v/tdet+w%vc
      w%vls=w%vls/tdet
      w%vlst=w%vlst/tdet
      tddet=tddet/tdet
      td2det=td2det/tdet
      w%dpsi=du+tddet
      w%d2psi=-d2u+td2det+sum((du+2.0_r8*tddet)*du)
   endif
   w%psi=tdet*exp(-u)
   phig=1.0_r8
   do i=1,npart
      phig=phig*sum(dconjg(smat(:,i))*smat(:,i))
   enddo
   w%psig=sqrt(dconjg(tdet)*tdet+a*phig)*exp(-u)
   ccvs=ccvs/tdet
   ccvt=ccvt/tdet
   ccvst=ccvst/tdet
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

   subroutine chkls(w)
   use stack
   use v6pot
   type (walker) :: w
   complex(kind=r8) :: spsave(4,npart),psi,dpsi(3,npart)
   complex(kind=r8) :: spx(4,15,npart),cvls,cvlst,p(3,3),pt(3,3)
   real(kind=r8) :: dx(3),r,vvls,vvlst
   integer(kind=i4) :: i,j,is,it,kc
   call hpsi(w,.true.)
   write(6,'(''hpsi,  vls= '',6e15.7)') w%vls
   write(6,'(''hpsi,  vlst= '',6e15.7)') w%vlst
   spsave=w%sp
   psi=w%psi
   dpsi=w%dpsi
   spx=opmult(w%sp)
   cvls=czero
   cvlst=czero
   do i=2,npart
      do j=1,i-1
         dx=w%x(:,i)-w%x(:,j)
         r=sqrt(sum(dx**2))
         call vls(r,vvls,vvlst)
         pt=0.0_r8
         do is=1,3
            w%sp(:,i)=spx(:,is,i)
            call hpsi(w,.true.)
            p(:,is)=(w%dpsi(:,i)-w%dpsi(:,j))*w%psi
            w%sp(:,i)=spsave(:,i)
            w%sp(:,j)=spx(:,is,j)
            call hpsi(w,.true.)
            p(:,is)=p(:,is)+(w%dpsi(:,i)-w%dpsi(:,j))*w%psi
            w%sp(:,j)=spsave(:,j)
            do it=1,3
               kc=it+3*(is-1)+6
               w%sp(:,i)=spx(:,kc,i)
               w%sp(:,j)=spx(:,it+3,j)
               call hpsi(w,.true.)
               pt(:,is)=pt(:,is)+(w%dpsi(:,i)-w%dpsi(:,j))*w%psi
               w%sp(:,i)=spx(:,it+3,i)
               w%sp(:,j)=spx(:,kc,j)
               call hpsi(w,.true.)
               pt(:,is)=pt(:,is)+(w%dpsi(:,i)-w%dpsi(:,j))*w%psi
               w%sp(:,i)=spsave(:,i)
               w%sp(:,j)=spsave(:,j)
            enddo
         enddo
         cvls=cvls+vvls*(dx(1)*(p(2,3)-p(3,2))+dx(2)*(p(3,1)-p(1,3))+dx(3)*(p(1,2)-p(2,1)))
         cvlst=cvlst+vvlst*(dx(1)*(pt(2,3)-pt(3,2))+dx(2)*(pt(3,1)-pt(1,3))+dx(3)*(pt(1,2)-pt(2,1)))
      enddo
   enddo
   cvls=0.25_r8*ci*cvls/psi
   cvlst=0.25_r8*ci*cvlst/psi
   write(6,'(''chkls, vls= '',6e15.7)') cvls
   write(6,'(''chkls, vlst= '',6e15.7)') cvlst
! cvlstau to be done!!!!!!!!!!!!!
   end subroutine chkls




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the following are subrotuines needed by the optimization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The subroutines are fake because in this wave function
! there are no parameters to optimize!

   subroutine setdetparam(params)
   real(kind=r8) :: params(:)
   end subroutine setdetparam

   subroutine getdetparam(nparam,params)
   integer(kind=i4) :: nparam
   real(kind=r8), pointer :: params(:)
   end subroutine getdetparam

   subroutine getderpsi(w,dpsi)
   use stack
   type (walker) :: w
   real(kind=r8) :: dpsi(:)
   end subroutine getderpsi

   subroutine getjasder(w,dpsi)
   use stack
   type (walker) :: w
   real(kind=r8) :: dpsi(:)
   end subroutine getjasder

end module wavefunction
