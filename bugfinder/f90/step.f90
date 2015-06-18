module step
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, save ::  hbar,dt,sigma,etrial,weightm,dts
   integer(kind=i4), private, save :: idmc,npart
   real(kind=r8), private, save :: nfac,el,eli
   complex(kind=r8), private, parameter :: czero = (0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci = (0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone = (1.0_r8,0.0_r8)
   integer(kind=i4), private, save :: wtisto(500),wtistotot(500)
   logical, private, save :: pls,vsym,noprot,pcoul,eucl
   integer(kind=i4), private, save :: nobs
contains
   subroutine setstep(idmcin,npartin,hbarin,dtin,etrialin,wtcut,elin,dtsin,plsin,nobsin,noprotin,nprot,nneut,pcoulin)
   use mympi
   use gofr
   integer(kind=i4) :: idmcin,npartin,nobsin,nprot,nneut
   real(kind=r8) :: hbarin,dtin,etrialin,wtcut,elin,dtsin
   logical :: plsin,noprotin,pcoulin
   idmc=idmcin
   npart=npartin
   hbar=hbarin
   dt=dtin
   dts=dtsin
   sigma=sqrt(2.0_r8*hbar*dt)
   etrial=etrialin
   weightm=wtcut
   wtisto=0
   wtistotot=0
   pcoul=pcoulin
   if (elin.gt.0.0_r8) then
      nfac=1.0_r8/npart
      el=elin
      eli=1.0_r8/el
      pcoul=.false.
   else
      nfac=1.0_r8
      el=0.0_r8
      eli=0.0_r8
   endif
   pls=plsin
   noprot=noprotin
   if (noprot) pcoul=.false.
   vsym=.true. ! symmetrize V6 in the propagator
   if (myrank().eq.0.and.idmc.gt.0) write (6,'(''Symmetrize V6 in the propagator ='',t40,l10)') vsym
   if (myrank().eq.0) write (6,'(''Coulomb included in the energy (and propagation) ='',t50,l10)') pcoul
   call setupgofr(npart,abs(el),nprot,nneut)
   nobs=nobsin
 eucl=.true.
 eucl=.false.
   return
   end subroutine setstep

   subroutine setidmc(idmcin)
   integer(kind=i4) :: idmcin
   idmc=idmcin
   return
   end subroutine setidmc

   subroutine step1(istin,istout,tstep,topt)
   use stack
   use random
   use wavefunction
   use estimator
   use mympi
   integer(kind=i4) :: istin,istout
   real(kind=r8), dimension (3,npart) :: gauss
   real(kind=r8) :: prob,rn(1),wt
   logical :: empty
   real(kind=r8) :: wt1,wt2,wt3,wt4,wtt,wttt
!  real(kind=r8) :: xgauss(15*npart)
!  real(kind=r8) :: xgauss(33*npart)
   real(kind=r8) :: xgauss(35*npart)
   real :: tim1,tim0,tstep,topt
   real(kind=r8) :: wtfacc1,wtfacc2,wtfacc3,wtfacc4
   real(kind=r8) :: wtfacls1,wtfacls2,wtfacls3,wtfacls4
   real(kind=r8) :: phase
   tstep=0.0d0
   topt=0.0d0
   call cpu_time(tim0)
   do while (.true.)
      if (amaboss()) call pop(istin,w1,empty)
      call share(empty)
      if (empty) exit
      call share(w1)
      call setrn(w1%irn)
      gauss=sigma*reshape(gaussian(3*npart),(/3,npart/))
      select case (idmc)
         case(-1,0)
            w2%x=w1%x+gauss
            w2%sp(:,:)=w1%sp(:,:)
            call hpsi(w2,.false.)
            w2%weight=w1%weight
            prob=(abs(w2%psi)/abs(w1%psi))**2
            prob=min(1.0_r8,prob)
            rn=randn(1)
            if (rn(1).lt.prob) then
               call addval(nobs+2,1.0_r8,1.0_r8)
               w1=w2
            else
               w2=w1
               call addval(nobs+2,0.0_r8,1.0_r8)
            endif
            call hspropvmc(w1,w2)
            call hpsi(w2,.false.)
            prob=(abs(w2%psi)/abs(w1%psi))**2
            prob=min(1.0_r8,prob)
            rn=randn(1)
            if (rn(1).lt.prob) then
               call addval(nobs+3,1.0_r8,1.0_r8)
            else
               w2=w1
               call addval(nobs+3,0.0_r8,1.0_r8)
            endif
         case (1:4)
            ! not existing anymore...
         case (5) ! move and rotate without drift, use psig for sampling the paths
            xgauss=gaussian(35*npart)
            wtmp1%x=w1%x+gauss
            call hsprop(w1,wtmp1,xgauss,wtfacc1,wtfacls1,.true.,.true.)
            wtmp2%x=w1%x+gauss
            call hsprop(w1,wtmp2,-xgauss,wtfacc2,wtfacls2,.false.,.false.)
            wtmp3%x=w1%x-gauss
            call hsprop(w1,wtmp3,xgauss,wtfacc3,wtfacls3,.false.,.true.)
            wtmp4%x=w1%x-gauss
            call hsprop(w1,wtmp4,-xgauss,wtfacc4,wtfacls4,.false.,.false.)
            wt1=exp(-wtfacc1*dt+wtfacls1)*wtmp1%psig/w1%psig
            wt2=exp(-wtfacc2*dt+wtfacls2)*wtmp2%psig/w1%psig
            wt3=exp(-wtfacc3*dt+wtfacls3)*wtmp3%psig/w1%psig
            wt4=exp(-wtfacc4*dt+wtfacls4)*wtmp4%psig/w1%psig
            wtt=wt1+wt2+wt3+wt4
            if (real(wtmp1%psi)/real(w1%psi).le.0.0_r8) then
               call addval(nobs+1,1.0_r8,1.0_r8)
            else
               call addval(nobs+1,0.0_r8,1.0_r8)
            endif
            if (real(wtmp2%psi)/real(w1%psi).le.0.0_r8) then
               call addval(nobs+1,1.0_r8,1.0_r8)
            else
               call addval(nobs+1,0.0_r8,1.0_r8)
            endif
            if (real(wtmp3%psi)/real(w1%psi).le.0.0_r8) then
               call addval(nobs+1,1.0_r8,1.0_r8)
            else
               call addval(nobs+1,0.0_r8,1.0_r8)
            endif
            if (real(wtmp4%psi)/real(w1%psi).le.0.0_r8) then
               call addval(nobs+1,1.0_r8,1.0_r8)
            else
               call addval(nobs+1,0.0_r8,1.0_r8)
            endif
            wttt=wt1+wt2+wt3+wt4
            wt1=wt1/wttt
            wt2=wt2/wttt
            wt3=wt3/wttt
            wt4=wt4/wttt
            wt2=wt2+wt1
            wt3=wt3+wt2
            wt4=wt4+wt3
            rn=randn(1)
            if (rn(1).lt.wt1) then
               w2=wtmp1
            elseif (rn(1).lt.wt2) then
               w2=wtmp2
            elseif (rn(1).lt.wt3) then
               w2=wtmp3
            else
               w2=wtmp4
            endif
            wt=0.25_r8*wtt*exp(etrial*dt)
            wt=min(wt,weightm)
            if (real(w2%psi)/real(w1%psi).le.0.0_r8) then ! cp 1
               call addval(nobs+2,1.0_r8,1.0_r8)
            else
               call addval(nobs+2,0.0_r8,1.0_r8)
            endif
            w2%weight=wt
         case (6) ! move and rotate without drift, new weight term like Shiwei's method
            if (.not.pcoul) then
               xgauss=gaussian(33*npart)
            else
               xgauss=gaussian(35*npart)
            endif
! let's assume that we have one boss and three workers here
            if (amaboss()) then
               wtmp1%x=w1%x+gauss
               call hsprop(w1,wtmp1,xgauss,wtfacc1,wtfacls1,.true.,.true.)
            endif           
            if (amaworker().and.idworker().eq.1.or.nworker().eq.0) then
               wtmp2%x=w1%x+gauss
               if (nworker().eq.0) then
                  call hsprop(w1,wtmp2,-xgauss,wtfacc2,wtfacls2,.false.,.false.)
               else
                  call hsprop(w1,wtmp2,-xgauss,wtfacc2,wtfacls2,.true.,.true.)
               endif
            endif
            if (amaworker().and.idworker().eq.2.or.nworker().eq.0) then
               wtmp3%x=w1%x-gauss
               if (nworker().eq.0) then
                  call hsprop(w1,wtmp3,xgauss,wtfacc3,wtfacls3,.false.,.true.)
               else
                  call hsprop(w1,wtmp3,xgauss,wtfacc3,wtfacls3,.true.,.true.)
               endif
            endif
            if (amaworker().and.idworker().eq.3.or.nworker().eq.0) then
               wtmp4%x=w1%x-gauss
               if (nworker().eq.0) then
                  call hsprop(w1,wtmp4,-xgauss,wtfacc4,wtfacls4,.false.,.false.)
               else
                  call hsprop(w1,wtmp4,-xgauss,wtfacc4,wtfacls4,.true.,.true.)
               endif
            endif
            if (nworker().ne.0) then
               call getw(wtmp2,1)
               call getw(wtmp3,2)
               call getw(wtmp4,3)
            endif
            if (amaboss()) then
               wt1=exp(-wtfacc1*dt+wtfacls1)*real(wtmp1%psi/w1%psi)
               wt2=exp(-wtfacc2*dt+wtfacls2)*real(wtmp2%psi/w1%psi)
               wt3=exp(-wtfacc3*dt+wtfacls3)*real(wtmp3%psi/w1%psi)
               wt4=exp(-wtfacc4*dt+wtfacls4)*real(wtmp4%psi/w1%psi)
               wtt=wt1+wt2+wt3+wt4
               if (wt1.le.0.0_r8) then
                  call addval(nobs+1,1.0_r8,1.0_r8)
                  wt1=0.0_r8
               else
                  call addval(nobs+1,0.0_r8,1.0_r8)
               endif
               if (wt2.le.0.0_r8) then
                  call addval(nobs+1,1.0_r8,1.0_r8)
                  wt2=0.0_r8
               else
                  call addval(nobs+1,0.0_r8,1.0_r8)
               endif
               if (wt3.le.0.0_r8) then
                  call addval(nobs+1,1.0_r8,1.0_r8)
                  wt3=0.0_r8
               else
                  call addval(nobs+1,0.0_r8,1.0_r8)
               endif
               if (wt4.le.0.0_r8) then
                  call addval(nobs+1,1.0_r8,1.0_r8)
                  wt4=0.0_r8
               else
                  call addval(nobs+1,0.0_r8,1.0_r8)
               endif
               if (wtt.le.0.0_r8) then
                  wtt=0.0_r8
                  call addval(nobs+2,1.0_r8,1.0_r8)
               else
                  call addval(nobs+2,0.0_r8,1.0_r8)
               endif
               wttt=wt1+wt2+wt3+wt4
               wt1=wt1/wttt
               wt2=wt2/wttt
               wt3=wt3/wttt
               wt4=wt4/wttt
               wt2=wt2+wt1
               wt3=wt3+wt2
               wt4=wt4+wt3
               rn=randn(1)
               if (rn(1).lt.wt1) then
                  w2=wtmp1
               elseif (rn(1).lt.wt2) then
                  w2=wtmp2
               elseif (rn(1).lt.wt3) then
                  w2=wtmp3
               else
                  w2=wtmp4
               endif
               wt=0.25_r8*wtt*exp(etrial*dt)
               wt=min(wt,weightm)
               if (real(w2%psi/w1%psi).le.0.0_r8) wt=0.0_r8
            endif
            call share(wt)
            w2%weight=wt
         case default
            write (6,'(''Illegal idmc value'')')
            stop
      end select
      w2%wt0=w1%wt0
      if (amaboss()) then
         select case (idmc)
            case (1:6)
               phase=atan2(aimag(w2%psi),real(w2%psi))
               call addval(nobs+3,phase,real(w2%weight))
               call addval(nobs+4,real(w2%psi),real(w2%weight))
               call addval(nobs+5,aimag(w2%psi),real(w2%weight))
               call addval(nobs+6,sign(1.0_r8,real(w2%psi)),1.0_r8)
               call addval(nobs+7,sign(1.0_r8,real(w2%psi*w2%wt0)),1.0_r8)
               call addval(nobs+8,real(w2%psi/w2%psig),1.0_r8)
               call addval(nobs+9,aimag(w2%psi/w2%psig),1.0_r8)
               call addval(nobs+14,((1.0_r8-real(w2%weight))/dt+etrial)*nfac,1.0_r8)
               call addval(nobs+16,real(w2%weight),1.0_r8)
         end select
      endif
      call savern(w2%irn)
      call spinnorm(w2)
      if (amaboss()) then
         if (idmc.le.0) then
            call push(istout,w2)
         else
            call branchone(istout,w2)
         endif
      endif
   enddo
   istout=istin
   istin=3-istin
   call cpu_time(tim1)
   tstep=tim1-tim0
   return
   end subroutine step1
 
   subroutine branch(istin,istout,factor)
   use stack
   use random2
   integer(kind=i4) :: istin,istout
   real(kind=r8) :: wt,rn,dummy,factor
   integer(kind=i4) :: i,iwt
   logical :: empty
   do while (.true.) 
      call pop(istin,w1,empty)
      if (empty) then
         istout=istin
         istin=3-istin
         return
      endif
      wt=abs(w1%weight)*factor
      w1%weight=w1%weight/abs(w1%weight)
      call ran1(rn,w1%irn)
      iwt=wt+rn
      do i=1,iwt
         call ran2(dummy,w1%irn)
         call push(istout,w1)
      enddo
      if (iwt.ge.0.and.iwt.lt.500) wtisto(iwt+1)=wtisto(iwt+1)+1
   enddo
   return
   end subroutine branch

   subroutine branchone(istout,w)
   use stack
   use random2
   real(kind=r8) :: wt,rn,dummy
   integer(kind=i4) :: i,iwt,istout
   type (walker) :: w 
   wt=abs(w%weight)
   w%weight=w%weight/abs(w%weight)
   call ran1(rn,w%irn)
   iwt=wt+rn
   do i=1,iwt
      call ran2(dummy,w%irn)
      call push(istout,w)
   enddo
   if (iwt.ge.0.and.iwt.lt.500) wtisto(iwt+1)=wtisto(iwt+1)+1
   return
   end subroutine branchone

   subroutine writewtisto
   use mympi
   integer(kind=i4) :: wttot(500),i
   if (idmc.lt.1) return
   call addall(wtisto,wttot)
   if (myrank().eq.0) then
      wtistotot=wtistotot+wttot
      write (6,*) ''
      write (6,*) 'walkers multiplicity, block, total'
      do i=1,500
         if (i.le.10) write (6,'(i4,2i15)') i-1,wttot(i),wtistotot(i)
         if (i.gt.10.and.wtistotot(i).ne.0) write (6,'(i4,2i15)') i-1,wttot(i),wtistotot(i)
      enddo
      write (6,*) ''
      write (6,*) ''
   endif
   wtisto=0
   end subroutine writewtisto

   subroutine spinnorm(w)
!
! normalize spinors and keep wave function consistent
! step always uses the ratio, so changing the norm this way is OK.
!
   use stack !to define walker type
   type (walker) w
   integer(kind=i4) :: i,j,k
   real(kind=r8) :: anorm
   do i=1,npart
      anorm=1.0_r8/sqrt(sum(w%sp(:,i)*conjg(w%sp(:,i))))
      w%sp(:,i)=anorm*w%sp(:,i)
      w%psi=w%psi*anorm
      w%psig=w%psig*anorm
      if (eucl) then
         do j=1,size(w%osp,1)
            do k=1,size(w%osp,3)
               w%osp(j,:,k,i)=anorm*w%osp(j,:,k,i)
            enddo
         enddo
      endif
   enddo
   end subroutine spinnorm
    
   subroutine hsprop(wold,wnew,xfields,wtfacc,wtfacls,dorold,dornew)
   use stack !to define walker type
   use v6pot
   use matrixmod
   use v3bpot
   use wavefunction
   type (walker) wold,wnew
   real(kind=r8) :: dummy
   real(kind=r8) :: xfields(:)
   real(kind=r8), dimension(npart,npart) :: v2,v3,v4
   real(kind=r8), dimension(3,npart,3,npart) :: v5,v6
   real(kind=r8), dimension(3,npart,3,npart) :: a3st,a3sttm,a3stvd1,a3stvd2,a3stvdc3,a3stvec3
   real(kind=r8), dimension(3,npart,3,npart) :: xpi,xd
   real(kind=r8), allocatable, save :: acsbo(:,:),acsbn(:,:)
   real(kind=r8), dimension(npart,npart) :: atau
   integer(kind=i4) :: xidx,ic,i,j
   logical :: dorold,dornew
   real(kind=r8) :: wtfacls
   real(kind=r8) :: dtfac,a2psc,a2pxdsc,a2pddsc
   real(kind=r8) :: wtfacc
   real(kind=r8), allocatable, save :: vvecso(:,:,:),vvalso(:)
   real(kind=r8), allocatable, save :: vvecsto(:,:,:),vvalsto(:)
   real(kind=r8), allocatable, save :: vvecto(:,:),vvalto(:)
   real(kind=r8), allocatable, save :: vveccsbo(:,:),vvalcsbo(:)
   real(kind=r8), allocatable, save :: vvecsn(:,:,:),vvalsn(:)
   real(kind=r8), allocatable, save :: vvecstn(:,:,:),vvalstn(:)
   real(kind=r8), allocatable, save :: vvectn(:,:),vvaltn(:)
   real(kind=r8), allocatable, save :: vveccsbn(:,:),vvalcsbn(:)
   real(kind=r8) :: rscal(3) ! match the dimension as in jastrowtabop
   if (.not.allocated(acsbo)) allocate(acsbo(npart,npart))
   if (.not.allocated(acsbn)) allocate(acsbn(npart,npart))
   if (.not.allocated(vvecso)) allocate(vvecso(3,npart,3*npart))
   if (.not.allocated(vvalso)) allocate(vvalso(3*npart))
   if (.not.allocated(vvecsto)) allocate(vvecsto(3,npart,3*npart))
   if (.not.allocated(vvalsto)) allocate(vvalsto(3*npart))
   if (.not.allocated(vvecto)) allocate(vvecto(npart,npart))
   if (.not.allocated(vvalto)) allocate(vvalto(npart))
   if (.not.allocated(vvecsn)) allocate(vvecsn(3,npart,3*npart))
   if (.not.allocated(vvalsn)) allocate(vvalsn(3*npart))
   if (.not.allocated(vvecstn)) allocate(vvecstn(3,npart,3*npart))
   if (.not.allocated(vvalstn)) allocate(vvalstn(3*npart))
   if (.not.allocated(vvectn)) allocate(vvectn(npart,npart))
   if (.not.allocated(vvaltn)) allocate(vvaltn(npart))
   if (.not.allocated(vveccsbo)) allocate(vveccsbo(npart,npart))
   if (.not.allocated(vvalcsbo)) allocate(vvalcsbo(npart))
   if (.not.allocated(vveccsbn)) allocate(vveccsbn(npart,npart))
   if (.not.allocated(vvalcsbn)) allocate(vvalcsbn(npart))
   rscal=1.0_r8
   dtfac=1.0_r8
   if (vsym) dtfac=0.5_r8
   wtfacls=0.0_r8
   xidx=1
! apply exp(-V6(R)*dt)
   if (dorold) then
      call hspot(wold%x,dummy,v2,v3,v4,v5,v6,.true.,1)
      call hstnimat(wold%x,dummy,a3st,a3sttm,a3stvd1,a3stvd2,&
           a3stvdc3,a3stvec3,atau,a2psc,a2pxdsc,a2pddsc,xpi,xd,rscal)
      v2=v2+atau
      v6=v6+a2psc*a3st+a3sttm+a3stvd1+a3stvd2+a2pxdsc*a3stvdc3+a2pddsc*a3stvec3
      if (noprot) then
         v2=0.0_r8
         v3=v3+v4
         v4=0.0_r8
         v5=v5+v6
         v6=0.0_r8
      endif
      do ic=1,3
         v5(ic,:,ic,:)=v5(ic,:,ic,:)+v3
         v6(ic,:,ic,:)=v6(ic,:,ic,:)+v4
      enddo
      vvecso=reshape(v5(:,:,:,:),shape(vvecso))
      call eigenrs(vvecso,vvalso,3*npart)
      if (.not.noprot) then
         vvecsto=reshape(v6(:,:,:,:),shape(vvecsto))
         vvecto=v2
         call eigenrs(vvecsto,vvalsto,3*npart)
         call eigenrs(vvecto,vvalto,npart)
      else
         vvecsto=0.0_r8
         vvecto=0.0_r8
         vvalsto=0.0_r8
         vvalto=0.0_r8
      endif
      if (pcoul) then
         call getvcsb(wold%x,acsbo)
         acsbo=acsbo/4.0_r8
         vveccsbo=acsbo
         call eigenrs(vveccsbo,vvalcsbo,npart)
      else
         vveccsbo=0.0_r8
         vvalcsbo=0.0_r8
      endif
   endif
   call propgv6(wold%sp,wnew%sp,dt*dtfac,xfields(xidx:xidx+15*npart), &
        vvecto,vvalto,vvecso,vvalso,vvecsto,vvalsto)
   if (eucl) then
      do i=1,size(wold%osp,1)
         do j=1,size(wold%osp,3)
            call propgv6(wold%osp(i,:,j,:),wnew%osp(i,:,j,:),dt*dtfac,xfields(xidx:xidx+15*npart), &
               vvecto,vvalto,vvecso,vvalso,vvecsto,vvalsto)
         enddo
      enddo
   endif
   xidx=xidx+15*npart
   if (pcoul) then
      call propgcoul(wnew%sp,wnew%sp,dt*dtfac,xfields(xidx:xidx+npart),vveccsbo,vvalcsbo,acsbo)
      if (eucl) then
         do i=1,size(wold%osp,1)
            do j=1,size(wold%osp,3)
               call propgcoul(wnew%osp(i,:,j,:),wnew%osp(i,:,j,:),dt*dtfac,xfields(xidx:xidx+npart),vveccsbo,vvalcsbo,acsbo)
            enddo
         enddo
      endif
      xidx=xidx+npart
   endif
   if (pls) then
      call propgls(wnew%sp,wnew%sp,wold%x,wnew%x,dt,xfields(xidx:xidx+3*npart),wtfacls)
      if (eucl) then
         do i=1,size(wold%osp,1)
            do j=1,size(wold%osp,3)
               call propgls(wnew%osp(i,:,j,:),wnew%osp(i,:,j,:),wold%x,wnew%x,dt,xfields(xidx:xidx+3*npart),wtfacls)
            enddo
         enddo
      endif
      xidx=xidx+3*npart
   endif
   if (vsym) then
! apply exp(-V6(R')*dt)
      if (dornew) then
         call hspot(wnew%x,dummy,v2,v3,v4,v5,v6,.true.,1)
         call hstnimat(wnew%x,dummy,a3st,a3sttm,a3stvd1,a3stvd2,&
              a3stvdc3,a3stvec3,atau,a2psc,a2pxdsc,a2pddsc,xpi,xd,rscal)
         v2=v2+atau
         v6=v6+a2psc*a3st+a3sttm+a3stvd1+a3stvd2+a2pxdsc*a3stvdc3+a2pddsc*a3stvec3
         if (noprot) then
            v2=0.0_r8
            v3=v3+v4
            v4=0.0_r8
            v5=v5+v6
            v6=0.0_r8
         endif
         do ic=1,3
            v5(ic,:,ic,:)=v5(ic,:,ic,:)+v3
            v6(ic,:,ic,:)=v6(ic,:,ic,:)+v4
         enddo
         vvecsn=reshape(v5(:,:,:,:),shape(vvecsn))
         call eigenrs(vvecsn,vvalsn,3*npart)
         if (.not.noprot) then
            vvecstn=reshape(v6(:,:,:,:),shape(vvecstn))
            vvectn=v2
            call eigenrs(vvecstn,vvalstn,3*npart)
            call eigenrs(vvectn,vvaltn,npart)
         else
            vvecstn=0.0_r8
            vvectn=0.0_r8
            vvalstn=0.0_r8
            vvaltn=0.0_r8
         endif
         if (pcoul) then
            call getvcsb(wnew%x,acsbn)
            acsbn=acsbn/4.0_r8
            vveccsbn=acsbn
            call eigenrs(vveccsbn,vvalcsbn,npart)
         else
            vveccsbn=0.0_r8
            vvalcsbn=0.0_r8
         endif
      endif
      if (pcoul) then
         call propgcoul(wnew%sp,wnew%sp,dt*dtfac,xfields(xidx:xidx+npart),vveccsbn,vvalcsbn,acsbn)
         if (eucl) then
            do i=1,size(wold%osp,1)
               do j=1,size(wold%osp,3)
                  call propgcoul(wnew%osp(i,:,j,:),wnew%osp(i,:,j,:),dt*dtfac,xfields(xidx:xidx+npart),vveccsbn,vvalcsbn,acsbn)
               enddo
            enddo
         endif
         xidx=xidx+npart
      endif
      call propgv6(wnew%sp,wnew%sp,dt*dtfac,xfields(xidx:xidx+15*npart), &
           vvectn,vvaltn,vvecsn,vvalsn,vvecstn,vvalstn)
      if (eucl) then
         do i=1,size(wold%osp,1)
            do j=1,size(wold%osp,3)
               call propgv6(wnew%osp(i,:,j,:),wnew%osp(i,:,j,:),dt*dtfac,xfields(xidx:xidx+15*npart), &
                  vvectn,vvaltn,vvecsn,vvalsn,vvecstn,vvalstn)
            enddo
         enddo
      endif
      xidx=xidx+15*npart
   endif
   call hpsi(wnew,.false.)
   if (vsym) then
      if (pcoul) then
         wtfacc=0.5_r8*(wnew%vc+wold%vc+wnew%tnic+wold%tnic+wnew%vcoulc+wold%vcoulc+wnew%vext+wold%vext)
      else
         wtfacc=0.5_r8*(wnew%vc+wold%vc+wnew%tnic+wold%tnic+wnew%vext+wold%vext)
      endif
   else
      if (pcoul) then
         wtfacc=wnew%vc+wnew%tnic+wnew%vcoulc+wnew%vext
      else
         wtfacc=wnew%vc+wnew%tnic+wnew%vext
      endif
   endif
   return
   end subroutine hsprop

   subroutine hspropvmc(wold,wnew)
   use stack !to define walker type
   use random
   type (walker) wold,wnew
   real(kind=r8) :: gauss(3)
   complex(kind=r8) :: rotation(3,npart),sptmp(4,npart)
   complex(kind=r8) :: arg,c,s
   integer(kind=i4) :: i,ic
! spin rotation
   do i=1,npart ! loop over eigenvectors
      gauss=gaussian(3)
      do ic=1,3
         rotation(ic,i)=dts*gauss(ic)
      enddo
   enddo
   do i=1,npart ! loop over particles rotating spins
      arg=sqrt(sum(rotation(:,i)**2))
      if (arg.ne.czero) rotation(:,i)=rotation(:,i)/arg
      c=cos(arg)
      s=sin(arg)
! Proton and Neutron rotatate the same way
      sptmp(1,i)=(c+s*rotation(3,i))*wold%sp(1,i) &
         +s*(rotation(1,i)-ci*rotation(2,i))*wold%sp(2,i)
      sptmp(3,i)=(c+s*rotation(3,i))*wold%sp(3,i) &
         +s*(rotation(1,i)-ci*rotation(2,i))*wold%sp(4,i)
      sptmp(2,i)=(c-s*rotation(3,i))*wold%sp(2,i) &
         +s*(rotation(1,i)+ci*rotation(2,i))*wold%sp(1,i)
      sptmp(4,i)=(c-s*rotation(3,i))*wold%sp(4,i) &
         +s*(rotation(1,i)+ci*rotation(2,i))*wold%sp(3,i)
   enddo
!isospin rotation
   do i=1,npart ! loop over eigenvectors
      gauss=gaussian(3)
      do ic=1,3
         rotation(ic,i)=dts*gauss(ic)
      enddo
   enddo
   do i=1,npart ! loop over particles rotating isospin
      arg=sqrt(sum(rotation(:,i)**2))
      if (arg.ne.czero) rotation(:,i)=rotation(:,i)/arg
      c=cos(arg)
      s=sin(arg)
! Spin up and down rotatate the same way
      wnew%sp(1,i)=(c+s*rotation(3,i))*sptmp(1,i) &
         +s*(rotation(1,i)-ci*rotation(2,i))*sptmp(3,i)
      wnew%sp(2,i)=(c+s*rotation(3,i))*sptmp(2,i) &
         +s*(rotation(1,i)-ci*rotation(2,i))*sptmp(4,i)
      wnew%sp(3,i)=(c-s*rotation(3,i))*sptmp(3,i) &
         +s*(rotation(1,i)+ci*rotation(2,i))*sptmp(1,i)
      wnew%sp(4,i)=(c-s*rotation(3,i))*sptmp(4,i) &
         +s*(rotation(1,i)+ci*rotation(2,i))*sptmp(2,i)
   enddo
   return
   end subroutine hspropvmc

   subroutine rotate(spold,spnew,rott,rots,rotst)
   use matrixmod
   complex(kind=r8) :: spold(4),spnew(4)
   complex(kind=r8) :: stmat(4,4)
   complex(kind=r8) :: rots(3),rotst(3,3),rott(3)
   stmat=0.0_r8
   stmat(1,1)=rots(3)+rott(3)+rotst(3,3)
   stmat(1,2)=rots(1)-ci*rots(2)+rotst(1,3)-ci*rotst(2,3)
   stmat(1,3)=rott(1)-ci*rott(2)+rotst(3,1)-ci*rotst(3,2)
   stmat(1,4)=rotst(1,1)-ci*rotst(2,1)-ci*rotst(1,2)-rotst(2,2)
   stmat(2,1)=rots(1)+ci*rots(2)+rotst(1,3)+ci*rotst(2,3)
   stmat(2,2)=-rots(3)+rott(3)-rotst(3,3)
   stmat(2,3)=rotst(1,1)-ci*rotst(1,2)+ci*rotst(2,1)+rotst(2,2)
   stmat(2,4)=rott(1)-ci*rott(2)-rotst(3,1)+ci*rotst(3,2)
   stmat(3,1)=rott(1)+ci*rott(2)+rotst(3,1)+ci*rotst(3,2)
   stmat(3,2)=rotst(1,1)+ci*rotst(1,2)-ci*rotst(2,1)+rotst(2,2)
   stmat(3,3)=rots(3)-rott(3)-rotst(3,3)
   stmat(3,4)=rots(1)-ci*rots(2)-rotst(1,3)+ci*rotst(2,3)
   stmat(4,1)=rotst(1,1)+ci*rotst(2,1)+ci*rotst(1,2)-rotst(2,2)
   stmat(4,2)=rott(1)+ci*rott(2)-rotst(3,1)-ci*rotst(3,2)
   stmat(4,3)=rots(1)+ci*rots(2)-rotst(1,3)-ci*rotst(2,3)
   stmat(4,4)=-rots(3)-rott(3)+rotst(3,3)
   spnew(:)=expmult(stmat(:,:),spold(:),4)
   end subroutine rotate

   subroutine rotatesp(spold,spnew,rots)
   use matrixmod
   complex(kind=r8) :: spold(4),spnew(4)
   complex(kind=r8) :: stmat(4,4)
   complex(kind=r8) :: rots(3)
   stmat=0.0_r8
   stmat(1,1)=rots(3)
   stmat(1,2)=rots(1)-ci*rots(2)
   stmat(2,1)=rots(1)+ci*rots(2)
   stmat(2,2)=-rots(3)
   stmat(3,3)=rots(3)
   stmat(3,4)=rots(1)-ci*rots(2)
   stmat(4,3)=rots(1)+ci*rots(2)
   stmat(4,4)=-rots(3)
   spnew(:)=expmult(stmat(:,:),spold(:),4)
   end subroutine rotatesp

   subroutine rotatetauz(spold,spnew,rottz)
   use matrixmod
   complex(kind=r8) :: spold(4),spnew(4)
   complex(kind=r8) :: stmat(4,4)
   complex(kind=r8) :: rottz
   stmat=0.0_r8
   stmat(1,1)=rottz
   stmat(2,2)=rottz
   stmat(3,3)=-rottz
   stmat(4,4)=-rottz
   spnew(:)=expmult(stmat(:,:),spold(:),4)
   end subroutine rotatetauz

   subroutine propgv6(spold,spnew,dtt,xfields,vvect,vvalt,vvecs,vvals,vvecst,vvalst)
   complex(kind=r8) :: spold(:,:),spnew(:,:)
   complex(kind=r8) :: cfac
   real(kind=r8) :: vvect(:,:),vvalt(:)
   real(kind=r8) :: vvecs(:,:,:),vvals(:)
   real(kind=r8) :: vvecst(:,:,:),vvalst(:)
   real(kind=r8) :: xfields(:),dtt
   integer(kind=i4) :: i,is,xidx
   complex(kind=r8) :: rots(3,npart),rotst(3,3,npart),rott(3,npart)
   rots=0.0_r8
   rott=0.0_r8
   rotst=0.0_r8
   xidx=1
   do i=1,3*npart ! loop over eigenvectors
      cfac=sqrt(-cmplx(vvals(i)*dtt,0.0_r8))
      rots=rots+xfields(xidx)*cfac*vvecs(:,:,i)
      xidx=xidx+1
   enddo
   do is=1,3 ! loop over taux, tauy, tauz
      do i=1,3*npart ! loop over eigenvectors
         cfac=sqrt(-cmplx(vvalst(i)*dtt,0.0_r8))
         rotst(:,is,:)=rotst(:,is,:)+xfields(xidx)*cfac*vvecst(:,:,i)
         xidx=xidx+1
      enddo
      do i=1,npart ! loop over eigenvectors
         cfac=sqrt(-cmplx(vvalt(i)*dtt,0.0_r8))
         rott(is,:)=rott(is,:)+xfields(xidx)*cfac*vvect(:,i)
         xidx=xidx+1
      enddo
   enddo
   do i=1,npart ! loop over particles rotating spins
      call rotate(spold(:,i),spnew(:,i),rott(:,i),rots(:,i),rotst(:,:,i))
   enddo
   end subroutine propgv6

   subroutine propgcoul(spold,spnew,dtt,xfields,vveccsb,vvalcsb,acsb)
   complex(kind=r8) :: spold(:,:),spnew(:,:)
   complex(kind=r8) :: cfac
   real(kind=r8) :: vveccsb(:,:),vvalcsb(:),acsb(:,:)
   real(kind=r8) :: xfields(:),dtt
   integer(kind=i4) :: i,j,xidx
   complex(kind=r8) :: rottz(npart)
   rottz=0.0_r8
   xidx=1
   do i=1,npart
      cfac=sqrt(-cmplx(vvalcsb(i)*dtt,0.0_r8))
      rottz(:)=rottz(:)+xfields(xidx)*cfac*vveccsb(:,i)
      xidx=xidx+1
   enddo
   do i=1,npart ! loop over particles rotating spins
      call rotatetauz(spold(:,i),spnew(:,i),rottz(i))
      do j=1,npart
         spnew(1,i)=exp(-acsb(i,j)*dtt)*spnew(1,i)
         spnew(2,i)=exp(-acsb(i,j)*dtt)*spnew(2,i)
         spnew(3,i)=exp(acsb(i,j)*dtt)*spnew(3,i)
         spnew(4,i)=exp(acsb(i,j)*dtt)*spnew(4,i)
      enddo
   enddo
   end subroutine propgcoul

   subroutine propgls(spold,spnew,xold,xnew,dt,xfields,wtfacls)
   use v6pot
   use matrixmod
   complex(kind=r8) :: spold(:,:),spnew(:,:)
   real(kind=r8) :: xold(:,:),xnew(:,:)
   complex(kind=r8) :: cfac
   real(kind=r8) :: xfields(:),dt
   integer(kind=i4) :: i,j,ic,jc,xidx
   real(kind=r8) :: dxold(3),dxnew(3),dx(3),cross(3),gls(3,npart)
   real(kind=r8) :: rnew,vvls,wtfacls,rold,vvls1
   complex(kind=r8) :: rots(3,npart)
   real(kind=r8) :: dummy
   real(kind=r8) :: eigvec(3,npart,3*npart),eigval(3*npart)
   rots=czero
   gls=0.0_r8
   wtfacls=0.0_r8
   eigvec=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         dxold=xold(:,i)-xold(:,j)
         dxold=dxold-el*nint(dxold*eli)
         dxnew=xnew(:,i)-xnew(:,j)
         dxnew=dxnew-el*nint(dxnew*eli)
         rnew=sqrt(sum(dxnew**2))
         call vls(rnew,vvls,dummy,1)
         rold=sqrt(sum(dxold**2))
         call vls(rold,vvls1,dummy,1)
         vvls=0.5_r8*(vvls+vvls1)
         dx=dxnew-dxold
         dx=dx-el*nint(dx*eli)
         cross(1)=0.5_r8*(dxold(2)+dxnew(2))*dx(3)-0.5_r8*(dxold(3)+dxnew(3))*dx(2)
         cross(2)=0.5_r8*(dxold(3)+dxnew(3))*dx(1)-0.5_r8*(dxold(1)+dxnew(1))*dx(3)
         cross(3)=0.5_r8*(dxold(1)+dxnew(1))*dx(2)-0.5_r8*(dxold(2)+dxnew(2))*dx(1)
         gls(:,i)=gls(:,i)+cross*vvls
         gls(:,j)=gls(:,j)+cross*vvls
      enddo
   enddo
   gls=-1.0_r8/(8.0_r8*hbar)*gls
! add LS counterterms (eq. 3.53 in the thesis)
   eigvec=0.0_r8
   do i=2,npart
      do j=1,i-1
         do ic=1,3
            do jc=1,3
               eigvec(ic,i,3*(j-1)+jc)=eigvec(ic,i,3*(j-1)+jc)-gls(jc,j)*gls(ic,i)/dt
               eigvec(jc,j,3*(i-1)+ic)=eigvec(jc,j,3*(i-1)+ic)-gls(jc,j)*gls(ic,i)/dt
            enddo
         enddo
      enddo
   enddo
   call eigenrs(eigvec,eigval,3*npart)
   rots=ci*gls
   wtfacls=0.5_r8*sum(gls*gls)
   xidx=1
   do i=1,3*npart ! loop over eigenvectors
      cfac=sqrt(-cmplx(eigval(i)*dt,0.0_r8))
      rots=rots+xfields(xidx)*cfac*eigvec(:,:,i)
      xidx=xidx+1
   enddo
   do i=1,npart ! loop over particles rotating spins
      call rotatesp(spold(:,i),spnew(:,i),rots(:,i))
   enddo
   end subroutine propgls
end module step
