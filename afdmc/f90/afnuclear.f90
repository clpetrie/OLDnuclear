   program afnuclear
   use stack
   use lattice
   use random
   use random2
   use estimator
   use mympi
   use v6pot
   use jastrow
   use step
   use wavefunction
   use optimizer
   use ioconfs
   use gofr
   use v3bpot
   use correlator
!  use euclidean
   use propv3
   use operators
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer, parameter :: i8=selected_int_kind(15)
   integer, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i8) :: irn,irn0
   integer(kind=i4) :: npart,nwalk,idmc,neq,nav,nstep,i,j,k,it,is,nsteptr
   integer(kind=i4) :: istin,istout,ntab,n1,ito,nt,lpot,lpotpr
   integer(kind=i4) :: ibox
   real(kind=r8) :: hbar,dt,dts,etrial,popup,popdn,stmax
   real(kind=r8) :: el,rn(1),xtemp(3),dcut
   real(kind=r8) :: vfact(8),dummy
   real(kind=r8) :: wtcut
   integer(kind=i4), allocatable :: numfig(:)
   real :: time0,time1,time2,timebl,timebl0,timeop,ts,to,tlb
   integer :: s_time0, s_time1, clock_rate
   real(kind=r8) :: timest,timeopc,timelb,tts,tto,ttb
   real(kind=r8) :: t0,t0t,tbl,tblt,top,topt,tot,tott
   character(len=90), dimension(:), allocatable :: answer
   character(len=20) :: potnam(32) = (/ &
      'malfliet-tjon       ', &
      'reid v8             ', &
      'urbana v14          ', &
      'argonne v8 (cut v14)',&
      'argonne v9          ', &
      'argonne v14         ', &
      'argonne v18-csbl    ', &
      'argonne v18-csbs    ', &
      'argonne v18         ', &
      'argonne v8''         ', &
      'argonne v6''         ', &
      'argonne v4''         ', &
      'argonne v2''         ', &
      'argonne v1''         ', &
      'argonne vx''         ', &
      'argonne v18p        ', &
      'argonne v18pq       ', &
      'super-soft core v14 ', &
      'super-soft core v8'' ', &
      'paris               ', &
      'argonne v18 v1.9    ', &
      'argonne v8'' v1.9    ', &
      'argonne v18 v1.7    ', &
      'argonne v8'' v1.7    ', &
      'argonne v7'' - Joe   ', &
      'Minnesota           ', &
      'LO (R0=1.0 fm)      ', &
      'NLO (R0=1.0 fm)     ', &
      'N2LO (R0=1.0 fm)    ', &
      'LO (R0=1.2 fm)      ', &
      'NLO (R0=1.2 fm)     ', &
      'N2LO (R0=1.2 fm)    '/)
   logical :: isite,empty,iout,operat
   real(kind=r8) :: ene
   integer(kind=i4) :: optiter,iters,ifile,nest,ifile2,nobs,dov3,iblk
   real(kind=r8) :: dobal
   real(kind=r8) :: tau,tau0
   real(kind=r8) :: nfac,fac
   complex(kind=r8) :: psi,psit,psi2,psi2t,psig,psigt,psig2,psig2t
   logical :: ivmc,iop,iwrt
   complex(kind=r8) :: j2,jz,t2,tz,j2t,jzt,t2t,tzt
   character(13) :: fileconfs
   logical :: pls,dovls,filexist,noprot,pcoul
   real(kind=r8) :: enow,eg,egerr,egp
   complex(kind=r8) :: v3,ek,vn
   integer(kind=i4) :: nw0,nw1
   integer(kind=i4), parameter :: nworkers=0 ! number of workers for each boss
   character(len=3) :: tni
   real(kind=r8) :: a2psc,a2pxdsc,a2pddsc
   integer(kind=i4) :: nprot,nneut
   integer(kind=i4) :: neucop,nspop
   character(len=10) :: filext
   logical :: fcsb,eucl
   pls=.true.! include LS in the propagator
   dovls=.true.! calculate LS energy
   pcoul=.false. ! include Coulomb in the propagator
!  eucl=.true.
   eucl=.false.
   call init0(nworkers) ! mpi initialization that needs to be done before reading
   if (myrank().eq.0) then
      read (5,*) irn     !random number seed for sites
      irn0=irn
      read (5,*) hbar    !hbar^2/2m
      read (5,*) idmc    !0=vmc
      read (5,*) isite   !if true take initial walkers from lattice sites
      read (5,*) nwalk   !default number of walkers
      read (5,*) neq     !equilibration blocks
! if nav is negative, then write out configurations after each block
      read (5,*) nav     !averaging blocks
      iwrt=.false.
      if (nav.lt.0) then
         nav=-nav
         iwrt=.true.
      endif
      read (5,*) nstep   !steps per block
! don't compute averages if nstep is negative, but write out configurations
      iop=.true.
      if (nstep.lt.0) then
         nstep=-nstep
         iop=.false.
         iwrt=.true.
      endif
      read (5,*) dt      !time step
      read (5,*) dts     !time step spin
      read (5,*) etrial  !trial energy
      read (5,*) dcut    !drift cutoff
      read (5,*) popup   !kill walkers if more than nwalk*popup
      read (5,*) popdn   !create walkers if less than nwalk*popdn
      read (5,*) dobal   !rebalance if walker number fluctuate more than this (per cpu)
      read (5,*) wtcut   !cutoff for weight
      read (5,*) ntab    !number of table points for jastrow and potential
      read (5,*) lpot    !Bob's potential number
      if (lpot.lt.0) then
         lpot=abs(lpot)
         pcoul=.true.
      endif
      read (5,*) lpotpr  !Bob's potential number
      do i=1,8
         read (5,*) vfact(i) !multiply potential
      enddo
      read (5,*) tni
      read (5,*) a2psc   !scale A2pi by this in the propagation
      read (5,*) a2pxdsc   !scale A2pixd by this in the propagation
      read (5,*) a2pddsc   !scale A2pixd by this in the propagation
      read (5,*) ibox    !lattice sum for potential
   endif
   call bcast(hbar)
   call bcast(idmc)
   call bcast(isite)
   call bcast(nwalk)
   call bcast(neq)
   call bcast(nav)
   call bcast(nstep)
   call bcast(iwrt)
   call bcast(iop)
   call bcast(dt)
   call bcast(dts)
   call bcast(etrial)
   call bcast(dcut)
   call bcast(popup)
   call bcast(popdn)
   call bcast(dobal)
   call bcast(wtcut)
   call bcast(ntab)
   call bcast(lpot)
   call bcast(pcoul)
   call bcast(lpotpr)
   call bcast(vfact)
   call bcast(tni)
   call bcast(a2psc)
   call bcast(a2pxdsc)
   call bcast(a2pddsc)
   call bcast(ibox)
   if (myrank().eq.0) then
      iout=.true.
      if (idmc.eq.-1) iout=.false.
      select case (idmc)
         case (-100)
            write (6,'(''Calculate observables from existing configurations'')')
         case (-1)
            write (6,'(''Variational Monte Carlo Run and optimization'')')
         case (0)
            write (6,'(''Variational Monte Carlo Run'')')
         case (5)
            write (6,'(''Unconstrained Diffusion Monte Carlo Run -- no drift'')')
         case (6)
            write (6,'(''Diffusion Monte Carlo Run -- no drift'')')
         case (7)
            write (6,'(''Hybrid Unconstrained Diffusion Monte Carlo Run'')')
         case default
            write (6,'(''Illegal idmc value'')')
            call abort
      end select
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''MPI parameters:'')')
      write (6,'(''Number of boss cpus ='',t40,i10)') nboss()
      write (6,'(''Number of worker cpus for each boss ='',t40,i10)') nworkers
      write (6,'(''Total number of cpus ='',t40,i10)') nproc()
      write (6,'(''Simulation parameters:'')')
      write (6,'(''random number seed ='',t30,i20)') irn
      write (6,'(''hbar^2/2m ='',t40,f10.5)') hbar
      write (6,'(''start from sites ='',t40,l10)') isite
      write (6,'(''nwalk ='',t40,i10)') nwalk
      if (idmc.ne.7) then
         write (6,'(''equilibration blocks ='',t40,i10)') neq
      else
         write (6,'(''number of unconstrained steps ='',t40,i10)') neq
      endif
      write (6,'(''averaging blocks ='',t40,i10)') nav
      write (6,'(''nsteps/blocks ='',t40,i10)') nstep
      write (6,'(''time step ='',t35,f15.10)') dt
      write (6,'(''time step rotations (vmc) ='',t35,f15.10)') dts
      write (6,'(''trial energy ='',t40,f10.5)') etrial
      write (6,'(''drift cutoff ='',t40,f10.5)') dcut
      write (6,'(''drastic population control up ='',t40,f10.5)') popup
      write (6,'(''drastic population control down ='',t40,f10.5)') popdn
      write (6,'(''rebalance if walkers fluctuate more than '',t45,f6.2,'' %'')') dobal
      write (6,'(''weight cut off ='',t40,f10.5)') wtcut
      write (6,'(''table size ='',t40,i10)') ntab
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''System parameters:'')')
      write (6,'(''lpot ='',t40,i10)') lpot
      write (6,'(''lpotpr ='',t40,i10)') lpotpr
      if (abs(lpot).le.26) then 
         write (6,'(''potential  ='',t40,a20)') potnam(abs(lpot))
      else
         if (lpot.eq.112) write (6,'(''potential  ='',t40,a20)') potnam(27)
         if (lpot.eq.113) write (6,'(''potential  ='',t40,a20)') potnam(28)
         if (lpot.eq.114) write (6,'(''potential  ='',t40,a20)') potnam(29)
         if (lpot.eq.122) write (6,'(''potential  ='',t40,a20)') potnam(30)
         if (lpot.eq.123) write (6,'(''potential  ='',t40,a20)') potnam(31)
         if (lpot.eq.124) write (6,'(''potential  ='',t40,a20)') potnam(32)
      endif
      if (abs(lpotpr).le.26) then 
         write (6,'(''potentialpr ='',t40,a20)') potnam(abs(lpotpr))
      else
         if (lpotpr.eq.112) write (6,'(''potential  ='',t40,a20)') potnam(27)
         if (lpotpr.eq.113) write (6,'(''potential  ='',t40,a20)') potnam(28)
         if (lpotpr.eq.114) write (6,'(''potential  ='',t40,a20)') potnam(29)
         if (lpotpr.eq.122) write (6,'(''potential  ='',t40,a20)') potnam(30)
         if (lpotpr.eq.123) write (6,'(''potential  ='',t40,a20)') potnam(31)
         if (lpotpr.eq.124) write (6,'(''potential  ='',t40,a20)') potnam(32)
      endif
      write (6,'(''multiply central v by ='',t40,f10.5)') vfact(1)
      write (6,'(''multiply tau v by ='',t40,f10.5)') vfact(2)
      write (6,'(''multiply sigma v by ='',t40,f10.5)') vfact(3)
      write (6,'(''multiply sigma tau v by ='',t40,f10.5)') vfact(4)
      write (6,'(''multiply tensor v by ='',t40,f10.5)') vfact(5)
      write (6,'(''multiply tensor tau v by ='',t40,f10.5)') vfact(6)
      write (6,'(''multiply spin-orbit v by ='',t40,f10.5)') vfact(7)
      write (6,'(''multiply spin-orbit tau v by ='',t40,f10.5)') vfact(8)
      write (6,'(''three-body potential ='',t40,a3)') tni
      write (6,'(''ibox = (0 is nearest image)'',t40,i10)') ibox
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      if (lpot.eq.11) pls=.false.
      if (lpot.eq.11.and.lpotpr.eq.11) dovls=.false.
      if (idmc.gt.1) write (6,'(''propagate LS ='',t40,l10)') pls
      write (6,'(''Calculate LS ='',t40,l10)') dovls
   endif
   call bcast(iout)
   call bcast(pls)
   call bcast(dovls)
   call setpsi(npart,el,hbar,dovls,noprot,nprot,nneut)
   if (myrank().eq.0) then
      write(6,'(''Number of protons ='',t40,i10)') nprot
      write(6,'(''Number of neutrons ='',t40,i10)') nneut
   endif
   etrial=etrial*npart
   if (el.gt.0.0_r8) then
      nfac=1.0_r8/npart
   else
      nfac=1.0_r8
   endif
   call v6potinit(npart,ibox,el,lpot,ntab,vfact,lpotpr)
   call v3bpotinit(npart,el,tni,dov3,a2psc,a2pxdsc,a2pddsc)
   call setdov3(dov3)
   call jasinit(el,ntab,2.0_r8*hbar,abs(lpot),fcsb)
   call setoptv3(dov3.ne.0)
   call setoptcsb(fcsb)
!call initgv3(npart,dt)
!stop
   optiter=1
   if (idmc.eq.-1) then
      call setupopt(npart,optiter)
      operat=.false.  ! do not compute operators during optimization
   endif
   stmax=20.0_r8
   if (nproc().gt.16) then
      stmax=0.2_r8
!     stmax=0.5_r8
   endif
   nobs=54 ! this is the number of common observables in VMC and AFDMC
   call setstep(idmc,npart,hbar,dt,etrial,wtcut,el,dts,pls,nobs,noprot,nprot,nneut,pcoul)
   call setop(idmc,npart,hbar,el,nobs,pcoul)
   select case (idmc) !initialize estimators
      case (-100)
         nest=nobs
         call setestnum(nest)
         if (myrank().eq.0) allocate(answer(nest))
         filext='exconf'
      case (-1:0)
         nest=nobs+3
         call setestnum(nest)
         if (myrank().eq.0) allocate(answer(nest))
         call addest(nobs+1,'jackson-feenberg energy')
         call addest(nobs+2,'acceptance')
         call addest(nobs+3,'acceptance-spin')
         filext='vmc'
      case (1:7)
         nest=nobs+16
         call setestnum(nest)
         if (myrank().eq.0) allocate(answer(nest))
         call addest(nobs+1,'single crossing/step')
         call addest(nobs+2,'global crossing/step')
         call addest(nobs+3,'phase psit')
         call addest(nobs+4,'real psit')
         call addest(nobs+5,'im psit')
         call addest(nobs+6,'sign real(psit)')
         call addest(nobs+7,'sign real(psit*wt0)')
         call addest(nobs+8,'real (psit/psig)')
         call addest(nobs+9,'aimag (psit/psig)')
         call addest(nobs+10,'energy numerator(real)')
         call addest(nobs+11,'denominator psit/psigi(real)')
         call addest(nobs+12,'energy numerator(aimag)')
         call addest(nobs+13,'denominator psit/psigi(aimag)')
         call addest(nobs+14,'growth(lin dt) energy')
         call addest(nobs+15,'growth(global)')
         call addest(nobs+16,'growth energy')
         nest=nest-2 ! the growth energy will need extra care
         filext='dmc'
         if (idmc.eq.5) filext='unc'
   end select
   call addest(1,'total energy-propagated')
   call addest(2,'kinetic energy')
   call addest(3,'potential energy')
   call addest(4,'im energy')
   call addest(5,'imkin energy')
   call addest(6,'impot energy')
   call addest(7,'v cent')
   call addest(8,'v tau')
   call addest(9,'v sigma')
   call addest(10,'v sigtau')
   call addest(11,'v tensor')
   call addest(12,'v tenstau')
   call addest(13,'v ls')
   call addest(14,'vls tau')
   call addest(15,'vpot v6')
   call addest(16,'vpot v8')
   call addest(17,'v coulomb')
   call addest(18,'tni 2pipr-anticommutator')
   call addest(19,'tni 2pixdpr-anticommutator')
   call addest(20,'tni 2piddpr-anticommutator')
   call addest(21,'tni 2pi-commutator')
   call addest(22,'tni 2pixd-commutator')
   call addest(23,'tni 2pidd-commutator')
   call addest(24,'tni 2pi-tm')
   call addest(25,'tni vd1')
   call addest(26,'tni vd2')
   call addest(27,'tni ve')
   call addest(28,'tni central')
   call addest(29,'vpottot v8+tni')
   call addest(30,'enetot v8+tni')
   call addest(31,'vpr cent')
   call addest(32,'vpr tau')
   call addest(33,'vpr sigma')
   call addest(34,'vpr sigtau')
   call addest(35,'vpr tensor')
   call addest(36,'vpr tenstau')
   call addest(37,'vpr ls')
   call addest(38,'vprls tau')
   call addest(39,'vpotpr v6')
   call addest(40,'vpotpr v8')
   call addest(41,'vpotv6-vpotprv6 diff')
   call addest(42,'vpotv8-vpotprv8 diff')
   call addest(43,'tni 2pi-anticommutator')
   call addest(44,'tni 2pixd-anticommutator')
   call addest(45,'tni 2pidd-anticommutator')
   call addest(46,'2piapr-2pia diff')
   call addest(47,'2piapr-2pia-2pic')
   call addest(48,'2piaprxd-2piaxd-2picxd')
   call addest(49,'2piaprdd-2piadd-2picdd')
   call addest(50,'proton rms')
   call addest(51,'neutron rms')
   call addest(52,'total rms')
   call addest(53,'neutron-proton rms')
   call addest(54,'external potential')
   ivmc=idmc.le.0
   call init1(npart,dobal,ivmc) ! additional mpi initialization here
   if (eucl) then
!     call initeuc(neucop,nspop,npart)
   else
      neucop=0
      nspop=0
   endif
   if (idmc.eq.7) then
      call create(3,int(nwalk*stmax),npart,neucop,nspop) !create stacks
   else
      call create(2,int(nwalk*stmax),npart,neucop,nspop) !create stacks
   endif
   allocate(numfig(0:nproc()-1))
   if (myrank().eq.0) call system_clock(s_time0)
   call cpu_time(time0)

   if (idmc.eq.-100) then
! calculate expectation values from configurations, then exit
! use nav to specify which file with name 'confs.xxx' should be the first
      call zerest
      do i=1,100000
         call cpu_time(time0)
         write(fileconfs,*) 100+i+nav
         fileconfs='confs.'//trim(adjustl(fileconfs))
         inquire(file=fileconfs,exist=filexist)
         if (filexist) then
            call readconfs(fileconfs,tau,nfac,hbar,npart)
            call update
            if (myrank().eq.0) then
               answer=resstring(tau)
               write (6,'(a90)') (answer(k),k=1,size(answer))
            endif
            call cpu_time(time1)
            t0=time1-time0
            call addall(t0,t0t)
            if (myrank().eq.0) then
               write (6,'(''Time to process file '',a15,'' = '',f10.2)') fileconfs,t0t/nproc()
               write (6,'(/)')
            endif
         endif
      enddo
      call done
      if (myrank().eq.0) write (6,*) 'Finished!'
      stop
   endif

   do iters=1,optiter
      if (optiter.gt.1.and.myrank().eq.0) write (6,'(''optimization iteration number '',i10)') iters
      if (optiter.gt.1) call updateparams
! careful, I change el to start from some better config for nuclei
!     if (el.lt.0.0_r8) el=-1.0_r8
      if (el.lt.0.0_r8) el=-10.0_r8
      if (iters.eq.1) then
         if (myrank().eq.0) then
            if (isite) then   !get initial walkers
               irn=irn0
               call setrn(irn)    !set seed for sites
               write (6,'(''Initial walkers from sites'')')
               w1%x=(bestcubic(npart)-0.5_r8)*abs(el)+.01_r8*abs(el)*reshape(randn(size(w1%x)),shape(w1%x))
               if (el.gt.0.0_r8) w1%x=(w1%x-el*nint(w1%x/el))
               w1%sp=(0.0_r8,0.0_r8)
               do i=1,npart
                  do is=1,4
                     rn=randn(1)
                     w1%sp(is,i)=rn(1)
                  enddo
               enddo
               n1=nwalk
               tau=0.0_r8
               write (6,'(''Configurations at tau ='',t35,f15.10)') tau
            else
               write (6,'(''Initial walkers from file'')')
               rewind 9
               read (9,'(i10,f15.8)') n1,tau
               if (idmc.le.0) tau=0.0_r8
               write (6,'(''Number of walkers ='',t40,i10)') n1
               write (6,'(''Configurations at tau ='',t35,f15.10)') tau
               if (idmc.eq.5) then
                  tau=0.0_r8
                  write (6,'(''Reset tau ='',t35,f15.10)') tau
               endif
            endif
         endif
         call bcast(n1)
         call bcast(tau)
         tau0=0.0_r8
         do k=1,n1
            if (myrank().eq.0) then
               if (isite) then
                  do i=1,npart
                     rn=randn(1)
                     j=npart*rn(1)+1
                     xtemp(:)=w1%x(:,i)
                     w1%x(:,i)=w1%x(:,j)
                     w1%x(:,j)=xtemp(:)
                  enddo
                  w1%weight=1.0_r8
                  w1%irn=irn
                  call ran2(dummy,w1%irn)
                  irn=w1%irn
               elseif (iters.eq.1) then
                  read (9,'(6e15.7)') w1%x
                  read (9,'(6e15.7)') w1%sp
                  read (9,'(4e15.7)') w1%weight,w1%wt0
                  read (9,'(4e15.7)') w1%psi0,w1%psig0
                  read (9,'(i20)') w1%irn
               endif
               if (iters.eq.1) call push(1,w1)
            endif
            call barrier
            if (iters.eq.1) then
               ito=mod(k-1,nboss())
               call movewalkers(1,1,0,ito,k)
            endif
         enddo
         istin=1
         istout=2
         n1=numstack(istin)
         call cpu_time(time1)
         t0=time1-time0
         call addall(t0,t0t)
         if (myrank().eq.0) then
            write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
            write (6,'(''Time for reading configurations (sec.) = '',f10.2)') t0t/nproc()
            write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
         endif
         psi=0.0_r8
         psi2=0.0_r8
         psig=0.0_r8
         psig2=0.0_r8
         t2t=0.0_r8
         tzt=0.0_r8
         j2t=0.0_r8
         jzt=0.0_r8
         call zerest
         call barrier
         do i=1,n1
            call pop(istin,w1,empty)
            if (empty) exit
            if (tau.eq.0.0_r8.and.iop.and.idmc.ne.5.and..not.isite) then
               call hpsi(w1,.true.)
            else
               call hpsi(w1,.false.)
            endif
!write (6,*) ' calling chkder '
!call chkder(w1,7.e-5_r8,1.e-4_r8)
!write (6,*) ' calling chkpot'
!call chkpot(w1)
!call chkls(w1)
!write (6,*) ' calling chkop'
!call chkop(w1)
!write (6,*) ' calling chkderparams'
!if (idmc.eq.-1) call chkderp(w1,0.0025_r8)
!stop
            psi=psi+w1%psi*dconjg(w1%psi)
            psi2=psi2+(w1%psi*dconjg(w1%psi))**2
            psig=psig+w1%psig**2
            psig2=psig2+w1%psig**4
!           w1%wt0=1.0_r8 ! this is terrible
            w1%wt0=abs(w1%psi)/w1%psi ! this is 1/exp(i*theta)
!           w1%wt0=w1%psig/w1%psi ! this seems better
!           w1%wt0=w1%psig/abs(w1%psi) ! terrible
            w1%weight=1.0_r8
            if (idmc.le.0) w1%wt0=1.0_r8
            ek=-hbar*w1%d2psi
            vn=sum(w1%v8all(:))+w1%tnic+w2%tni2piapr+w1%tni2pitm+w1%tnivd1+w1%tnivd2+w1%tnive+w1%tni2piaxdpr+w1%tni2piaddpr
            call addval(1,real(ek+vn)*nfac,real(w1%weight))
            call addval(2,real(ek)*nfac,real(w1%weight))
            call addval(3,real(vn)*nfac,real(w1%weight))
            call addval(7,real(w1%v8all(1))*nfac,real(w1%weight))
            call addval(8,real(w1%v8all(2))*nfac,real(w1%weight))
            call addval(9,real(w1%v8all(3))*nfac,real(w1%weight))
            call addval(10,real(w1%v8all(4))*nfac,real(w1%weight))
            call addval(11,real(w1%v8all(5))*nfac,real(w1%weight))
            call addval(12,real(w1%v8all(6))*nfac,real(w1%weight))
            call addval(13,real(w1%v8all(7))*nfac,real(w1%weight))
            call addval(14,real(w1%v8all(8))*nfac,real(w1%weight))
            call addval(17,real(w1%vcoul)*nfac,real(w1%weight))
            call addval(18,real(w1%tni2piapr)*nfac,real(w1%weight))
            call addval(19,real(w1%tni2piaxdpr)*nfac,real(w1%weight))
            call addval(20,real(w1%tni2piaddpr)*nfac,real(w1%weight))
            call addval(21,real(w1%tni2pic)*nfac,real(w1%weight))
            call addval(22,real(w1%tni2picxd)*nfac,real(w1%weight))
            call addval(23,real(w1%tni2picdd)*nfac,real(w1%weight))
            call addval(24,real(w1%tni2pitm)*nfac,real(w1%weight))
            call addval(25,real(w1%tnivd1)*nfac,real(w1%weight))
            call addval(26,real(w1%tnivd2)*nfac,real(w1%weight))
            call addval(27,real(w1%tnive)*nfac,real(w1%weight))
            call addval(28,w1%tnic*nfac,real(w1%weight))
            vn=sum(w1%v8all(:))+w1%tnic+w1%tni2piapr+w1%tni2pic+w1%tni2pitm+w1%tnivd1+w1%tnivd2+w1%tnive+w1%tni2piaxdpr+w1%tni2piaddpr
            call addval(30,real(ek+vn)*nfac,real(w1%weight))
            if (el.lt.0.0_r8) then
               call checkj(w1,1.e-2_r8,t2,tz,j2,jz)
               t2t=t2t+t2
               tzt=tzt+tz
               j2t=j2t+j2
               jzt=jzt+jz
            endif
!           if (eucl) call setupeuc(w1)
!call overlap(w1)
            call push(istout,w1)
         enddo
         istout=istin
         istin=3-istin
         call barrier
         call addall(n1,nt)
         call addall(psi,psit)
         call addall(psi2,psi2t)
         call addall(psig,psigt)
         call addall(psig2,psig2t)
         call update   !collect block averages
         if (el.lt.0.0_r8) then
            call addall(t2t,t2)
            call addall(tzt,tz)
            call addall(j2t,j2)
            call addall(jzt,jz)
            t2=t2/nt
            tz=tz/nt
            j2=j2/nt
            jz=jz/nt
         endif
         if (myrank().eq.0) then
            write(6,'(''Initial configurations:'')')
            if (tau.eq.0.0_r8.and.iop.and..not.isite) then
               answer=resstring(tau)
               write (6,'(a90)') answer(1)
               write (6,'(a90)') answer(2)
               write (6,'(a90)') answer(3)
               write (6,'(a90)') answer(7)
               write (6,'(a90)') answer(8)
               write (6,'(a90)') answer(9)
               write (6,'(a90)') answer(10)
               write (6,'(a90)') answer(11)
               write (6,'(a90)') answer(12)
               write (6,'(a90)') answer(13)
               write (6,'(a90)') answer(14)
               write (6,'(a90)') answer(17)
               write (6,'(a90)') answer(18)
               write (6,'(a90)') answer(19)
               write (6,'(a90)') answer(20)
               write (6,'(a90)') answer(21)
               write (6,'(a90)') answer(22)
               write (6,'(a90)') answer(23)
               write (6,'(a90)') answer(24)
               write (6,'(a90)') answer(25)
               write (6,'(a90)') answer(26)
               write (6,'(a90)') answer(27)
               write (6,'(a90)') answer(28)
               write (6,'(a90)') answer(30)
            endif
            write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
            psi=psit/nt
            psi2=psi2t/nt
            psi2=sqrt(abs(psi2-psi*psi)/nt)
            psig=psigt/nt
            psig2=psig2t/nt
            psig2=sqrt(abs(psig2-psig*psig)/nt)
            write(6,'(''Psi**2 of initial configurations  = '',2e15.7,''  +-  '',2e15.7)') psi,psi2
            write(6,'(''Psig**2 of initial configurations = '',2e15.7,''  +-  '',2e15.7)') psig,psig2
            write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
            if (el.lt.0.0_r8) then
               write(6,'(''T**2 of initial configurations  = '',2e15.7)') t2
               write(6,'(''Tz of initial configurations    = '',2e15.7)') tz
               write(6,'(''J**2 of initial configurations  = '',2e15.7)') j2
               write(6,'(''Jz of initial configurations    = '',2e15.7)') jz
               write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
            endif
         endif
      endif
      call zerest
      call cpu_time(time2)
      t0=time2-time1
      call addall(t0,t0t)
      if (myrank().eq.0) then
         write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
         write (6,'(''Time to process initial configurations (sec.) = '',f10.2)') t0t/nproc()
         write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      endif
      call barrier
      it=0
      call zerotimebal
      ifile=0
      ifile2=0
      nw0=nt
      iblk=0
      if (idmc.ne.7) then
         do i=1,nav+neq  !do blocks
            if (i.eq.neq+1) then
               if (myrank().eq.0.and.iout) write (6,'(/,''equilibration done'')')
               call zerest
               it=0
            endif
            it=it+1
            call cpu_time(timebl0)
            timest=0.0_r8
            timeopc=0.0_r8
            timelb=0.0_r8
            do j=1,nstep
               call step1(istin,istout,ts,to)
               timest=timest+ts
               n1=numstack(istin)
               call addall(n1,nt)
               call bcast(nt)
               nw1=nt
               if (amaboss().and.idmc.gt.0) then
                  if (myrank().eq.0) then
                     fac=real(nw1)/real(nw0)
                     call addval(nobs+15,fac,1.0_r8)
                  endif
                  if (nt.gt.popup*nwalk.or.nt.lt.popdn*nwalk) call branch(istin,istout,1.0_r8*nwalk/nt)
                  tau=tau+dt
               endif
               tau0=tau0+dt
               call barrier
               call checkpop(istin,tlb)
               timelb=timelb+tlb
               n1=numstack(istin)
               call addall(n1,nt)
               call bcast(nt)
               nw0=nt  ! update nw0 after the drastic population control
            enddo
            if (amaboss()) then
               if (iop.and.i.gt.neq) then
                  call calcop(istin,istout,to)
! call writeall(istin,istout,tau,1)
               endif
            else
               to=0.0_r8
            endif
            timeopc=timeopc+to
            call cpu_time(timebl)
            n1=numstack(istin)
            call addall(n1,nt)
            call bcast(nt)
            call update   !collect block averages
            if (iout.or.i.eq.nav+neq) then
               if (myrank().eq.0) then
                  if (idmc.ne.-1) then
                     write (6,'(/,''iteration ='',t30,i14)') it
                     write (6,'(''walker number = '',t30,i14)') nt
                  endif
                  if (iop) then
                     if (idmc.gt.1.or.i.gt.neq) then
                        answer=resstring(tau)
                        write (6,'(a90)') (answer(k),k=1,nest)
                     else
                        answer=resstring(tau)
                        write (6,'(a90)') answer(nobs+2)
                        write (6,'(a90)') answer(nobs+3)
                     endif
                  endif
               endif
               if (.not.iop) then ! write out everything but local energy
                  if (idmc.eq.0.and.myrank().eq.0) then
                     answer=resstring(tau)
                     write (6,'(a90)') answer(12)
                     write (6,'(a90)') answer(13)
                  else if (myrank().eq.0) then
                     answer=resstring(tau)
                     write (6,'(a90)') answer(nest)
                  endif
               endif
               if (idmc.ge.1.and.myrank().eq.0) then ! write out growth energies
                  call result(nest+2,enow,egerr)
                  egp=enow+egerr
                  eg=(etrial-log(enow)/dt)*nfac
                  egp=(etrial-log(egp)/dt)*nfac
                  egerr=abs(egp-eg)
                  write (6,'(''growth energy '',t30,1x,1p,e15.8,e19.10,'' +- '',e19.10)') tau,eg,egerr
                  call result(nest+1,enow,egerr)
                  egp=enow+egerr
                  eg=(etrial-log(enow)/dt)*nfac
                  egp=(etrial-log(egp)/dt)*nfac
                  egerr=abs(egp-eg)
                  write (6,'(''growth(global) energy '',t30,1x,1p,e15.8,e19.10,'' +- '',e19.10)') tau,eg,egerr
               endif
               if (iwrt) then ! write out configurations
                  iblk=iblk+1
                  if (mod(iblk,5).eq.0) then ! write out configurations every 20 blocks
                     do 
                        ifile=ifile+1
                        write(fileconfs,*) 100+ifile
                        fileconfs='confs.'//trim(adjustl(fileconfs))
                        inquire(file=fileconfs,exist=filexist)
                        if (.not.filexist) exit
                     enddo
                     call barrier
                     n1=numstack(istin)
                     call gather(n1,numfig)
                     call bcast(numfig)
                     nt=sum(numfig)
                     call writeconfs(fileconfs,istin,istout,nt,tau,numfig,lpot,idmc)
                     call cpu_time(timeop)         
                     t0=timeop-timebl
                     call addall(t0,t0t)
                     if (myrank().eq.0) then
                        write (6,'(''Time to write file '',a15,'' = '',f10.2)') fileconfs,t0t/nproc()
                     endif
                  endif
               endif
! write out other operators
               if (iop.and.i.gt.neq) then
                  call writegofr(tau,filext)
               endif
            endif
            call cpu_time(timeop)         
            tbl=timebl-timebl0
            call addall(tbl,tblt)
            top=timeop-timebl
            call addall(top,topt)
            tot=timeop-time2
            call addall(tot,tott)
            call addall(timest,tts)
            call addall(timeopc,tto)
            call addall(timelb,ttb)
            if (myrank().eq.0.and.iout) then
               write (6,'(''Time for propagation (sec) (min) = '',t40,2f11.3)') tts/nproc(),tts/nproc()/60.0_r8
               write (6,'(''Time for energies (sec) (min) = '',t40,2f11.3)') tto/nproc(),tto/nproc()/60.0_r8
               write (6,'(''Time for load balancing (sec) (min) = '',t40,2f11.3)') ttb/nproc(),ttb/nproc()/60.0_r8
               write (6,'(''Time for block (sec) (min) = '',t40,2f11.3)') tblt/nproc(),tblt/nproc()/60.0_r8
               if (iop) write (6,'(''Time for operators (sec) (min) = '',t40,2f11.3)') topt/nproc(),topt/nproc()/60.0_r8
               write (6,'(''Elapsed block time (sec) (min) = '',t40,2f11.3)') tott/nproc(),tott/nproc()/60.0_r8
            endif
            call writewtisto
            if (idmc.eq.-1.and.i.gt.neq) then
               do while(.true.)
                  call pop(istin,w1,empty)
                  if (empty) exit
                  call hpsi(w1,.true.)
                  v3=w1%tnic+w1%tni2piapr+w1%tni2piaxd+w1%tni2piadd+w1%tni2pic+w1%tni2picxd+w1%tni2picdd &
                    +w1%tni2pitm+w1%tnivd1+w1%tnivd2+w1%tnive
                  ene=(-hbar*real(w1%d2psi)+real(sum(w1%v8all(:))+w1%vcoul+w1%vext+v3))*nfac
                  call calder(w1,ene)
                  call push(istout,w1)
               enddo
               istout=istin
               istin=3-istout
            endif
         enddo
      else ! idmc=7
 nsteptr=neq
         it=0
         do i=1,nav
            it=it+1
            call setidmc(6)
            tau0=0.0_r8
            call cpu_time(timebl0)
            do j=1,nstep ! do some equilibration
               call step1(istin,istout,ts,to)
               n1=numstack(istin)
               call addall(n1,nt)
               call bcast(nt)
               if (nt.gt.popup*nwalk.or.nt.lt.popdn*nwalk) call branch(istin,istout,1.0_r8*nwalk/nt)
               call checkpop(istin,tlb)
            enddo
            call calcop(istin,istout,to)
            n1=numstack(istin)
            call addall(n1,nt)
            call bcast(nt)
            call update   !collect block averages
            if (myrank().eq.0) then
               write (6,'(/,''iteration ='',t30,i14)') it
               write (6,'(''walker number = '',t30,i14)') nt
               write (6,'(''unconstrained steps = '',t30,i14)') nsteptr+1 ! one for tau=0
               answer=resstring(tau0)
               write (6,'(a90)') answer(1)
               write (6,'(a90)') answer(2)
            endif
            do while(.true.)
               call pop(3,w1,empty)
               if (empty) exit
            enddo
            do while(.true.)
               call pop(istin,w1,empty)
               if (empty) exit
               w1%wt0=abs(w1%psi)/w1%psi ! reset wt0 to start transient
               call push(istout,w1)
               call push(3,w1)
            enddo
            istout=istin
            istin=3-istout
            call cpu_time(timest)
            timest=timest-timebl0
            call addall(timest,tts)
            call cpu_time(timebl)
            tau0=tau0+dt
            n1=numstack(istin)
            call zerest
            call setidmc(5)
            do j=1,nsteptr
               call step1(istin,istout,ts,to)
               n1=numstack(istin)
               call addall(n1,nt)
               call bcast(nt)
               if (nt.gt.popup*nwalk.or.nt.lt.popdn*nwalk) call branch(istin,istout,1.0_r8*nwalk/nt)
               call checkpop(istin,tlb)
               call calcop(istin,istout,to)
               call update   !collect block averages
               if (myrank().eq.0) then
                  answer=resstring(tau0)
                  write (6,'(a90)') answer(1)
                  write (6,'(a90)') answer(2)
               endif
               tau0=tau0+dt
            enddo
            do while(.true.)
               call pop(istin,w1,empty)
               if (empty) exit
            enddo
            do while(.true.)
               call pop(3,w1,empty)
               if (empty) exit
               call hpsi(w1,.false.)
               call push(istin,w1)
            enddo !
            if (myrank().eq.0) write (6,'(''Time for CP steps (sec) (min) = '',t40,2f11.3)') tts/nproc(),tts/nproc()/60.0_r8
            call cpu_time(timest)
            timest=timest-timebl
            call addall(timest,tts)
            if (myrank().eq.0) write (6,'(''Time for UP steps (sec) (min) = '',t40,2f11.3)') tts/nproc(),tts/nproc()/60.0_r8
            call cpu_time(timest)
            timest=timest-timebl0
            call addall(timest,tts)
            if (myrank().eq.0) write (6,'(''Time for block (sec) (min) = '',t40,2f11.3)') tts/nproc(),tts/nproc()/60.0_r8
            call cpu_time(timest)
            timest=timest-time2
            call addall(timest,tts)
            if (myrank().eq.0) write (6,'(''Elapsed block time (sec) (min) = '',t40,2f11.3)') tts/nproc(),tts/nproc()/60.0_r8
         enddo
      endif
      if (idmc.eq.-1) call updatepar(ene) ! parameters optimization
      call cpu_time(timeop)
      t0=timeop-time2
      call addall(t0,t0t)
      if (myrank().eq.0) write (6,'(/,''Time for steps ='',f10.3,'' seconds, '',f10.3,'' minutes'')') t0t/nproc(),t0t/60.0/nproc()
      n1=numstack(istin)
      call gather(n1,numfig)
      call bcast(numfig)
      nt=sum(numfig)
   enddo
   if (nav+neq.ne.0.and.idmc.ne.5) then
      rewind 9
      if (myrank().eq.0) write(9,'(i10,f15.8)') nt,tau
      close(9)
      call barrier
      do k=0,nboss()-1
         if (k.eq.myrank()) then
            open(unit=9,position='append')
            do i=1,numfig(k)
               call pop(istin,w1,empty)
               if (.not.empty) then
                  write (9,'(6e15.7)') w1%x
                  write (9,'(6e15.7)') w1%sp
                  write (9,'(4e15.7)') w1%weight,w1%wt0
                  write (9,'(4e15.7)') w1%psi,w1%psig
                  write (9,'(i20)') w1%irn
               endif
            enddo
            close(9)
         endif
         call barrier
      enddo
   endif
   call cpu_time(time1)
   t0=time1-time0
   call addall(t0,t0t)
   if (myrank().eq.0) then
      write (6,'(''Total job time ='',f15.8,'' seconds, '',f15.8,'' minutes, ''  &
       & ,f15.8,''hours'')') t0t/nproc(),t0t/nproc()/60.0_r8,t0t/nproc()/3600.0_r8
      write (6,'(''Total cpu-time ='',f15.8,'' minutes, '',f15.8,'' hours'')') t0t/60.0_r8,t0t/3600.0_r8
      write (6,'(''Number of cpus ='',i10)') nproc()
   endif
   call barrier
   call printmpilog
   call done
   if (myrank().eq.0) then
    call system_clock(s_time1,clock_rate)
    write(6,*) 'Total time= ',real(s_time1-s_time0)/real(clock_rate),'seconds'
    write (6,*) 'Finished!'
   endif
   end program afnuclear
