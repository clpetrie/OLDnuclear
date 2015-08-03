module jastrow
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, save :: el2,dc,dtn,dr,dr2,hb,el,eli
   integer(kind=i4), private, save :: ntab,lpot
   real(kind=r8), private, save, allocatable :: uoptab(:,:,:),duoptab(:,:,:),d2uoptab(:,:,:),ddtn(:)
   real(kind=r8), private, save :: acn,at,as,ast,atn,att,als,bt,bs,bst,btn,btt,bls,rho_snm
   real(kind=r8), private, save :: rscpp,rscnn,bpp,bnn
   real(kind=r8), private, parameter :: dp=0.0025_r8
   real(kind=r8), private, save :: q1c,q2c,q1p,q2p
   integer(kind=i4), private, save :: ijas
   real(kind=r8), private, save, allocatable :: rsctni(:),ptni(:)
   logical, private, save :: fcsb
   logical, private, save :: fixdc
contains
   subroutine jasinit(elin,ntabin,hbin,lpotin,fcsbin,f3in)
   use mympi
   use nucma
   real(kind=r8) :: elin,hbin
   integer(kind=i4) :: ntabin,lpotin
   integer(kind=i4) :: i,j
   logical :: fcsbin,f3in,f3
   if (elin.lt.0.0_r8) then
      el=0.0_r8
      eli=0.0_r8
      el2=-elin
   else
      el=elin
      eli=1.0_r8/el
      el2=0.5_r8*el
   endif
   ntab=ntabin
   hb=hbin
   lpot=lpotin
   f3=.true.
   fcsb=.true.
   fixdc=.true.  !readjust dc if larger than dtn
   allocate(rsctni(3),ptni(11))
   if (myrank().eq.0) then
      read (5,*) dc      !jastrow healing distance
      read (5,*) dtn
      rho_snm=0.01_r8
      read (5,*) acn
      read (5,*) at
      read (5,*) as
      read (5,*) ast
      read (5,*) atn
      read (5,*) att
      read (5,*) als
      read (5,*) bt
      read (5,*) bs
      read (5,*) bst
      read (5,*) btn
      read (5,*) btt
      read (5,*) bls
      read (5,*) q1c
      read (5,*) q2c
      read (5,*) q1p
      read (5,*) q2p
      read (5,*) rsctni(1) ! rscal cent
      read (5,*) rsctni(2) ! rscal xpi
      read (5,*) rsctni(3) ! rscal delta
      read (5,*) ptni(1) ! eps_cent
      read (5,*) ptni(2) ! eps_anti
      read (5,*) ptni(3) ! eps_tm
      read (5,*) ptni(4) ! eps vd1
      read (5,*) ptni(5) ! eps_vd2
      read (5,*) ptni(6) ! eps_vdc3
      read (5,*) ptni(7) ! eps_vec3
      read (5,*) ptni(8) ! eps_ve !CODY DELETE good to here
      read (5,*) ptni(9) ! eps_comm
      read (5,*) ptni(10) ! eps_comm-xd
      read (5,*) ptni(11) ! eps_comm-dd
      if (ptni(9).eq.0.0_r8.and.ptni(10).eq.0.0_r8.and.ptni(11).eq.0.0_r8) f3=.false.
      read (5,*) rscpp
      read (5,*) rscnn
      read (5,*) bpp
      read (5,*) bnn
      if (bpp.eq.0.0_r8.and.bnn.eq.0.0_r8) fcsb=.false.
      write (6,'(''jastrow healing distance dc ='',t40,f10.5)') dc
      write (6,'(''jastrow healing distance dtn ='',t40,f10.5)') dtn
      write (6,'(''jastrow rho ='',t40,f10.5)') rho_snm
      write (6,'(''potential quencher acn (v1)='',t40,f10.5)') acn
      write (6,'(''potential quencher at  (v2)='',t40,f10.5)') at
      write (6,'(''potential quencher as  (v3)='',t40,f10.5)') as
      write (6,'(''potential quencher ast (v4)='',t40,f10.5)') ast
      write (6,'(''potential quencher atn (v5)='',t40,f10.5)') atn
      write (6,'(''potential quencher atn (v6)='',t40,f10.5)') att
      write (6,'(''potential quencher als (v7 v8)='',t40,f10.5)') als
      write (6,'(''correlation quencher bt (f2)='',t40,f10.5)') bt
      write (6,'(''correlation quencher bs (f3)='',t40,f10.5)') bs
      write (6,'(''correlation quencher bst (f4)='',t40,f10.5)') bst
      write (6,'(''correlation quencher btn (f5)='',t40,f10.5)') btn
      write (6,'(''correlation quencher btn (f6)='',t40,f10.5)') btt
      write (6,'(''correlation quencher bls (f7)='',t40,f10.5)') bls
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''q1c parameter in f3c ='',t40,f10.5)') q1c
      write (6,'(''q2c parameter in f3c ='',t40,f10.5)') q2c
      write (6,'(''q1p parameter in f3p ='',t40,f10.5)') q1p
      write (6,'(''q2p parameter in f3p ='',t40,f10.5)') q2p
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''rscal parameter in Uijk cent ='',t40,f10.5)') rsctni(1)
      write (6,'(''rscal parameter in Uijk xpi ='',t40,f10.5)') rsctni(2)
      write (6,'(''rscal parameter in Uijk delta ='',t40,f10.5)') rsctni(3)
      write (6,'(''eps_c parameter in Uijk ='',t40,f10.5)') ptni(1)
      write (6,'(''eps_anti parameter in Uijk ='',t40,f10.5)') ptni(2)
      write (6,'(''eps_tm parameter in Uijk ='',t40,f10.5)') ptni(3)
      write (6,'(''eps_vd1 parameter in Uijk ='',t40,f10.5)') ptni(4)
      write (6,'(''eps_vd2 parameter in Uijk ='',t40,f10.5)') ptni(5)
      write (6,'(''eps_vdc3 parameter in Uijk ='',t40,f10.5)') ptni(6)
      write (6,'(''eps_vec3 parameter in Uijk ='',t40,f10.5)') ptni(7)
      write (6,'(''eps_ve parameter in Uijk ='',t40,f10.5)') ptni(8)
      write (6,'(''eps_comm parameter in Uijk ='',t40,f10.5)') ptni(9)
      write (6,'(''eps_comm-XD parameter in Uijk ='',t40,f10.5)') ptni(10)
      write (6,'(''eps_comm-DD parameter in Uijk ='',t40,f10.5)') ptni(11)
      write (6,'(''TNI-comm correlation included ='',t40,l10)') f3
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''CSB correlation rscalpp ='',t40,f10.5)') rscpp
      write (6,'(''CSB correlation rscalnn ='',t40,f10.5)') rscnn
      write (6,'(''CSB correlation quencher pp ='',t40,f10.5)') bpp
      write (6,'(''CSB correlation quencher nn ='',t40,f10.5)') bnn
      write (6,'(''CSB correlation included ='',t40,l10)') fcsb
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
   endif
   call bcast(dc)
   call bcast(dtn)
   if (dc.ge.el2) then
      dc=el2
      if (myrank().eq.0) write (6,'(''dc reset to '',f15.10)') dc
   endif
   if (dtn.ge.el2) then
      dtn=el2
      if (myrank().eq.0) write (6,'(''dtn reset to '',f15.10)') dtn
   endif
   if (dc.ge.dtn-0.2_r8) then
      if (fixdc) then
         dc=dtn-0.2_r8
         if (myrank().eq.0) write (6,'(''dc larger than dtn, reset to '',f15.10)') dc
      else
         dtn=dc+0.2_r8
         if (myrank().eq.0) write (6,'(''dtn smaller than dc, reset to '',f15.10)') dtn
      endif
   endif
   call bcast(rho_snm)
   call bcast(acn)
   call bcast(at)
   call bcast(as)
   call bcast(ast)
   call bcast(atn)
   call bcast(att)
   call bcast(als)
   call bcast(bt)
   call bcast(bs)
   call bcast(bst)
   call bcast(btn)
   call bcast(btt)
   call bcast(bls)
   call bcast(q1c)
   call bcast(q2c)
   call bcast(q1p)
   call bcast(q2p)
   call bcast(rsctni)
   call bcast(ptni)
   call bcast(rscpp)
   call bcast(rscnn)
   call bcast(bpp)
   call bcast(bnn)
   call bcast(f3)
   call bcast(fcsb)
   call bcast(rscpp)
   call bcast(rscpp)
   allocate(ddtn(18)) ! number of parameters for two-body correlations plus one
   call jasinit_nucma()
   ddtn=dtn
   if (myrank().eq.0) then
      write (6,'(''jastrow initialized'')')
      do i=1,ntab
         write(78,'(11e25.10)') i*dr,(uoptab(1,i,j),j=1,10)
      enddo 
   endif
   f3in=f3
   fcsbin=fcsb
   end subroutine jasinit


   subroutine jasinit_nucma()
   use nucma
   use mympi
   integer(kind=i4) :: j
   integer(kind=i4) :: lc,ls,lt,ll,nm
! We choose symmetric nuclear matter. For neutron drops, nm=2
   nm=1
   if (.not.allocated(uoptab)) allocate(uoptab(1,ntab,10),duoptab(1,ntab,10),d2uoptab(1,ntab,10))
   ddtn(1)=dtn
   lt=ntab
   lc=dble(lt)*dc/dtn
   ll=lc
   ls=lc
   dr=dtn/dble(ntab)
   call nmfts(lc,ls,lt,ll,acn,at,as,ast,atn,att,als,dtn,nm,rho_snm,uoptab(1,:,1:8),&
     &        duoptab(1,:,1:8),d2uoptab(1,:,1:8),ntab,dr,lpot,hb)

! Normalize non central correlations

   uoptab(1,:,2)=bt*uoptab(1,:,2)
   uoptab(1,:,3)=bs*uoptab(1,:,3)
   uoptab(1,:,4)=bst*uoptab(1,:,4)
   uoptab(1,:,5)=btn*uoptab(1,:,5)
   uoptab(1,:,6)=btt*uoptab(1,:,6)
   uoptab(1,:,7)=bls*uoptab(1,:,7)
   if (fcsb) then
      uoptab(1,:,9)=bpp/bt*uoptab(1,:,2)
      uoptab(1,:,10)=bnn/bt*uoptab(1,:,2)
   else
      uoptab(1,:,9)=0.0_r8
      uoptab(1,:,10)=0.0_r8
   endif
   do j=2,10
    uoptab(1,:,j)=uoptab(1,:,j)/uoptab(1,:,1)
   enddo
   uoptab(1,:,1)=-log(uoptab(1,:,1))
   end subroutine jasinit_nucma

   subroutine jastro_nucma(r,uop,itab)
   integer(kind=i4) :: index,itab
   real(kind=r8) :: r,uop(8),ddr,c1,c2,c3,c4
   ddr=r*ntab/ddtn(itab)
   index=ddr
   index=max(2,min(index,ntab-2))
   ddr=ddr-index
   c1=-ddr*(ddr-1.0_r8)*(ddr-2.0_r8)/6.0_r8
   c2=(ddr+1.0_r8)*(ddr-1.0_r8)*(ddr-2.0_r8)/2.0_r8
   c3=-(ddr+1.0_r8)*ddr*(ddr-2.0_r8)/2.0_r8
   c4=(ddr+1.0_r8)*ddr*(ddr-1.0_r8)/6.0_r8
   uop(1:8)=c1*uoptab(itab,index-1,1:8)+c2*uoptab(itab,index,1:8) &
           +c3*uoptab(itab,index+1,1:8)+c4*uoptab(itab,index+2,1:8)
   if (r.gt.dc) then
      uop(1:4)=0.0_r8
   endif
   if (r.gt.dtn) then
      uop(5:6)=0.0_r8
   endif
   end subroutine jastro_nucma

   subroutine jastro_nucmacsb(r,fpp,fnn,itab)
   integer(kind=i4) :: index,itab
   real(kind=r8) :: r,fpp,fnn,ddr,c1,c2,c3,c4
   if (fcsb) then
      ddr=r*rscpp*ntab/ddtn(itab)
      index=ddr
      index=max(2,min(index,ntab-2))
      ddr=ddr-index
      c1=-ddr*(ddr-1.0_r8)*(ddr-2.0_r8)/6.0_r8
      c2=(ddr+1.0_r8)*(ddr-1.0_r8)*(ddr-2.0_r8)/2.0_r8
      c3=-(ddr+1.0_r8)*ddr*(ddr-2.0_r8)/2.0_r8
      c4=(ddr+1.0_r8)*ddr*(ddr-1.0_r8)/6.0_r8
      fpp=c1*uoptab(itab,index-1,9)+c2*uoptab(itab,index,9) &
         +c3*uoptab(itab,index+1,9)+c4*uoptab(itab,index+2,9)
      if (r*rscpp.gt.dc) then
         fpp=0.0_r8
      endif
      ddr=r*rscnn*ntab/ddtn(itab)
      index=ddr
      index=max(2,min(index,ntab-2))
      ddr=ddr-index
      c1=-ddr*(ddr-1.0_r8)*(ddr-2.0_r8)/6.0_r8
      c2=(ddr+1.0_r8)*(ddr-1.0_r8)*(ddr-2.0_r8)/2.0_r8
      c3=-(ddr+1.0_r8)*ddr*(ddr-2.0_r8)/2.0_r8
      c4=(ddr+1.0_r8)*ddr*(ddr-1.0_r8)/6.0_r8
      fnn=c1*uoptab(itab,index-1,10)+c2*uoptab(itab,index,10) &
         +c3*uoptab(itab,index+1,10)+c4*uoptab(itab,index+2,10)
      if (r*rscnn.gt.dc) then
         fnn=0.0_r8
      endif
   else
      fpp=0.0_r8
      fnn=0.0_r8
   endif
   end subroutine jastro_nucmacsb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the following are subrotuines needed by the optimization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine f3p(n)
   integer(kind=i4) :: n
   select case (n)
      case(1)
         q1c=q1c+dp
      case(2)
         q2c=q2c+dp
      case(3)
         q1p=q1p+dp
      case(4)
         q2p=q2p+dp
      case(5:7)
         rsctni(n-4)=rsctni(n-4)+dp
      case(8:18)
         ptni(n-7)=ptni(n-7)+dp
      case default 
         write(6,'(''Error in f3p!!!'')')
   end select
   end subroutine f3p

   subroutine f3m(n)
   integer(kind=i4) :: n
   select case (n)
      case(1)
         q1c=q1c-dp
      case(2)
         q2c=q2c-dp
      case(3)
         q1p=q1p-dp
      case(4)
         q2p=q2p-dp
      case(5:7)
         rsctni(n-4)=rsctni(n-4)-dp
      case(8:18)
         ptni(n-7)=ptni(n-7)-dp
      case default 
         write(6,'(''Error in f3m!!!'')')
   end select
   end subroutine f3m

   subroutine setcsb(rscppin,rscnnin)
   real(kind=r8) :: rscppin,rscnnin
   rscpp=rscppin
   rscnn=rscnnin
   end subroutine setcsb

   subroutine getcsb(rscppin,rscnnin)
   real(kind=r8) :: rscppin,rscnnin
   rscppin=rscpp
   rscnnin=rscnn
   end subroutine getcsb

   subroutine jasinitopt
   use nucma
   use mympi
   integer(kind=i4) :: j,i
   integer(kind=i4) :: lc,lt,nm,ipar,npar
   real(kind=r8) :: dc0,dtn0,acn0,at0,as0,ast0,atn0,att0,als0,bt0,bs0,bst0,btn0,btt0,bls0
   real(kind=r8) :: rscpp0,rscnn0,bpp0,bnn0
   real(kind=r8), parameter :: tiny=1.0e-10_r8
   nm=1
! the variational parameters are:
! dc,dtn,acn,at,as,ast,atn,att,als,bt,bs,bst,btn,btt,bls
! bpp,bnn
! build a table for each variational parameter
   deallocate(uoptab,duoptab,d2uoptab)
   npar=18  ! number of different tables to construct, add one for the original parameters
   allocate(uoptab(npar,ntab,10),duoptab(npar,ntab,10),d2uoptab(npar,ntab,10))
   lt=ntab
   do ipar=1,npar
      dc0=dc
      dtn0=dtn
      acn0=acn
      at0=at
      as0=as
      ast0=ast
      atn0=atn
      att0=att
      als0=als
      bt0=bt
      bs0=bs
      bst0=bst
      btn0=btn
      btt0=btt
      bls0=bls
      rscpp0=rscpp
      rscnn0=rscnn
      bpp0=bpp
      bnn0=bnn
      lc=dble(lt)*dc/dtn
      dr=dtn/dble(ntab)
      select case (ipar)
      case (1) 
         lc=dble(lt)*dc/dtn
         dr=dtn/dble(ntab)
      case (2) ! dc
         dc=dc0+dp
         lc=dble(lt)*dc/dtn
         dr=dtn/dble(ntab)
      case (3) ! dtn
         dtn=dtn0+dp
         lc=dble(lt)*dc/dtn
         dr=dtn/dble(ntab)
      case (4) ! acn
         acn=acn0+dp
      case (5) ! at
         at=at0+dp
      case (6) ! as
         as=as0+dp
      case (7) ! ast
         ast=ast0+dp
      case (8) ! atn
         atn=atn0+dp
      case (9) ! att
         att=att0+dp
      case (10) ! als
         als=als0+dp
      case (11) ! bt
         bt=bt0+dp
      case (12) ! bs
         bs=bs0+dp
      case (13) ! bst
         bst=bst0+dp
      case (14) ! btn
         btn=btn0+dp
      case (15) ! btt
         btt=btt0+dp
      case (16) ! bls
         bls=bls0+dp
      case (17) ! bpp
         bpp=bpp0+dp
      case (18) ! bnn
         bnn=bnn0+dp
      end select
      ddtn(ipar)=dtn
      if (ipar.le.10) then
!        if (myrank().eq.0) write (6,'(''Calculating table '',i5)') ipar
         call nmfts(lc,lc,lt,lc,acn,at,as,ast,atn,att,als,dtn,nm,rho_snm,uoptab(ipar,:,1:8),&
              duoptab(ipar,:,1:8),d2uoptab(ipar,:,1:8),ntab,dr,lpot,hb)
         uoptab(ipar,:,9)=uoptab(ipar,:,2) ! start with pp correlation equal to ftau
         uoptab(ipar,:,10)=uoptab(ipar,:,2) ! start with nn correlation equal to ftau
!        if (myrank().eq.0) write (6,'(''Done table '',i5)') ipar
      else
         uoptab(ipar,:,1)=uoptab(1,:,1)
         uoptab(ipar,:,2)=uoptab(1,:,2)/(bt0+tiny)
         uoptab(ipar,:,3)=uoptab(1,:,3)/(bs0+tiny)
         uoptab(ipar,:,4)=uoptab(1,:,4)/(bst0+tiny)
         uoptab(ipar,:,5)=uoptab(1,:,5)/(btn0+tiny)
         uoptab(ipar,:,6)=uoptab(1,:,6)/(btt0+tiny)
         uoptab(ipar,:,7)=uoptab(1,:,7)/(bls0+tiny)
         uoptab(ipar,:,9)=uoptab(1,:,9)/(bpp0+tiny)
         uoptab(ipar,:,10)=uoptab(1,:,10)/(bnn0+tiny)
      endif
      dc=dc0
      dtn=dtn0
      acn=acn0
      at=at0
      as=as0
      ast=ast0
      atn=atn0
      att=att0
      als=als0
      uoptab(ipar,:,2)=bt*uoptab(ipar,:,2)
      bt=bt0
      uoptab(ipar,:,3)=bs*uoptab(ipar,:,3)
      bs=bs0
      uoptab(ipar,:,4)=bst*uoptab(ipar,:,4)
      bst=bst0
      uoptab(ipar,:,5)=btn*uoptab(ipar,:,5)
      btn=btn0
      uoptab(ipar,:,6)=btt*uoptab(ipar,:,6)
      btt=btt0
      uoptab(ipar,:,7)=bls*uoptab(ipar,:,7)
      bls=bls0
      if (fcsb) then
         uoptab(ipar,:,9)=bpp*uoptab(ipar,:,9)
         uoptab(ipar,:,10)=bnn*uoptab(ipar,:,10)
      else
         uoptab(ipar,:,9:10)=0.0_r8
      endif
      bpp=bpp0
      bnn=bnn0
   enddo
   do j=2,10
      uoptab(:,:,j)=uoptab(:,:,j)/uoptab(:,:,1)
   enddo
   uoptab(:,:,1)=-log(uoptab(:,:,1))
!  do ipar=1,npar
!     do i=1,ntab
!        if (myrank().eq.0) then
!           write(79,'(2i10,10f25.15)') ipar,i,(uoptab(ipar,i,j),j=1,10)
!        endif
!     enddo
!     write(79,*) ' '
!     write(79,*) ' '
!  enddo
   do i=1,ntab
      do j=1,7
         do ipar=2,npar
!           if (abs(uoptab(ipar,i,j)-uoptab(1,i,j)).gt.0.1_r8) then
!              write(6,'(''big change in correlations in jasinitopt'')')
!              print*,i,j,ipar,uoptab(ipar,i,j),uoptab(1,i,j)
!              stop
!           endif
            if (isnan(uoptab(ipar,i,j)).or.isnan(uoptab(1,i,j))) then
               write(6,'(''NaN in correlations in jasinitopt'')')
               print*,i,j,ipar,uoptab(ipar,i,j),uoptab(1,i,j)
               stop
            endif
         enddo
      enddo
   enddo
   end subroutine jasinitopt

   subroutine setjasparam(params)
   use mympi
   real(kind=r8) :: params(:)
   dc=params(1)
   dtn=params(2)
   acn=params(3)
   at=params(4)
   as=params(5)
   ast=params(6)
   atn=params(7)
   att=params(8)
   als=params(9)
   bt=params(10)
   bs=params(11)
   bst=params(12)
   btn=params(13)
   btt=params(14)
   bls=params(15)
   if (dc.ge.el2) then
      dc=el2
      params(1)=dc
      if (myrank().eq.0) write (6,'(''dc reset to '',f15.10)') dc
   endif
   if (dtn.ge.el2) then
      dtn=el2
      params(2)=dtn
      if (myrank().eq.0) write (6,'(''dtn reset to '',f15.10)') dtn
   endif
   if (dc.ge.dtn-0.2_r8) then
      if (fixdc) then
         dc=dtn-0.2_r8
         params(1)=dc
         if (myrank().eq.0) write (6,'(''dc larger than dtn, reset to '',f15.10)') dc
      else
         dtn=dc+0.2_r8
         params(2)=dtn
         if (myrank().eq.0) write (6,'(''dtn smaller than dc, reset to '',f15.10)') dtn
      endif
   endif
   q1c=params(16)
   q2c=params(17)
   q1p=params(18)
   q2p=params(19)
   rsctni(1:3)=params(20:22)
   ptni(1:11)=params(23:33)
   rscpp=params(34)
   rscnn=params(35)
   if (rscpp*dc.ge.el2) then
      rscpp=el2/dc
      params(34)=rscpp
      if (myrank().eq.0) write (6,'(''rscalpp reset to '',f15.10)') rscpp
   endif
   if (rscnn*dc.ge.el2) then
      rscnn=el2/dc
      params(35)=rscnn
      if (myrank().eq.0) write (6,'(''rscalnn reset to '',f15.10)') rscnn
   endif
   bpp=params(36)
   bnn=params(37)
   end subroutine setjasparam

   subroutine getjasparam(nparam,params)
   use mympi
   integer(kind=i4) :: nparam
   real(kind=r8), pointer :: params(:)
   if (dc.ge.el2) then
      dc=el2
      if (myrank().eq.0) write (6,'(''dc reset to '',f15.10)') dc
   endif
   if (dtn.ge.el2) then
      dtn=el2
      if (myrank().eq.0) write (6,'(''dtn reset to '',f15.10)') dtn
   endif
   if (dc.ge.dtn-0.2_r8) then
      if (fixdc) then
         dc=dtn-0.2_r8
         if (myrank().eq.0) write (6,'(''dc larger than dtn, reset to '',f15.10)') dc
      else
         dtn=dc+0.2_r8
         if (myrank().eq.0) write (6,'(''dtn smaller than dc, reset to '',f15.10)') dtn
      endif
   endif
   if (rscpp*dc.ge.el2) then
      rscpp=el2/dc
      if (myrank().eq.0) write (6,'(''rscalpp reset to '',f15.10)') rscpp
   endif
   if (rscnn*dc.ge.el2) then
      rscnn=el2/dc
      if (myrank().eq.0) write (6,'(''rscalnn reset to '',f15.10)') rscnn
   endif
   nparam=37
   allocate(params(nparam))
   params(1)=dc
   params(2)=dtn
   params(3)=acn
   params(4)=at
   params(5)=as
   params(6)=ast
   params(7)=atn
   params(8)=att
   params(9)=als
   params(10)=bt
   params(11)=bs
   params(12)=bst
   params(13)=btn
   params(14)=btt
   params(15)=bls
   params(16)=q1c
   params(17)=q2c
   params(18)=q1p
   params(19)=q2p
   params(20:22)=rsctni(1:3)
   params(23:33)=ptni(1:11)
   params(34)=rscpp
   params(35)=rscnn
   params(36)=bpp
   params(37)=bnn
   end subroutine getjasparam

   subroutine hscor(npart,x,uc,sigma,sigmatau,tau,fls,fpp,fnn)
   integer(kind=i4) :: i,j,k,npart
   real(kind=r8) :: x(3,npart)
   real(kind=r8) :: sigma(3,npart,3,npart),sigmatau(3,npart,3,npart)
   real(kind=r8) :: tau(npart,npart),fpp(npart,npart),fnn(npart,npart)
   real(kind=r8) :: dx(3),r
   real(kind=r8) :: uu(8),uc,us(3,3),ust(3,3),fls(npart,npart)
   real(kind=r8) :: u3c,dot1,dot2,dot3,u3p
   real(kind=r8), dimension(3) :: dxij,dxjk,dxik
   real(kind=r8) :: rij,rik,rjk
   real(kind=r8) :: dxp(3,npart,npart),rp(npart,npart),exq2p(npart,npart)
   real(kind=r8) :: exq2c(npart,npart)
   dxp(:,npart,npart)=0.0_r8
   rp(npart,npart)=0.0_r8
   exq2p(npart,npart)=0.0_r8
   exq2c(npart,npart)=0.0_r8
   do i=1,npart-1
      dxp(:,i,i)=0.0_r8
      rp(i,i)=0.0_r8
      exq2p(i,i)=0.0_r8
      exq2c(i,i)=0.0_r8
      do j=i+1,npart
         dx=x(:,i)-x(:,j)
         dx=dx-el*nint(dx*eli)
         r=sqrt(sum(dx**2))
         dx=dx/r
         dxp(:,i,j)=dx
         dxp(:,j,i)=-dx
         rp(i,j)=r
         rp(j,i)=r
         exq2p(i,j)=exp(-q2p*r)
         exq2p(j,i)=exq2p(i,j)
         exq2c(i,j)=exp(-q2c*r)
         exq2c(j,i)=exq2c(i,j)
      enddo
   enddo
   uc=0.0_r8
   sigma=0.0_r8
   sigmatau=0.0_r8
   tau=0.0_r8
   u3p=1.0_r8
   fls=0.0_r8
   fpp=0.0_r8
   fnn=0.0_r8
   do j=2,npart
      do i=1,j-1
         us=0.0_r8
         ust=0.0_r8
         dx=dxp(:,i,j)
         rij=rp(i,j)
         call jastro_nucma(rij,uu,ijas) ! the use of ijas is ugly but easy
         if (fcsb) then
            call jastro_nucmacsb(rij,fpp(i,j),fnn(i,j),ijas)
            fpp(j,i)=fpp(i,j)
            fnn(j,i)=fnn(i,j)
         endif
         u3p=1.0_r8
         if (q1p.ne.0.0_r8) then
            do k=1,npart
               dxjk=dxp(:,j,k)
               dxik=dxp(:,i,k)
               dot1=sum(dxik*dxjk)
               u3p=u3p*(1.0_r8-q1p*(1.0_r8-dot1) &
                  *exq2p(i,j)*exq2p(j,k)*exq2p(i,k))
            enddo
         endif
         uu(2:8)=uu(2:8)*u3p
         uc=uc+uu(1)
         us(1,1)=us(1,1)+uu(3)+uu(5)*(3.0_r8*dx(1)*dx(1)-1.0_r8)
         us(1,2)=us(1,2)+uu(5)*3.0_r8*dx(1)*dx(2)
         us(1,3)=us(1,3)+uu(5)*3.0_r8*dx(1)*dx(3)
         us(2,2)=us(2,2)+uu(3)+uu(5)*(3.0_r8*dx(2)*dx(2)-1.0_r8)
         us(2,3)=us(2,3)+uu(5)*3.0_r8*dx(2)*dx(3)
         us(3,3)=us(3,3)+uu(3)+uu(5)*(3.0_r8*dx(3)*dx(3)-1.0_r8)
         ust(1,1)=ust(1,1)+uu(4)+uu(6)*(3.0_r8*dx(1)*dx(1)-1.0_r8)
         ust(1,2)=ust(1,2)+uu(6)*3.0_r8*dx(1)*dx(2)
         ust(1,3)=ust(1,3)+uu(6)*3.0_r8*dx(1)*dx(3)
         ust(2,2)=ust(2,2)+uu(4)+uu(6)*(3.0_r8*dx(2)*dx(2)-1.0_r8)
         ust(2,3)=ust(2,3)+uu(6)*3.0_r8*dx(2)*dx(3)
         ust(3,3)=ust(3,3)+uu(4)+uu(6)*(3.0_r8*dx(3)*dx(3)-1.0_r8)
         tau(i,j)=tau(i,j)+uu(2)
         us(2,1)=us(1,2)
         us(3,1)=us(1,3)
         us(3,2)=us(2,3)
         ust(2,1)=ust(1,2)
         ust(3,1)=ust(1,3)
         ust(3,2)=ust(2,3)
         tau(j,i)=tau(i,j)
         sigma(:,i,:,j)=us(:,:)
         sigma(:,j,:,i)=us(:,:)
         sigmatau(:,i,:,j)=ust(:,:)
         sigmatau(:,j,:,i)=ust(:,:)
         fls(i,j)=uu(7)
         fls(j,i)=uu(7)
      enddo
   enddo
   u3c=1.0_r8
   if (q1c.ne.0.0_r8) then
      do i=1,npart-2
         do j=i+1,npart-1
            dxij=dxp(:,i,j)*rp(i,j)
            do k=j+1,npart
               dxjk=dxp(:,j,k)*rp(j,k)
               dxik=dxp(:,i,k)*rp(i,k)
               dot1=sum(dxij*dxik)
               dot2=-sum(dxij*dxjk)
               dot3=sum(dxik*dxjk)
               u3c=u3c*(1.0_r8+q1c*dot1*dot2*dot3 &
                  *exq2c(i,j)*exq2c(j,k)*exq2c(i,k))
            enddo
         enddo
      enddo
   endif
   u3c=-log(u3c)
   uc=uc+u3c
   end subroutine hscor

   subroutine tniacor(npart,x,uc,sigmatau,tau,xpi,xd,cf3xx,cf3xd,cf3dd)
   use v3bpot
   integer(kind=i4) :: npart
   real(kind=r8) :: x(3,npart),uc
   real(kind=r8) :: sigmatau(3,npart,3,npart),tau(npart,npart),xpi(3,npart,3,npart),xd(3,npart,3,npart)
   real(kind=r8), dimension(3,npart,3,npart) :: a3st,a3sttm,a3stvd1,a3stvd2,a3stvdc3,a3stvec3
   real(kind=r8) :: dummy,cf3xx,cf3xd,cf3dd
   real(kind=r8), dimension(npart,npart) :: atau
   call hstnimat(x,uc,a3st,a3sttm,a3stvd1,a3stvd2,a3stvdc3,a3stvec3,atau,dummy,dummy,dummy,xpi,xd,rsctni)
   uc=ptni(1)*uc
   sigmatau=ptni(2)*a3st+ptni(3)*a3sttm+ptni(4)*a3stvd1+ptni(5)*a3stvd2+ptni(6)*a3stvdc3+ptni(7)*a3stvec3
   tau=ptni(8)*atau
   cf3xx=ptni(9)
   cf3xd=ptni(10)
   cf3dd=ptni(11)
   end subroutine tniacor

   subroutine setijas(n) ! use a different jastrow
   integer(kind=i4) :: n
   ijas=n
   end subroutine setijas
end module jastrow
