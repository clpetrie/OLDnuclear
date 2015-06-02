module operators
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, save ::  hbar
   integer(kind=i4), private, save :: idmc,npart,nobs
   real(kind=r8), private, save :: nfac,el
   logical, private, save :: pcoul
contains
   subroutine setop(idmcin,npartin,hbarin,elin,nobsin,pcoulin)
   real(kind=r8) :: hbarin,elin
   integer(kind=i4) :: idmcin,npartin,nobsin
   logical :: pcoulin
   idmc=idmcin
   npart=npartin
   hbar=hbarin
   if (elin.gt.0.0_r8) then
      nfac=1.0_r8/npart
      el=elin
      pcoul=.false.
   else
      nfac=1.0_r8
      el=0.0_r8
   endif
   nobs=nobsin
   pcoul=pcoulin
   end subroutine setop
 
   subroutine calcop(istin,istout,topt)
   use stack
   use random
   use wavefunction
   use estimator
   use mympi
   use gofr
!  use euclidean
   integer(kind=i4) :: istin,istout
   real :: tim1,tim0,topt
   real(kind=r8) :: ejf,rad,radp,radn
   complex(kind=r8) :: vn,ek,diff,v3n,en,wtg
   logical :: empty
   call cpu_time(tim0)
   topt=0.0d0
   do while (.true.)
      call pop(istin,w2,empty)
      if (empty) exit
      call hpsi(w2,.true.)
      ek=-hbar*w2%d2psi
! only terms that are propagated are included in the energy
      v3n=w2%tnic+w2%tni2piapr+w2%tni2piaxdpr+w2%tni2piaddpr+w2%tni2pitm+w2%tnivd1+w2%tnivd2+w2%tnive
      vn=sum(w2%v8all(1:7))+v3n+w2%vext
      if (pcoul) vn=vn+w2%vcoul
      wtg=w2%weight
      if (idmc.eq.5) then
         wtg=w2%psi/w2%psig*w2%wt0
      endif
      en=ek+vn
      ejf=-0.5_r8*hbar*(real(w2%d2psi)-sum(conjg(w2%dpsi)*w2%dpsi))+real(vn)
      call addval(1,rcst(en)*nfac,wtg)
      call addval(2,rcst(ek)*nfac,wtg)
      call addval(3,rcst(vn)*nfac,wtg)
      call addval(4,aimag(en)*nfac,wtg)
      call addval(5,aimag(ek)*nfac,wtg)
      call addval(6,aimag(vn)*nfac,wtg)
      call addval(7,rcst(w2%v8all(1))*nfac,wtg)
      call addval(8,rcst(w2%v8all(2))*nfac,wtg)
      call addval(9,rcst(w2%v8all(3))*nfac,wtg)
      call addval(10,rcst(w2%v8all(4))*nfac,wtg)
      call addval(11,rcst(w2%v8all(5))*nfac,wtg)
      call addval(12,rcst(w2%v8all(6))*nfac,wtg)
      call addval(13,rcst(w2%v8all(7))*nfac,wtg)
      call addval(14,rcst(w2%v8all(8))*nfac,wtg)
      call addval(15,rcst(sum(w2%v8all(1:6)))*nfac,wtg)
      call addval(16,rcst(sum(w2%v8all(:)))*nfac,wtg)
      call addval(17,rcst(w2%vcoul)*nfac,wtg)
      call addval(18,rcst(w2%tni2piapr)*nfac,wtg)
      call addval(19,rcst(w2%tni2piaxdpr)*nfac,wtg)
      call addval(20,rcst(w2%tni2piaddpr)*nfac,wtg)
      call addval(21,rcst(w2%tni2pic)*nfac,wtg)
      call addval(22,rcst(w2%tni2picxd)*nfac,wtg)
      call addval(23,rcst(w2%tni2picdd)*nfac,wtg)
      call addval(24,rcst(w2%tni2pitm)*nfac,wtg)
      call addval(25,rcst(w2%tnivd1)*nfac,wtg)
      call addval(26,rcst(w2%tnivd2)*nfac,wtg)
      call addval(27,rcst(w2%tnive)*nfac,wtg)
      call addval(28,w2%tnic*nfac,wtg)
! include everything to the energy
      v3n=v3n+w2%tni2pic+w2%tni2picxd+w2%tni2picdd
      vn=sum(w2%v8all(:))+v3n+w2%vcoul+w2%vext
      call addval(29,rcst(vn)*nfac,wtg)
      call addval(30,rcst(ek+vn)*nfac,wtg)
      call addval(31,rcst(w2%v8allpr(1))*nfac,wtg)
      call addval(32,rcst(w2%v8allpr(2))*nfac,wtg)
      call addval(33,rcst(w2%v8allpr(3))*nfac,wtg)
      call addval(34,rcst(w2%v8allpr(4))*nfac,wtg)
      call addval(35,rcst(w2%v8allpr(5))*nfac,wtg)
      call addval(36,rcst(w2%v8allpr(6))*nfac,wtg)
      call addval(37,rcst(w2%v8allpr(7))*nfac,wtg)
      call addval(38,rcst(w2%v8allpr(8))*nfac,wtg)
      call addval(39,rcst(sum(w2%v8allpr(1:6)))*nfac,wtg)
      call addval(40,rcst(sum(w2%v8allpr(:)))*nfac,wtg)
      call addval(41,rcst((sum(w2%v8all(1:6))-sum(w2%v8allpr(1:6))))*nfac,wtg)
      call addval(42,rcst((sum(w2%v8all(:))-sum(w2%v8allpr(:))))*nfac,wtg)
      call addval(43,rcst(w2%tni2pia)*nfac,wtg)
      call addval(44,rcst(w2%tni2piaxd)*nfac,wtg)
      call addval(45,rcst(w2%tni2piadd)*nfac,wtg)
      call addval(46,rcst((w2%tni2piapr+w2%tni2piaxdpr+w2%tni2piaddpr-w2%tni2pia-w2%tni2piaxd-w2%tni2piadd))*nfac,wtg)
      diff=w2%tni2piapr-w2%tni2pia-w2%tni2pic
      call addval(47,rcst(diff)*nfac,wtg)
      diff=w2%tni2piaxdpr-w2%tni2piaxd-w2%tni2picxd
      call addval(48,rcst(diff)*nfac,wtg)
      diff=w2%tni2piaddpr-w2%tni2piadd-w2%tni2picdd
      call addval(49,rcst(diff)*nfac,wtg)
      if (el.eq.0.0_r8) then 
         call getradii(w2,rad,radp,radn,npart)
      else
         rad=0.0_r8
         radp=0.0_r8
         radn=0.0_r8
      endif 
      call addval(50,radp,wtg)
      call addval(51,radn,wtg)
      call addval(52,rad,wtg)
      call addval(53,(radn-radp),wtg)
      call addval(54,w2%vext,wtg)
      select case (idmc)
         case (-1:0)
            call addval(nobs+1,ejf*nfac,real(w2%weight))
         case(5,6)
            v3n=w2%tnic+w2%tni2piapr+w2%tni2piaxdpr+w2%tni2piaddpr+w2%tni2pitm+w2%tnivd1+w2%tnivd2+w2%tnive
            vn=sum(w2%v8all(1:7))+v3n+w2%vext
            if (pcoul) vn=vn+w2%vcoul
            en=ek+vn
            call addval(nobs+10,en*wtg,1.0_r8)
            call addval(nobs+11,wtg,1.0_r8)
            call addval(nobs+12,aimag(en*wtg),1.0_r8)
            call addval(nobs+13,aimag(wtg),1.0_r8)
      end select
! calculate other operators here
      call addgofr(w2,real(wtg))
!     if (eucl) call overlap(w2)
!
      call push(istout,w2)
   enddo
   istout=istin
   istin=3-istin
   call cpu_time(tim1)
   topt=topt+tim1-tim0
   contains
      function rcst(x)
      complex(kind=r8) :: rcst,x
      if (idmc.eq.5) then
         rcst=x
      else
         rcst=real(x)
      endif
      return
      end function rcst
   end subroutine calcop

   subroutine writeall(istin,istout,tau,nfiles)
   use stack
   use mympi
   integer(kind=i4) :: istin,istout
   real(kind=r8) :: tau
   integer(kind=i4) :: n1,i,j
   logical :: empty,master,filex
   integer(kind=i4) :: nfiles,idfil,ncpus,bossid
   character(len=20) :: filename
   complex(kind=r8) :: psi,psig,wt0,ek,v3n,vn,en
! write all the observables into a file
   master=.false.
   ncpus=nproc()/nfiles
   idfil=int(myrank()/ncpus)
   if (mod(myrank(),ncpus).eq.0) master=.true.
   if (master) then
      write(filename,*) idfil
      filename='output.'//trim(adjustl(filename))
      inquire(file=filename,exist=filex)
!     open(unit=9,file=filename,position='append',form='unformatted')
      open(unit=9,file=filename,position='append')
      if (.not.filex.and.myrank().eq.0) then
! list what is being printed only on output.0
         write(9,'(''tau'')')
         write(9,'(''psi'')')
         write(9,'(''psig'')')
         write(9,'(''wt0'')')
         write(9,'(''totan el'')')
         write(9,'(''#'')')
      endif
      do i=myrank(),myrank()+ncpus-1
!        write(6,'(''master id, receiving from i '',2i10)') myrank(),i
         if (i.eq.myrank()) then
            n1=numstack(istin)
            do j=1,n1
               call pop(istin,w2,empty)
               ek=-hbar*w2%d2psi
! only terms that are propagated are included in the energy
               v3n=w2%tnic+w2%tni2piapr+w2%tni2piaxdpr+w2%tni2piaddpr+w2%tni2pitm+w2%tnivd1+w2%tnivd2+w2%tnive
               vn=sum(w2%v8all(1:7))+v3n+w2%vext
               if (pcoul) vn=vn+w2%vcoul
               en=ek+vn
               write(9,*) tau
               write(9,*) w2%psi
               write(9,*) w2%psig
               write(9,*) w2%wt0
               write(9,*) en*nfac
               call push(istout,w2)
            enddo
         else
            call recv(n1,i,1)
            do j=1,n1
               call recv(psi,i,1)
               call recv(psig,i,1)
               call recv(wt0,i,1)
               call recv(en,i,1)
               write(9,*) tau
               write(9,*) psi
               write(9,*) psig
               write(9,*) wt0
               write(9,*) en*nfac
            enddo
         endif
      enddo
      close(9)
   else
      bossid=int(myrank()/ncpus)*ncpus
!     write(6,'(''worker id, sending to i '',2i10)') myrank(),bossid
      n1=numstack(istin)
      call send(n1,bossid,1)
      do j=1,n1
         call pop(istin,w2,empty)
         ek=-hbar*w2%d2psi
! only terms that are propagated are included in the energy
         v3n=w2%tnic+w2%tni2piapr+w2%tni2piaxdpr+w2%tni2piaddpr+w2%tni2pitm+w2%tnivd1+w2%tnivd2+w2%tnive
         vn=sum(w2%v8all(1:7))+v3n+w2%vext
         if (pcoul) vn=vn+w2%vcoul
         en=ek+vn
         call send(w2%psi,bossid,1)
         call send(w2%psig,bossid,1)
         call send(w2%wt0,bossid,1)
         call send(en,bossid,1)
         call push(istout,w2)
      enddo
   endif
   istout=istin
   istin=3-istin
   end subroutine writeall
 


end module operators
