module ioconfs
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
contains
   subroutine readconfs(filename,tau,nfac,hbar,npart)
   use stack
   use mympi
   use estimator
   use wavefunction
   use gofr
   integer(kind=i4) :: n1,k,ndis,ito,istin,istout
   integer :: lpot,idmc,npart
   real(kind=r8) :: tau,nfac,wtg,hbar
   real(kind=r8) :: en,vn,topt
   character(len=*) :: filename
   real :: time,time0,tim0,tim1
   logical :: empty
   real(kind=r8) :: rad,radp,radn,wtt
   complex(kind=r8) :: enc,wtc,wt0
   call cpu_time(time0)
   if (myrank().eq.0) then
      open(unit=9,file=filename)
      write (6,'(''readconfs: Reading walkers from file '',a20)') filename
      rewind 9
      read (9,'(i10,f15.8,2i10)') n1,tau,lpot,idmc
      write (6,'(''readconfs: Number of configurations ='',t60,i10)') n1
      write (6,'(''readconfs: Configurations at tau ='',t55,f15.10)') tau
      write (6,'(''readconfs: lpot ='',t60,i10)') lpot
      write (6,'(''readconfs: idmc ='',t60,i10)') idmc
! make sure that the number of confs is a multiple of the number of cpus
      ndis=mod(n1,nboss())
!ndis=0
      n1=n1-ndis
      write (6,'(''readconfs: Number of configurations discarded ='',t60,i10)') ndis
      write (6,'(''readconfs: New number of configurations ='',t60,i10)') n1
   endif
   call bcast(n1)
   call bcast(tau)
   do k=1,n1
      if (myrank().eq.0) then
        read (9,'(6e15.7)') w1%x
        read (9,'(6e15.7)') w1%sp
        read (9,'(4e15.7)') w1%weight,w1%wt0
        read (9,'(4e15.7)') w1%psi0,w1%psig0
        read (9,'(i20)') w1%irn
        call push(1,w1)
      endif
      call barrier
      ito=mod(k-1,nboss())
      call movewalkers(1,1,0,ito,k)
   enddo
   istin=1
   istout=2
   n1=numstack(istin)
   call barrier
!  write (6,'(''myrank, n1 '',2i10)') myrank(),n1
   if (myrank().eq.0) write (6,'(''readconfs: Done!'')')
   call cpu_time(time)
   time=time-time0
! now each cpu compute averages
   topt=0.0
   do k=1,n1
      call cpu_time(tim0)
      call pop(istin,w1,empty)
      if (empty) exit         
      call hpsi(w1,.true.)
      wt0=1.0_r8
!     en=-hbar*real(w1%d2psi)+real(w1%v)
!     vn=real(w1%v)
      wtg=w1%weight ! not workign for transient, fix it !!!
      wtt=1.0_r8
      wtc=1.0_r8
! not working for idmc=5 and/or transient, fix it !!!
!     enc=-hbar*w1%d2psi+w1%v
!     call addval(1,enc*nfac,wtc)
!     call addval(2,wtc,(1.0_r8,0.0_r8))
!     call addval(1,en*nfac*wtt,wtg)
!     call addval(2,(en-vn)*nfac*wtt,wtg)
!     call addval(3,vn*nfac*wtt,wtg)
!     en=-hbar*aimag(w1%d2psi)+aimag(w1%v)
!     vn=aimag(w1%v)
!     call addval(4,en*nfac*wtt,wtg)
!     call addval(5,(en-vn)*nfac*wtt,wtg)
!     call addval(6,vn*nfac*wtt,wtg)
!     call addval(7,real(w1%v8all(1))*nfac*wtt,wtg)
!     call addval(8,real(w1%v8all(2))*nfac*wtt,wtg)
!     call addval(9,real(w1%v8all(3))*nfac*wtt,wtg)
!     call addval(10,real(w1%v8all(4))*nfac*wtt,wtg)
!     call addval(11,real(w1%v8all(5))*nfac*wtt,wtg)
!     call addval(12,real(w1%v8all(6))*nfac*wtt,wtg)
!     call addval(13,real(w1%v8all(7))*nfac*wtt,wtg)
!     call addval(14,real(w1%v8all(8))*nfac*wtt,wtg)   
      en=-hbar*real(w1%d2psi)+real(sum(w1%v8all(1:6)))
!     call addval(15,en*nfac*wtt,wtg)   
!     call addval(16,real(w1%vcoul)*nfac*wtt,wtg)
!     call getradii(w1,rad,radp,radn,npart)
!     call addval(17,radp*wtt,wtg)
!     call addval(18,radn*wtt,wtg)
      call cpu_time(tim1)
      topt=topt+tim1-tim0
   enddo
   end subroutine readconfs

   subroutine writeconfs(filename,istin,istout,nt,tau,numfig,lpot,idmc)
   use stack
   use mympi
   character(len=*) :: filename
   integer(kind=i4) :: istin,istout,k,i,numfig(:),nt,lpot,idmc
   real(kind=r8) :: tau
   logical :: empty
   call barrier
   if (myrank().eq.0) then
      open(unit=9,file=filename)
      rewind 9
      write(9,'(i10,f15.8,2i10)') nt,tau,lpot,idmc
      close(9)
   endif
   call barrier
   do k=0,nboss()-1
      if (k.eq.myrank()) then
         open(unit=9,file=filename,position='append')
         do i=1,numfig(k+1)
            call pop(istin,w1,empty)
            if (.not.empty) then
               write (9,'(6e15.7)') w1%x
               write (9,'(6e15.7)') w1%sp
               write (9,'(4e15.7)') w1%weight,w1%wt0
               write (9,'(4e15.7)') w1%psi,w1%psig
               write (9,'(i20)') w1%irn
               call push(istout,w1)
            endif
         enddo
         close(9)
      endif
      call barrier
   enddo
   istout=istin
   istin=3-istout
   end subroutine writeconfs

   subroutine printout(filename,istin,istout,nt,tau,tau0,numfig,hbar)
   use stack
   use mympi
   character(len=*) :: filename
   integer(kind=i4) :: istin,istout,k,i,numfig(:),nt
   real(kind=r8) :: tau,tau0,hbar
   logical :: empty
   call barrier
   if (myrank().eq.0) then
      open(unit=9,file=filename)
      rewind 9
      write(9,'(i10,2f15.8)') nt,tau,tau0
      close(9)
   endif
   call barrier
   do k=0,nboss()-1
      if (k.eq.myrank()) then
         open(unit=9,file=filename,position='append')
         do i=1,numfig(k+1)
            call pop(istin,w1,empty)
            if (.not.empty) then
               write(9,*) w1%psi 
               write(9,*) w1%psig
               write(9,*) -hbar*w1%d2psi
!              write(9,*) w1%v
               write(9,*) w1%wt0
               call push(istout,w1)
            endif
         enddo
         close(9)
      endif
      call barrier
   enddo
   istout=istin
   istin=3-istout
   end subroutine printout

   subroutine writewalk(id,nunc,enc,wtc)
   integer(kind=i4) :: id
   integer(kind=i4) :: nunc
   complex(kind=r8) :: enc,wtc
   open(unit=100+id,position='append')
   write(100+id,'(''nunc='',i5,4e19.10)') nunc,enc,wtc
   close(100+id)
   end subroutine writewalk

   subroutine readfiles(nfiles)
   use stack
   use mympi
   integer(kind=i4) :: nfiles,ncpus,idfil,n1,i
   real(kind=r8) :: tau
   character(len=20) :: filename
   logical :: master
   real :: tlb
   call barrier
   master=.false.
   ncpus=nproc()/nfiles
   idfil=int(myrank()/ncpus)+101
   if (mod(myrank(),ncpus).eq.0) master=.true.
   if (master) then
      write(filename,*) idfil
      filename='confs.'//trim(adjustl(filename))
      open(unit=9,file=filename)
      read (9,'(i10,f15.8)') n1,tau
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write(6,'(''Master id, file, nwalks, tau '',i6,'' '',a,i10,f15.10)') myrank(),filename,n1,tau
      do i=1,n1
         read (9,'(6e15.7)') w1%x
         read (9,'(6e15.7)') w1%sp
         read (9,'(4e15.7)') w1%weight,w1%wt0
         read (9,'(4e15.7)') w1%psi0,w1%psig0
         read (9,'(i20)') w1%irn
         call push(1,w1)
      enddo
   endif
   call barrier
   call checkpop(1,tlb)
   end subroutine readfiles
end module
