!$Id: nompi.f90,v 1.1.1.1 2013/08/26 20:57:32 nuclear Exp $
module mympi
!
! fake versions that do nothing to run on single process w/o mpi
! load balancing mpi routines -- these are rewrites by K.E. Schmidt
! of the mpi routines written by Michael A. Lee and I. Lomonosov for the
! parallel version of the Schmidt and Lee electronic structure GFMC code
! 
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: irank,iproc,nleft,npart
   integer(kind=i4), private, allocatable, save :: numfig(:)
   logical, private, save :: nobal

interface bcast ! broadcast from process 0
   module procedure bcasti1,bcasti1d,bcasti2d,bcasti3d
   module procedure bcastr1,bcastr1d,bcastr2d,bcastr3d
   module procedure bcastc1,bcastc1d,bcastc2d,bcastc3d
   module procedure bcastl1,bcastl1d
end interface

interface addall ! return sum to process 0
   module procedure addalli1,addalli1d
   module procedure addallr1,addallr1d,addallr2d
   module procedure addallc1,addallc1d
end interface

interface gather ! gather to process 0
   module procedure gatheri1,gatheri1d
   module procedure gatherr1,gatherr1d
end interface

interface share ! a boss send something to all his workers
   module procedure sharew,sharei,sharer,sharel
end interface

contains
   subroutine init0(iworkerin) ! call this before anything else
   integer(kind=i4) :: ierror,iworkerin
   irank=0
   iproc=1
   nobal=.true.
   end subroutine init0

   subroutine init1(npartin,dobalin,ivmcin) !call after everyone has these
   integer(kind=i4) :: npartin
   real(kind=r8) :: dobalin
   logical :: ivmcin
   npart=npartin
!   dobal=dobalin
!   ivmc=ivmcin
   end subroutine init1

   subroutine done ! wrapper for finalize routine
   integer(kind=i4) :: ierror
   end subroutine done

   subroutine bcasti1(i)
   integer(kind=i4) :: i,ierror
   return
   end subroutine bcasti1

   subroutine bcasti1d(i)
   integer(kind=i4) :: i(:),ierror
   return
   end subroutine bcasti1d

   subroutine bcasti2d(i)
   integer(kind=i4) :: i(:,:),ierror
   return
   end subroutine bcasti2d

   subroutine bcasti3d(i)
   integer(kind=i4) :: i(:,:,:),ierror
   return
   end subroutine bcasti3d

   subroutine bcastr1d(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:)
   return
   end subroutine bcastr1d

   subroutine bcastr2d(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:,:)
   return
   end subroutine bcastr2d

   subroutine bcastr3d(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:,:,:)
   return
   end subroutine bcastr3d

   subroutine bcastr1(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r
   return
   end subroutine bcastr1

   subroutine bcastc1d(r)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: r(:)
   return
   end subroutine bcastc1d

   subroutine bcastc2d(r)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: r(:,:)
   return
   end subroutine bcastc2d

   subroutine bcastc3d(r)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: r(:,:,:)
   return
   end subroutine bcastc3d

   subroutine bcastc1(r)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: r
   return
   end subroutine bcastc1

   subroutine bcastl1(l)
   integer(kind=i4) :: ierror
   logical :: l
   return
   end subroutine bcastl1

   subroutine bcastl1d(l)
   integer(kind=i4) :: ierror
   logical :: l(:)
   end subroutine bcastl1d

   function myrank() ! which process am I?
   integer(kind=i4) :: myrank
   myrank=irank
   end function myrank

   function nproc() ! How many of use are there anyway?
   integer(kind=i4) :: nproc
   nproc=iproc
   end function nproc

   subroutine barrier ! wrapper for mpi_barrier
   integer(kind=i4) :: ierror
   end subroutine barrier

   subroutine addalli1(i,isum)
   integer(kind=i4) :: ierror,i,isum
   isum=i
   return
   end subroutine addalli1

   subroutine addalli1d(i,isum)
   integer(kind=i4) :: ierror,i(:),isum(:)
   isum=i
   return
   end subroutine addalli1d

   subroutine addallr1(r,rsum)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r,rsum
   rsum=r
   return
   end subroutine addallr1

   subroutine addallr1d(r,rsum)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:),rsum(:)
   rsum=r
   return
   end subroutine addallr1d

   subroutine addallr2d(r,rsum)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:,:),rsum(:,:)
   rsum=r
   return
   end subroutine addallr2d

   subroutine addallc1(c,csum)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c,csum
   csum=c
   end subroutine addallc1

   subroutine addallc1d(c,csum)
   complex(kind=r8) :: c(:),csum(:)
   integer(kind=i4) :: ierror
   csum=c
   end subroutine addallc1d

   subroutine gatheri1(i,igather)
   integer(kind=i4) :: i,igather(1:),ierror
   igather(1)=i
   return
   end subroutine gatheri1

   subroutine gatheri1d(i,igather)
   integer(kind=i4) :: i(:),igather(:,1:),ierror
   igather(:,1)=i
   return
   end subroutine gatheri1d

   subroutine gatherr1(r,rgather)
   real(kind=r8) :: r,rgather(1:)
   integer(kind=i4) :: ierror
   rgather(1)=r
   return
   end subroutine gatherr1

   subroutine gatherr1d(r,rgather)
   real(kind=r8) :: r(:),rgather(:,1:)
   integer(kind=i4) :: ierror
   rgather(:,1)=r
   return
   end subroutine gatherr1d

   subroutine movewalkers(istack,nwalk,ifrom,ito,id)
!
! routine to move nwalk walkers on istack from process ifrom to
! process ito. id is an arbitrary integer identifier
!
   integer(kind=i4) :: nwalk,ifrom,ito,id,i,istack
   return
   end subroutine movewalkers

   subroutine loadme(done,istack)
!
! processes should call this routine when they want more work
!
   logical :: done
   integer(kind=i4) :: istack
   done=.true.
   return
   end subroutine loadme

   subroutine loadcheck(istack)
!
! processes should call this routine to see if others need work
!
   integer(kind=i4) :: istack
   return
   end subroutine loadcheck

   subroutine load(istack)
   integer(kind=i4) :: istack
   return
   end subroutine load

   subroutine loadon
   return
   end subroutine loadon

   subroutine loadoff
   return
   end subroutine loadoff

   subroutine abort
   stop
   return
   end subroutine abort

   subroutine checkpop(istack,tlb)
   use stack
   integer(kind=i4) :: istack
   integer(kind=i4) :: n1,nt,ierror,nave
   integer(kind=i4) :: nmax,nmin
   real(kind=r8) :: frac
   real(kind=r8) :: tbal0,tbal1
   real :: tlb
   end subroutine checkpop

   subroutine printmpilog
   end subroutine printmpilog

   subroutine zerotimebal
!   tbal=0.0_r8
   end subroutine zerotimebal
   
   function nboss() ! How many of use are there anyway?
   integer(kind=i4) :: nboss
   nboss=1
   end function nboss

   function nworker() ! How many of use are there anyway?
   integer(kind=i4) :: nworker
   nworker=0
   end function nworker

   subroutine sharew(w)
   use stack
   type(walker) :: w
   end subroutine sharew

   subroutine sharel(cond)
   logical :: cond
   end subroutine sharel

   subroutine sharei(n)
   integer(kind=i4) :: n
   end subroutine sharei

   subroutine sharer(r)
   real(kind=r8) :: r
   end subroutine sharer

   subroutine getpsi(w,ifrom)
   use stack
   type(walker) :: w
   integer(kind=i4) :: ifrom
   end subroutine getpsi

   subroutine getvc(w,ifrom)
   use stack
   type(walker) :: w
   integer(kind=i4) :: ifrom
   end subroutine getvc

   subroutine getw(w,ifrom)
   use stack
   type(walker) :: w
   integer(kind=i4) :: ifrom
   end subroutine getw

   function amaboss()
   logical :: amaboss
   amaboss=.true.
   end function amaboss

   function amaworker()
   logical :: amaworker
   amaworker=.false.
   end function amaworker

   function idworker()
   integer(kind=i4) :: idworker
   idworker=0
   end function idworker

   function idmyboss()
   integer(kind=i4) :: idmyboss
   idmyboss=-1
   end function idmyboss

end module mympi
