!$Id: stack.f90,v 1.4 2013/12/10 21:21:53 nuclear Exp $
module stack
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   type :: walker
      real(kind=r8), pointer :: x(:,:)
      complex(kind=r8), pointer :: dpsi(:,:),sp(:,:)
      complex(kind=r8), pointer :: sigma(:,:),tau(:,:),sigmatau(:,:,:)
      complex(kind=r8) :: psi,d2psi,v,weight,psi0,psig,psig0,wt0
      complex(kind=r8) :: vcoul,v8all(8)
      real(kind=r8) :: vc
      integer(kind=i8) :: irn
   end type

!
! define 2 temporary walker variables that can be used by any subroutine
! the save is only to ensure that the allocated memory is not lost
! maybe it is not needed
!
   type (walker), public, save :: w1,w2
   type (walker), public, save :: wtmp,wtmp1,wtmp2,wtmp3,wtmp4

   type (walker), private, allocatable, save :: s(:,:)
   integer(kind=i4), private, allocatable, save :: ist(:) 
   integer(kind=i4), private, save :: nstack,nwalk

   interface assignment (=)
      module procedure copywalker
   end interface

   contains
      subroutine copywalker(wl,wr)
         type (walker), intent(out) :: wl
         type (walker), intent(in) :: wr
         wl%x(:,:)=wr%x(:,:)
         wl%dpsi(:,:)=wr%dpsi(:,:)
         wl%sp(:,:)=wr%sp(:,:)
         wl%sigma(:,:)=wr%sigma(:,:)
         wl%tau(:,:)=wr%tau(:,:)
         wl%sigmatau(:,:,:)=wr%sigmatau(:,:,:)
         wl%psi=wr%psi
         wl%d2psi=wr%d2psi
         wl%v=wr%v
         wl%weight=wr%weight
         wl%wt0=wr%wt0
         wl%vcoul=wr%vcoul
         wl%v8all(:)=wr%v8all(:)
         wl%vc=wr%vc
         wl%psi0=wr%psi0
         wl%psig=wr%psig
         wl%psig0=wr%psig0
         wl%irn=wr%irn
      end subroutine copywalker

      subroutine create(nst,nw,npart)
      integer(kind=i4) :: nst,nw,npart
      integer(kind=i4) :: i,j
      nwalk=nw
      nstack=nst
      allocate(s(nwalk,nstack),ist(nstack))
      do i=1,nstack
         ist(i)=0
         do j=1,nwalk
            allocate(s(j,i)%x(3,npart),s(j,i)%dpsi(3,npart),s(j,i)%sp(4,npart))
            allocate(s(j,i)%sigma(3,npart))
            allocate(s(j,i)%tau(3,npart),s(j,i)%sigmatau(3,3,npart))
         enddo
      enddo
      allocate(w1%x(3,npart),w1%dpsi(3,npart),w1%sp(4,npart),w1%sigma(3,npart))
      allocate(w1%tau(3,npart),w1%sigmatau(3,3,npart))
      allocate(w2%x(3,npart),w2%dpsi(3,npart),w2%sp(4,npart),w2%sigma(3,npart))
      allocate(w2%tau(3,npart),w2%sigmatau(3,3,npart))
      allocate(wtmp%x(3,npart),wtmp%dpsi(3,npart),wtmp%sp(4,npart),wtmp%sigma(3,npart))
      allocate(wtmp%tau(3,npart),wtmp%sigmatau(3,3,npart))
      allocate(wtmp1%x(3,npart),wtmp1%dpsi(3,npart),wtmp1%sp(4,npart),wtmp1%sigma(3,npart))
      allocate(wtmp1%tau(3,npart),wtmp1%sigmatau(3,3,npart))
      allocate(wtmp2%x(3,npart),wtmp2%dpsi(3,npart),wtmp2%sp(4,npart),wtmp2%sigma(3,npart))
      allocate(wtmp2%tau(3,npart),wtmp2%sigmatau(3,3,npart))
      allocate(wtmp3%x(3,npart),wtmp3%dpsi(3,npart),wtmp3%sp(4,npart),wtmp3%sigma(3,npart))
      allocate(wtmp3%tau(3,npart),wtmp3%sigmatau(3,3,npart))
      allocate(wtmp4%x(3,npart),wtmp4%dpsi(3,npart),wtmp4%sp(4,npart),wtmp4%sigma(3,npart))
      allocate(wtmp4%tau(3,npart),wtmp4%sigmatau(3,3,npart))
      return
      end subroutine create

      subroutine push(i,w)
      integer(kind=i4) :: i
      type(walker) :: w
      if (ist(i).ge.nwalk) then
         write (6,'(1x,'' stack overflow'',i10)') ist(i)
         stop
         endif
      ist(i)=ist(i)+1
      s(ist(i),i)=w
      return
      end subroutine push

      subroutine pop(i,w,empty)
      integer(kind=i4) :: i
      type(walker) :: w
      logical :: empty
      empty=ist(i).eq.0
      if (.not.empty) then
         w=s(ist(i),i)
         ist(i)=ist(i)-1
         endif
      return
      end subroutine pop

      function numstack(i)
      integer(kind=i4) :: i,numstack
      numstack=ist(i)
      return
      end function numstack
end module stack
