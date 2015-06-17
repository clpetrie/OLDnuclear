module stack
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   type :: walker
      real(kind=r8), pointer :: x(:,:)
      complex(kind=r8), pointer :: dpsi(:,:),sp(:,:)
      complex(kind=r8) :: psi,d2psi,weight,psi0,psig,psig0,wt0
      complex(kind=r8) :: vcoul,v8all(8),v8allpr(8)
      complex(kind=r8) :: tni2pia,tni2piaxd,tni2piadd
      complex(kind=r8) :: tni2pic,tni2picxd,tni2picdd
      complex(kind=r8) :: tni2pitm,tnivd1,tnivd2,tnive
      complex(kind=r8) :: tni2piapr,tni2piaxdpr,tni2piaddpr
      real(kind=r8) :: vc,vcoulc,tnic,vext
      integer(kind=i8) :: irn
      real(kind=r8), pointer :: x0(:,:)
      complex(kind=r8), pointer :: sp0(:,:),osp(:,:,:,:)
   end type

!
! define 2 temporary walker variables that can be used by any subroutine
! the save is only to ensure that the allocated memory is not lost
! maybe it is not needed
!
   type (walker), public, save :: w1,w2
   type (walker), public, save :: wtmp1,wtmp2,wtmp3,wtmp4

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
         wl%psi=wr%psi
         wl%d2psi=wr%d2psi
         wl%weight=wr%weight
         wl%wt0=wr%wt0
         wl%vcoul=wr%vcoul
         wl%v8all(:)=wr%v8all(:)
         wl%v8allpr(:)=wr%v8allpr(:)
         wl%tni2pia=wr%tni2pia
         wl%tni2piaxd=wr%tni2piaxd
         wl%tni2piadd=wr%tni2piadd
         wl%tni2pic=wr%tni2pic
         wl%tni2picxd=wr%tni2picxd
         wl%tni2picdd=wr%tni2picdd
         wl%tni2pitm=wr%tni2pitm
         wl%tnivd1=wr%tnivd1
         wl%tnivd2=wr%tnivd2
         wl%tnive=wr%tnive
         wl%tni2piapr=wr%tni2piapr
         wl%tni2piaxdpr=wr%tni2piaxdpr
         wl%tni2piaddpr=wr%tni2piaddpr
         wl%vc=wr%vc
         wl%vcoulc=wr%vcoulc
         wl%tnic=wr%tnic
         wl%vext=wr%vext
         wl%psi0=wr%psi0
         wl%psig=wr%psig
         wl%psig0=wr%psig0
         wl%irn=wr%irn
         wl%x0(:,:)=wr%x0(:,:)
         wl%sp0(:,:)=wr%sp0(:,:)
         wl%osp(:,:,:,:)=wr%osp(:,:,:,:)
      end subroutine copywalker

      subroutine create(nst,nw,npart,nop,nspop)
      integer(kind=i4) :: nst,nw,npart,nop,nspop
      integer(kind=i4) :: i,j
      nwalk=nw
      nstack=nst
      allocate(s(nwalk,nstack),ist(nstack))
      do i=1,nstack
         ist(i)=0
         do j=1,nwalk
            allocate(s(j,i)%x(3,npart),s(j,i)%dpsi(3,npart),s(j,i)%sp(4,npart))
            if (nop.eq.0) then
               allocate(s(j,i)%x0(1,1),s(j,i)%sp0(1,1),s(j,i)%osp(1,1,1,1))
            else
               allocate(s(j,i)%x0(3,npart),s(j,i)%sp0(4,npart),s(j,i)%osp(nop,4,nspop,npart))
            endif
         enddo
      enddo
      allocate(w1%x(3,npart),w1%dpsi(3,npart),w1%sp(4,npart))
      allocate(w2%x(3,npart),w2%dpsi(3,npart),w2%sp(4,npart))
      allocate(wtmp1%x(3,npart),wtmp1%dpsi(3,npart),wtmp1%sp(4,npart))
      allocate(wtmp2%x(3,npart),wtmp2%dpsi(3,npart),wtmp2%sp(4,npart))
      allocate(wtmp3%x(3,npart),wtmp3%dpsi(3,npart),wtmp3%sp(4,npart))
      allocate(wtmp4%x(3,npart),wtmp4%dpsi(3,npart),wtmp4%sp(4,npart))
      if (nop.eq.0) then
         allocate(w1%x0(1,1),w1%sp0(1,1),w1%osp(1,1,1,1))
         allocate(w2%x0(1,1),w2%sp0(1,1),w2%osp(1,1,1,1))
         allocate(wtmp1%x0(1,1),wtmp1%sp0(1,1),wtmp1%osp(1,1,1,1))
         allocate(wtmp2%x0(1,1),wtmp2%sp0(1,1),wtmp2%osp(1,1,1,1))
         allocate(wtmp3%x0(1,1),wtmp3%sp0(1,1),wtmp3%osp(1,1,1,1))
         allocate(wtmp4%x0(1,1),wtmp4%sp0(1,1),wtmp4%osp(1,1,1,1))
      else
         allocate(w1%x0(3,npart),w1%sp0(4,npart),w1%osp(nop,4,nspop,npart))
         allocate(w2%x0(3,npart),w2%sp0(4,npart),w2%osp(nop,4,nspop,npart))
         allocate(wtmp1%x0(3,npart),wtmp1%sp0(4,npart),wtmp1%osp(nop,4,nspop,npart))
         allocate(wtmp2%x0(3,npart),wtmp2%sp0(4,npart),wtmp2%osp(nop,4,nspop,npart))
         allocate(wtmp3%x0(3,npart),wtmp3%sp0(4,npart),wtmp3%osp(nop,4,nspop,npart))
         allocate(wtmp4%x0(3,npart),wtmp4%sp0(4,npart),wtmp4%osp(nop,4,nspop,npart))
      endif
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
