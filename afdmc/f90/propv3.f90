module propv3
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, parameter :: levi(2,3)=reshape((/2,3, 3,1, 1,2/),(/2,3/))
   complex(kind=r8), private, parameter :: czero = (0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci = (0.0_r8,1.0_r8)
   integer(kind=i4), private, save :: npart,ntrip
   integer(kind=i4), private, save, allocatable :: ijktab(:,:)
   real(kind=r8), private, save :: dt,dtt
contains
   subroutine initgv3(npartin,dtin)
   use mympi
   integer(kind=i4) :: npartin
   real(kind=r8) :: dtin
   integer(kind=i4) :: i,j,k,ijk
   npart=npartin
   ntrip=npart*(npart-1)*(npart-2)/6
   write(6,'(''number of triples ='',t40,i10)') ntrip
   allocate(ijktab(3,ntrip))
   ijktab=0
   ijk=0
   do i=1,npart-2
      do j=i+1,npart-1
         do k=j+1,npart
            ijk=ijk+1
            ijktab(1,ijk)=i
            ijktab(2,ijk)=j
            ijktab(3,ijk)=k
         enddo
      enddo
   enddo
!  do i=1,ntrip
!     write(6,*) i,ijktab(1,i),ijktab(2,i),ijktab(3,i)
!  enddo
   dt=dtin
   dtt=dt**(1.0_r8/3.0_r8)
   end subroutine initgv3
            
   subroutine gv3(w)
   use stack
   use v3bpot
   use random
   type(walker) :: w
   real(kind=r8) :: vc,xpi(3,npart,3,npart),g2s3b(3,npart,npart),delta(npart,npart)
   real(kind=r8) :: dummy,a2pscalin,a2pxdscalin,a2pddscalin
   real(kind=r8) :: xgauss(3)
   complex(kind=r8) :: sp(4,3),spnew(4,3),spn(8,4,3)
   integer(kind=i4) :: ijk,i,j,k
   real(kind=r8) :: wttot,wt(8),wtt,wttt,rn(1),wtfacc,wtc
   real(kind=r8) :: xx(3**3*6),yy(3**3*6),zz(3**3*6)
   complex(kind=r8) :: rat
   call hstnipot(w%x,vc,xpi,g2s3b,delta,dummy,dummy,dummy,dummy,dummy,a2pscalin,a2pxdscalin,a2pddscalin,dummy,1.0_r8)
   wttot=1.0_r8
   do ijk=1,ntrip
      i=ijktab(1,ijk)
      j=ijktab(2,ijk)
      k=ijktab(3,ijk)
      sp(:,1)=w%sp(:,i)
      sp(:,2)=w%sp(:,j)
      sp(:,3)=w%sp(:,k)
      xx=gaussian(3**3*6)
      yy=gaussian(3**3*6)
      zz=gaussian(3**3*6)
      call g3op(xx,yy,zz,i,j,k,xpi,sp,spn(1,:,:))
      call wfrat(sp,spn,i,j,k,rat,wtfacc)
      wt(1)=exp(-wtfacc*dt)*real(rat)
      call g3op(xx,yy,-zz,i,j,k,xpi,sp,spn(2,:,:))
      call wfrat(sp,spn,i,j,k,rat,wtfacc)
      wt(2)=exp(-wtfacc*dt)*real(rat)
      call g3op(xx,-yy,zz,i,j,k,xpi,sp,spn(3,:,:))
      call wfrat(sp,spn,i,j,k,rat,wtfacc)
      wt(3)=exp(-wtfacc*dt)*real(rat)
      call g3op(xx,-yy,-zz,i,j,k,xpi,sp,spn(4,:,:))
      call wfrat(sp,spn,i,j,k,rat,wtfacc)
      wt(4)=exp(-wtfacc*dt)*real(rat)
      call g3op(-xx,yy,zz,i,j,k,xpi,sp,spn(5,:,:))
      call wfrat(sp,spn,i,j,k,rat,wtfacc)
      wt(5)=exp(-wtfacc*dt)*real(rat)
      call g3op(-xx,yy,-zz,i,j,k,xpi,sp,spn(6,:,:))
      call wfrat(sp,spn,i,j,k,rat,wtfacc)
      wt(6)=exp(-wtfacc*dt)*real(rat)
      call g3op(-xx,-yy,zz,i,j,k,xpi,sp,spn(7,:,:))
      call wfrat(sp,spn,i,j,k,rat,wtfacc)
      wt(7)=exp(-wtfacc*dt)*real(rat)
      call g3op(-xx,-yy,-zz,i,j,k,xpi,sp,spn(8,:,:))
      call wfrat(sp,spn,i,j,k,rat,wtfacc)
      wt(8)=exp(-wtfacc*dt)*real(rat)
      wtt=sum(wt(:))
      wt=max(0.0_r8,wt)
      wttt=sum(wt)
      wt=wt/wttt
      do i=2,8
         wt(i)=wt(i)+wt(i-1)
      enddo
      rn=randn(1)
      if (rn(1).lt.wt(1)) then
         spnew=spn(1,:,:)
      else if (rn(1).lt.wt(2)) then
         spnew=spn(2,:,:)
      else if (rn(1).lt.wt(3)) then
         spnew=spn(3,:,:)
      else if (rn(1).lt.wt(4)) then
         spnew=spn(4,:,:)
      else if (rn(1).lt.wt(5)) then
         spnew=spn(5,:,:)
      else if (rn(1).lt.wt(6)) then
         spnew=spn(6,:,:)
      else if (rn(1).lt.wt(7)) then
         spnew=spn(7,:,:)
      else if (rn(1).lt.wt(8)) then
         spnew=spn(8,:,:)
      else 
         write(6,'(''Error in gv3!!!'')')
      endif
      wtc=wtt/8.0_r8
      w%sp(:,i)=spnew(:,1)
      w%sp(:,j)=spnew(:,2)
      w%sp(:,k)=spnew(:,3)
      wttot=wttot*wtc
   enddo
   end subroutine gv3

   subroutine g3op(xx,yy,zz,i,j,k,xpi,spold,spnew)
   real(kind=r8) :: xx(:),yy(:),zz(:)
   complex(kind=r8) :: spold(:,:),spnew(:,:) ! spins of the triple
   integer(kind=i4) :: ic,jc,kc,it,i,j,k,idx
   real(kind=r8) :: opi,opj,opk
   real(kind=r8) :: sti(3,3),stj(3,3),stk(3,3)
   real(kind=r8) :: xpi(3,npart,3,npart)
   sti=0.0_r8
   stj=0.0_r8
   stk=0.0_r8
   idx=1
   do ic=1,3
      do jc=1,3
         do kc=1,3
            opi=-xx(idx)*yy(idx)*dtt*xpi(ic,i,levi(1,jc),j)
            opj=-xx(idx)*zz(idx)*dtt
            opk=-yy(idx)*zz(idx)*dtt*xpi(levi(2,jc),j,kc,k)
            idx=idx+1
            do it=1,3
               sti(ic,levi(1,it))=sti(ic,levi(1,it))+opi
               stj(jc,it)=stj(jc,it)+opj
               stk(kc,levi(2,it))=stk(kc,levi(2,it))+opk
            enddo
            opi=xx(idx)*yy(idx)*dtt*xpi(ic,i,levi(2,jc),j)
            opj=-xx(idx)*zz(idx)*dtt
            opk=yy(idx)*zz(idx)*dtt*xpi(levi(1,jc),j,kc,k)
            idx=idx+1
            do it=1,3
               sti(ic,levi(2,it))=sti(ic,levi(2,it))-opi
               stj(jc,it)=stj(jc,it)+opj
               stk(kc,levi(1,it))=stk(kc,levi(1,it))-opk
            enddo
            opi=-xx(idx)*yy(idx)*dtt*xpi(jc,j,levi(1,kc),k)
            opj=-xx(idx)*zz(idx)*dtt
            opk=-yy(idx)*zz(idx)*dtt*xpi(levi(2,kc),k,ic,i)
            idx=idx+1
            do it=1,3
               sti(ic,levi(1,it))=sti(ic,levi(1,it))+opi
               stj(jc,it)=stj(jc,it)+opj
               stk(kc,levi(2,it))=stk(kc,levi(2,it))+opk
            enddo
            opi=xx(idx)*yy(idx)*dtt*xpi(jc,j,levi(2,kc),k)
            opj=-xx(idx)*zz(idx)*dtt
            opk=yy(idx)*zz(idx)*dtt*xpi(levi(1,kc),k,ic,i)
            idx=idx+1
            do it=1,3
               sti(ic,levi(2,it))=sti(ic,levi(2,it))-opi
               stj(jc,it)=stj(jc,it)+opj
               stk(kc,levi(1,it))=stk(kc,levi(1,it))-opk
            enddo
            opi=-xx(idx)*yy(idx)*dtt*xpi(kc,k,levi(1,ic),i)
            opj=-xx(idx)*zz(idx)*dtt
            opk=-yy(idx)*zz(idx)*dtt*xpi(levi(2,ic),i,jc,j)
            idx=idx+1
            do it=1,3
               sti(ic,levi(1,it))=sti(ic,levi(1,it))+opi
               stj(jc,it)=stj(jc,it)+opj
               stk(kc,levi(2,it))=stk(kc,levi(2,it))+opk
            enddo
            opi=xx(idx)*yy(idx)*dtt*xpi(kc,k,levi(2,ic),i)
            opj=-xx(idx)*zz(idx)*dtt
            opk=yy(idx)*zz(idx)*dtt*xpi(levi(1,ic),i,jc,j)
            idx=idx+1
            do it=1,3
               sti(ic,levi(2,it))=sti(ic,levi(2,it))-opi
               stj(jc,it)=stj(jc,it)+opj
               stk(kc,levi(1,it))=stk(kc,levi(1,it))-opk
            enddo
         enddo
      enddo
   enddo
   call rotate(spold(:,1),spnew(:,1),sti)
   call rotate(spold(:,2),spnew(:,2),stj)
   call rotate(spold(:,3),spnew(:,3),stk)
   end subroutine g3op

   subroutine rotate(spold,spnew,rotst)
   use matrixmod
   complex(kind=r8) :: spold(4),spnew(4)
   complex(kind=r8) :: stmat(4,4)
   real(kind=r8) :: rotst(3,3)
   stmat=0.0_r8
   stmat(1,1)=rotst(3,3)
   stmat(1,2)=rotst(1,3)-ci*rotst(2,3)
   stmat(1,3)=rotst(3,1)-ci*rotst(3,2)
   stmat(1,4)=rotst(1,1)-ci*rotst(2,1)-ci*rotst(1,2)-rotst(2,2)
   stmat(2,1)=rotst(1,3)+ci*rotst(2,3)
   stmat(2,2)=-rotst(3,3)
   stmat(2,3)=rotst(1,1)-ci*rotst(1,2)+ci*rotst(2,1)+rotst(2,2)
   stmat(2,4)=-rotst(3,1)+ci*rotst(3,2)
   stmat(3,1)=rotst(3,1)+ci*rotst(3,2)
   stmat(3,2)=rotst(1,1)+ci*rotst(1,2)-ci*rotst(2,1)+rotst(2,2)
   stmat(3,3)=-rotst(3,3)
   stmat(3,4)=-rotst(1,3)+ci*rotst(2,3)
   stmat(4,1)=rotst(1,1)+ci*rotst(2,1)+ci*rotst(1,2)-rotst(2,2)
   stmat(4,2)=-rotst(3,1)-ci*rotst(3,2)
   stmat(4,3)=-rotst(1,3)-ci*rotst(2,3)
   stmat(4,4)=rotst(3,3)
   spnew(:)=expmult(stmat(:,:),spold(:),4)
   end subroutine rotate

   subroutine wfrat(sp,spn,i,j,k,rat,wtfacc)
   complex(kind=r8) :: sp(:,:),spn(:,:,:),rat
   real(kind=r8) :: wtfacc
   integer(kind=i4) :: i,j,k
! fake, just to compile
   end subroutine
 
end module propv3
