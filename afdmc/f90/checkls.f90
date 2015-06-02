module checkls

   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   integer (kind=i4), private, parameter :: levi(2,3)= &
      reshape((/2,3, 3,1, 1,2/),(/2,3/))

contains


   subroutine chkls(w,dphi,npart)
   use stack
   type (walker) :: w
   real(kind=r8) :: dphi,dx(3),r
   complex(kind=r8) :: opls,oplst
   complex(kind=r8) :: op(0:3,0:3,0:3,0:3,0:3,0:3)
   integer(kind=i4) :: l1,l2,isi,isj,iti,itj,it,ic1,i,j,npart
   do i=2,npart
   do j=1,i-1
   do l1=0,3
!     do l2=0,3
l2=0
         do isi=0,3
            do isj=0,3
               do iti=0,3
                  do itj=0,3
                     op(l1,l2,isi,isj,iti,itj)= &
                        chkopij(w,dphi,i,j,l1,l2,isi,isj,iti,itj)
                  enddo
               enddo
            enddo
         enddo
!     enddo
   enddo
   dx=w%x(:,i)-w%x(:,j)
   r=sqrt(sum(dx**2))
   dx=dx/r
   opls=czero
   oplst=czero
   do ic1=1,3
      opls=opls+0.5_r8*(op(ic1,0,ic1,0,0,0)+op(ic1,0,0,ic1,0,0))
      do it=1,3
         oplst=oplst+0.5_r8*(op(ic1,0,ic1,0,it,it)+op(ic1,0,0,ic1,it,it))
      enddo
   enddo
   write(6,'(''chk , i,j,opls  '',2i3,2e25.15)') i,j,opls
   write(6,'(''chk , i,j,oplst '',2i3,2e25.15)') i,j,oplst
!write (6,*) 'op(ic1,0,ic1,0,it,it)',((op(ic1,0,ic1,0,it,it),ic1=0,3),it=0,3)
   enddo
   enddo
   end subroutine chkls

   function chkopij(w,dphi,i,j,l1,l2,isi,isj,iti,itj)
   use stack
   use wavefunction
   type (walker) :: w
   real(kind=r8) :: dphi,dx(3),cm(3),c,s,xr1(3),xr2(3),r
   integer(kind=i4) :: i,j,isi,isj,l1,l2,iti,itj,ll
   complex(kind=r8) :: chkopij,psi0,sisave(4),sjsave(4),psi,stemp(4)
   complex(kind=r8) :: psipp,psimp,psipm,psimm
call setout(.false.)
   call hpsi(w)
   psi0=w%psi
   dx=w%x(:,i)-w%x(:,j)
   r=sqrt(sum(dx**2))
   dx=0.5_r8*dx
   cm=dx+w%x(:,j)
   sisave=w%sp(:,i)
   sjsave=w%sp(:,j)
   select case (isi)
      case (0)
      case (1)
         w%sp(1,i)=sisave(2)
         w%sp(2,i)=sisave(1)
         w%sp(3,i)=sisave(4)
         w%sp(4,i)=sisave(3)
     case (2)
         w%sp(1,i)=-ci*sisave(2)
         w%sp(2,i)=ci*sisave(1)
         w%sp(3,i)=-ci*sisave(4)
         w%sp(4,i)=ci*sisave(3)
     case (3)
         w%sp(1,i)=sisave(1)
         w%sp(2,i)=-sisave(2)
         w%sp(3,i)=sisave(3)
         w%sp(4,i)=-sisave(4)
    end select
   select case (isj)
      case (0)
      case (1)
         w%sp(1,j)=sjsave(2)
         w%sp(2,j)=sjsave(1)
         w%sp(3,j)=sjsave(4)
         w%sp(4,j)=sjsave(3)
      case (2)
         w%sp(1,j)=-ci*sjsave(2)
         w%sp(2,j)=ci*sjsave(1)
         w%sp(3,j)=-ci*sjsave(4)
         w%sp(4,j)=ci*sjsave(3)
      case (3)
         w%sp(1,j)=sjsave(1)
         w%sp(2,j)=-sjsave(2)
         w%sp(3,j)=sjsave(3)
         w%sp(4,j)=-sjsave(4)
   end select
   stemp=w%sp(:,i)
   select case (iti)
      case (0)
      case (1)
         w%sp(1,i)=stemp(3)
         w%sp(2,i)=stemp(4)
         w%sp(3,i)=stemp(1)
         w%sp(4,i)=stemp(2)
      case (2)
         w%sp(1,i)=-ci*stemp(3)
         w%sp(2,i)=-ci*stemp(4)
         w%sp(3,i)=ci*stemp(1)
         w%sp(4,i)=ci*stemp(2)
      case (3)
         w%sp(1,i)=stemp(1)
         w%sp(2,i)=stemp(2)
         w%sp(3,i)=-stemp(3)
         w%sp(4,i)=-stemp(4)
   end select
   stemp=w%sp(:,j)
   select case (itj)
      case (0)
      case (1)
         w%sp(1,j)=stemp(3)
         w%sp(2,j)=stemp(4)
         w%sp(3,j)=stemp(1)
         w%sp(4,j)=stemp(2)
      case (2)
         w%sp(1,j)=-ci*stemp(3)
         w%sp(2,j)=-ci*stemp(4)
         w%sp(3,j)=ci*stemp(1)
         w%sp(4,j)=ci*stemp(2)
      case (3)
         w%sp(1,j)=stemp(1)
         w%sp(2,j)=stemp(2)
         w%sp(3,j)=-stemp(3)
         w%sp(4,j)=-stemp(4)
   end select
    c=cos(dphi)
    s=sin(dphi)
    call hpsi(w)
    psi=w%psi
    if (l1.eq.0.and.l2.eq.0) then
       chkopij=psi/psi0
    else if (l1.eq.0.or.l2.eq.0) then
       ll=max(l1,l2)
       xr1(levi(1,ll))=c*dx(levi(1,ll))+s*dx(levi(2,ll))
       xr1(levi(2,ll))=c*dx(levi(2,ll))-s*dx(levi(1,ll))
       xr1(ll)=dx(ll)
       w%x(:,i)=cm+xr1
       w%x(:,j)=cm-xr1
       call hpsi(w)
       psipp=w%psi
       xr1(levi(1,ll))=c*dx(levi(1,ll))-s*dx(levi(2,ll))
       xr1(levi(2,ll))=c*dx(levi(2,ll))+s*dx(levi(1,ll))
       xr1(ll)=dx(ll)
       w%x(:,i)=cm+xr1
       w%x(:,j)=cm-xr1
       call hpsi(w)
       psimm=w%psi
       chkopij=(psipp-psimm)/(2.0_r8*dphi*psi0*ci)
    else
      xr1(levi(1,l1))=c*dx(levi(1,l1))+s*dx(levi(2,l1))
      xr1(levi(2,l1))=c*dx(levi(2,l1))-s*dx(levi(1,l1))
      xr1(l1)=dx(l1)
      xr2(levi(1,l2))=c*xr1(levi(1,l2))+s*xr1(levi(2,l2))
      xr2(levi(2,l2))=c*xr1(levi(2,l2))-s*xr1(levi(1,l2))
      xr2(l2)=xr1(l2)
      w%x(:,i)=cm+xr2
      w%x(:,j)=cm-xr2
      call hpsi(w)
      psipp=w%psi
      xr2(levi(1,l2))=c*xr1(levi(1,l2))-s*xr1(levi(2,l2))
      xr2(levi(2,l2))=c*xr1(levi(2,l2))+s*xr1(levi(1,l2))
      xr2(l2)=xr1(l2)
      w%x(:,i)=cm+xr2
      w%x(:,j)=cm-xr2
      call hpsi(w)
      psipm=w%psi
      xr1(levi(1,l1))=c*dx(levi(1,l1))-s*dx(levi(2,l1))
      xr1(levi(2,l1))=c*dx(levi(2,l1))+s*dx(levi(1,l1))
      xr1(l1)=dx(l1)
      xr2(levi(1,l2))=c*xr1(levi(1,l2))+s*xr1(levi(2,l2))
      xr2(levi(2,l2))=c*xr1(levi(2,l2))-s*xr1(levi(1,l2))
      xr2(l2)=xr1(l2)
      w%x(:,i)=cm+xr2
      w%x(:,j)=cm-xr2
      call hpsi(w)
      psimp=w%psi
      xr2(levi(1,l2))=c*xr1(levi(1,l2))-s*xr1(levi(2,l2))
      xr2(levi(2,l2))=c*xr1(levi(2,l2))+s*xr1(levi(1,l2))
      xr2(l2)=xr1(l2)
      w%x(:,i)=cm+xr2
      w%x(:,j)=cm-xr2
      call hpsi(w)
      psimm=w%psi
      chkopij=-(psipp-psipm-psimp+psimm)/(4.0_r8*dphi**2*psi0)
   endif
   w%x(:,i)=cm+dx
   w%x(:,j)=cm-dx
   w%sp(:,i)=sisave
   w%sp(:,j)=sjsave
call setout(.true.)
   end function chkopij

end module checkls

